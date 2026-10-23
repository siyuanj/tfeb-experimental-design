#!/usr/bin/env python3
"""Verify local source-file hashes and document the S14 history correction."""

import csv
import hashlib
import io
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "docs/协作/主线来源核验.json"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def audit():
    ledger_path = "资料/B来源索引.tsv"
    ledger = (ROOT / ledger_path).read_bytes()
    snapshots = {ledger_path: sha(ledger)}
    records = []
    rows = list(csv.DictReader(io.StringIO(ledger.decode()), delimiter="\t"))
    assert len({row['source_id'] for row in rows}) == len(rows)
    for row in rows:
        expected = row["raw_file_sha256"]
        if not re.fullmatch(r"[0-9a-f]{64}", expected):
            assert expected == "not applicable", row['source_id']
            continue
        relative = row["raw_file_path"]
        assert not Path(relative).is_absolute(), relative
        path = (ROOT / relative).resolve(strict=True)
        assert path.is_relative_to(ROOT), relative
        actual = sha(path.read_bytes())
        assert actual == expected, f"{row['source_id']}: {relative}"
        snapshots[relative] = actual
        records.append({"source_id": row['source_id'], "path": relative, "sha256": actual})

    course = "docs/10-23 课程要求与待办.md"
    original = subprocess.check_output(["git", "show", f"8bdd898:{course}"], cwd=ROOT)
    assert sha(original) == "d3800629bfe0f3a80f55f1e4e86be9e9b5915e7d5e2f275ecb416690cdd78843"
    current = (ROOT / course).read_bytes()
    assert original.split("## 已确定".encode())[0] == current.split("## 已确定".encode())[0]
    gates = {}
    for name in ("sequence", "methods", "report"):
        relative = f"docs/协作/{name}-ready.json"
        data = (ROOT / relative).read_bytes()
        gate = json.loads(data)
        assert ledger_path not in {**gate.get("inputs", {}), **gate.get("files", {})}
        gates[relative] = sha(data)
        snapshots[relative] = sha(data)
    for relative, expected in snapshots.items():
        assert sha((ROOT / relative).read_bytes()) == expected, f"Changed during review: {relative}"
    return {
        "status": "passed", "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Local files explicitly hashed in the source ledger; no URL-only source re-fetch or scientific regeneration",
        "ledger": {"path": ledger_path, "sha256": sha(ledger), "rows": len(rows)},
        "local_sources_checked": records, "url_only_sources_not_refetched": len(rows) - len(records),
        "S002": {"status": "closed", "source_id": "S14", "historical_commit": "8bdd898",
                 "historical_sha256": sha(original), "current_sha256": sha(current),
                 "requirements_summary_unchanged": True,
                 "change_scope": "Work scope, final construct selection and five B checklist items; requirement-source section unchanged"},
        "unchanged_ready_files": gates,
        "files": snapshots, "script_sha256": sha(Path(__file__).read_bytes()),
        "limitations": ["S05/S06 full sources remain local and ignored; not redistributed", "Does not certify wet experiments, final course submission or all external URLs"]}


if __name__ == "__main__":
    if OUTPUT.is_symlink():
        raise ValueError("Output must be an ordinary project file")
    try:
        result = audit()
    except Exception as error:
        OUTPUT.write_text(json.dumps({"status": "failed", "error": str(error)}, ensure_ascii=False, indent=2) + "\n")
        raise
    OUTPUT.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": result['status'], "source_rows": result['ledger']['rows'],
                      "local_sources_checked": len(result['local_sources_checked']), "S002": "closed"}))
