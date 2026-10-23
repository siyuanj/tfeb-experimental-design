#!/usr/bin/env python3
"""Check the final handoff chain without rewriting scientific outputs."""

import hashlib
import json
import re
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "docs/协作/主线验收核查.json"
GATES = [f"docs/协作/{name}-ready.json" for name in ("sequence", "methods", "report")]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    OUTPUT.write_text(json.dumps({"status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat()}) + "\n")
    snapshots = {}
    checks = {}
    gate_counts = {}

    def read(relative):
        path = (ROOT / relative).resolve(strict=True)
        if not path.is_relative_to(ROOT):
            raise ValueError(f"Input escapes repository: {relative}")
        data = path.read_bytes()
        key = str(path.relative_to(ROOT))
        current = digest(data)
        if key in snapshots and current != snapshots[key]:
            raise ValueError(f"Input changed during audit: {key}")
        snapshots[key] = current
        return data

    def record(name, passed):
        checks[name] = bool(passed)
        if not passed:
            raise ValueError(f"Check failed: {name}")

    for relative in GATES:
        gate = json.loads(read(relative))
        record(f"{relative}:ready", gate.get("status") == "ready")
        files = gate.get("files")
        record(f"{relative}:nonempty_bindings", isinstance(files, dict) and bool(files))
        for filename, expected in files.items():
            record(f"sha256:{filename}", digest(read(filename)) == expected)
        for filename, expected in gate.get("inputs", {}).items():
            record(f"input_sha256:{filename}", digest(read(filename)) == expected)
        gate_counts[relative] = len(files)

    primary = json.loads(read("设计/verification_summary.json"))
    independent = json.loads(read("设计/independent_verification.json"))
    record("primary_scientific_checks", primary.get("status") == "passed" and all(primary["checks"].values()))
    record("independent_scientific_checks", independent.get("status") == "pass")
    record("unchanged_ElegantPaper", digest(read("报告/elegantpaper.cls")) == "29370455a0cdf0af94039a6483ba1433f5b067666771b0f96c4cca0159830173")
    record("project_rules_equal", read("AGENTS.md") == read("CLAUDE.md"))

    tex_paths = [ROOT / "报告/main.tex", *sorted((ROOT / "报告/sections").glob("*.tex"))]
    tex = "\n".join(read(str(p.relative_to(ROOT))).decode() for p in tex_paths)
    bib = read("报告/references.bib").decode()
    keys = re.findall(r"@\w+\s*\{\s*([^,\s]+)\s*,", bib)
    record("unique_bibliography_keys", all(n == 1 for n in Counter(keys).values()))
    cited = set()
    for group in re.findall(r"\\(?:[A-Za-z]*cite[A-Za-z]*)(?:\[[^\]]*\])*\{([^}]+)\}", tex):
        cited.update(k.strip() for k in group.split(","))
    record("all_citation_keys_exist", bool(cited) and cited.issubset(set(keys)))
    includes = re.findall(r"\\input\{([^}]+)\}", tex)
    for name in includes:
        path = ROOT / "报告" / name
        if not path.suffix:
            path = path.with_suffix(".tex")
        read(str(path.relative_to(ROOT)))
    graphics = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
    for name in graphics:
        path = ROOT / "报告" / name
        if not path.suffix:
            matches = [path.with_suffix(ext) for ext in (".pdf", ".png", ".jpg") if path.with_suffix(ext).is_file()]
            if not matches:
                raise ValueError(f"Missing figure: {name}")
            path = matches[0]
        read(str(path.relative_to(ROOT)))
    record("core_figures_in_report", len(graphics) >= 3)

    for name in ("methods-cloning.tex", "results-cloning.tex"):
        text = read(f"报告/sections/{name}").decode()
        record(f"no_substantive_placeholder:{name}", not re.search(r"\b(?:TBD|TODO|FIXME|placeholder)\b", text, re.I))
    record("working_draft_identity_retained", "Working draft" in tex)

    pdf = "报告/main.pdf"
    read(pdf)
    info = subprocess.run(["pdfinfo", str(ROOT / pdf)], check=True, capture_output=True, text=True).stdout
    match = re.search(r"^Pages:\s+(\d+)", info, re.M)
    record("nonempty_readable_PDF", match is not None and int(match[1]) > 0)
    extracted = subprocess.run(["pdftotext", "-layout", str(ROOT / pdf), "-"], check=True, capture_output=True, text=True).stdout
    record("PDF_has_extractable_content", len(extracted.strip()) > 1000)
    record("PDF_no_undefined_citation_marker", "[?]" not in extracted and "(??)" not in extracted)
    record("PDF_no_replacement_character", "\ufffd" not in extracted)

    # A second read rejects concurrent edits rather than binding a mixed snapshot.
    for filename, expected in list(snapshots.items()):
        if digest((ROOT / filename).read_bytes()) != expected:
            raise ValueError(f"Input changed before audit publication: {filename}")
    payload = {
        "status": "passed",
        "scope": "Read-only final structural/hash/reference audit; manual science and page review are recorded separately in docs/B验收记录.md.",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "gate_file_counts": gate_counts,
        "checks": checks,
        "citation_keys_used": sorted(cited),
        "bibliography_entries": len(keys),
        "graphics_includes": graphics,
        "pdf_pages": int(match[1]),
        "files": snapshots,
        "audit_script_sha256": digest(Path(__file__).read_bytes()),
        "limitations": ["No wet-lab data", "Does not validate independent A protocols", "Text extraction does not replace visual PDF inspection", "Does not assert Git/PR delivery"],
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": payload["status"], "checks": len(checks), "bound_files": len(snapshots), "pdf_pages": payload["pdf_pages"]}))


if __name__ == "__main__":
    main()
