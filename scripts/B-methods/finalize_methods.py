#!/usr/bin/env python3
"""Publish the P6 ready manifest only for reviewed, hash-bound B-v1 outputs."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

ROOT = Path(__file__).resolve().parents[2]
READY = ROOT / "docs/协作/methods-ready.json"
TEMP = Path(__file__).parent / ".methods-ready.json.tmp"


def sha(path: Path) -> str:
    """Hash exact bytes."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(relative: str):
    return json.loads((ROOT / relative).read_text())


def verify(mapping: dict):
    """Fail closed if an artifact was edited after its recorded check."""
    for name, expected in mapping.items():
        if sha(ROOT / name) != expected:
            raise RuntimeError(f"Changed artifact: {name}")


def publish(value: dict):
    """Replace the marker atomically on the same filesystem."""
    try:
        TEMP.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")
        TEMP.replace(READY)
    finally:
        if TEMP.exists():
            TEMP.unlink()


def main():
    if READY.exists():
        publish({"status": "revalidation_required",
                 "previous_ready_sha256": sha(READY),
                 "reason": "P6 finalization invoked; readiness requires fresh hash checks"})
    primary = load("docs/协作/方法独立核验.json")
    cross = load("docs/协作/方法交叉核验.json")
    parameters = load("设计/methods_parameters.json")
    if primary["status"] != "passed" or len(primary["checks"]) != 28 or not all(primary["checks"].values()):
        raise RuntimeError("The 28 P6 generator checks did not all pass")
    if cross["status"] != "passed" or cross["checks"] != 48 or cross["passed"] != 48 or not all(cross["check_results"].values()):
        raise RuntimeError("The 48 independent checks did not all pass")
    if primary["outputs"] != cross["method_hashes"]:
        raise RuntimeError("The checks did not evaluate the same method outputs")
    if primary["inputs"] != parameters["inputs"]:
        raise RuntimeError("The method parameters refer to another input snapshot")
    for mapping in (primary["inputs"], primary["outputs"], primary["scripts"], cross["input_hashes"]):
        verify(mapping)
    if sha(ROOT / "scripts/B-methods/independent_review.py") != cross["script_sha256"]:
        raise RuntimeError("Independent checker changed after its recorded run")
    if load("docs/协作/sequence-ready.json")["status"] != "ready":
        raise RuntimeError("P5 is no longer ready")
    with (ROOT / "docs/协作/方法来源.tsv").open() as stream:
        sources = list(csv.DictReader(stream, delimiter="\t"))
    ids = [source["source_id"] for source in sources]
    if len(ids) != len(set(ids)) or set(ids) != {f"M{i:02d}" for i in range(1, 21)}:
        raise RuntimeError("Missing or duplicate method sources")
    if set(parameters["source_ids"]) - set(ids):
        raise RuntimeError("Unresolved source keys")
    artifacts = list(primary["outputs"]) + [
        "docs/协作/方法正文.md", "docs/协作/方法来源.tsv",
        "docs/协作/方法独立核验.json", "docs/协作/方法交叉核验.json",
        "docs/协作/方法验收.md", "docs/协作/方法进度.md", "docs/协作/方法问题.md",
    ] + [str(path.relative_to(ROOT)) for path in sorted(Path(__file__).parent.glob("*.py"))]
    file_hashes = {name: sha(ROOT / name) for name in artifacts}
    interpreter = "'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' -B"
    manifest = {
        "status": "ready", "stage": "P6", "design_version": primary["design_version"],
        "ready_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "scope": "Method design and computational checks only; not wet-lab completion or overall B acceptance",
        "inputs": primary["inputs"], "files": file_hashes,
        "sources": {"path": "docs/协作/方法来源.tsv", "source_count": len(sources),
                    "historical_only": ["M03"], "used": parameters["source_ids"]},
        "validation": {"generator_checks_passed": len(primary["checks"]),
                       "independent_checks_passed": cross["passed"],
                       "primary_record": "docs/协作/方法独立核验.json",
                       "independent_record": "docs/协作/方法交叉核验.json",
                       "review_record": "docs/协作/方法验收.md",
                       "software": primary["software"]},
        "metrics": primary["metrics"],
        "working_directory": str(ROOT),
        "commands": [
            f"{interpreter} scripts/B-methods/build_methods.py",
            f"{interpreter} scripts/B-methods/write_methods_text.py",
            f"{interpreter} scripts/B-methods/independent_review.py > docs/协作/方法交叉核验.json",
            f"{interpreter} scripts/B-methods/finalize_methods.py",
        ],
        "unresolved_P4_P5_issues": [], "user_decisions_required": [],
        "limitations": parameters["limitations"] + [
            "Template/vector procurement, stock concentrations and actual reagent inventory have not been verified",
            "PCR temperatures, digestion recovery and ligation/transformation efficiency are proposed, not optimized measurements",
            "Sanger service submission requirements must be checked again before actual ordering",
        ],
    }
    # Recheck inputs and all published files after assembling the manifest.
    verify(manifest["inputs"])
    verify(file_hashes)
    publish(manifest)
    stored = load("docs/协作/methods-ready.json")
    try:
        verify(stored["inputs"])
        verify(stored["files"])
    except Exception:
        publish({"status": "revalidation_required",
                 "reason": "An artifact changed during the post-publication hash check"})
        raise
    print(json.dumps({"status": stored["status"], "stage": stored["stage"],
                      "files": len(stored["files"]), "inputs": len(stored["inputs"]),
                      "ready_sha256": sha(READY)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
