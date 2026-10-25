#!/usr/bin/env python3
"""Bind the reviewed report after compilation, consistency and page inspection.

This script verifies existing review records; it cannot perform or assert a new
manual visual review. A changed PDF requires a new per-page review record.
"""

import json
import shlex
import sys
from datetime import datetime
from zoneinfo import ZoneInfo

from handoff import ROOT, assert_unchanged, sha256, verify_gate


def main():
    ready_path = ROOT / "docs/协作/report-ready.json"
    ready_path.write_text(json.dumps({
        "status": "validating", "stage": "P7/P8",
        "started_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
    }) + "\n")
    _, sequence = verify_gate("docs/协作/sequence-ready.json")
    _, methods = verify_gate("docs/协作/methods-ready.json")
    audit = json.loads((ROOT / "docs/协作/报告一致性核验.json").read_text())
    review = json.loads((ROOT / "docs/协作/报告逐页核验.json").read_text())
    build_path = ROOT / "docs/协作/报告编译核验.json"
    build = json.loads(build_path.read_text())
    record_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in (
        ROOT / "docs/协作/报告一致性核验.json", ROOT / "docs/协作/报告逐页核验.json", build_path,
    )}
    pdf_hash = sha256(ROOT / "报告/main.pdf")
    assert audit["status"] == "passed" and all(audit["checks"].values())
    assert build["status"] == "passed"
    assert sha256(build_path) == audit["compilation_record"]["sha256"]
    assert sha256(build_path) == review["compilation_record_sha256"]
    assert build["pdf"]["sha256"] == pdf_hash
    assert build["sequence_gate"]["sha256"] == sequence["sha256"]
    assert build["methods_gate"]["sha256"] == methods["sha256"]
    assert build["script_sha256"] == sha256(ROOT / "scripts/B-figures/compile_report.py")
    for relative, expected in build["inputs"].items():
        assert sha256(ROOT / relative) == expected, f"Changed compilation input: {relative}"
    for item in (build["pdf_text"], build["compiler_log"], *build["page_images"]):
        assert sha256(ROOT / item["path"]) == item["sha256"], item["path"]
    assert audit["pdf"]["sha256"] == review["reviewed_pdf_sha256"] == pdf_hash
    assert audit["script_sha256"] == sha256(ROOT / "scripts/B-figures/check_report_consistency.py")
    assert audit["sequence_gate"]["sha256"] == sequence["sha256"]
    assert audit["methods_gate"]["sha256"] == methods["sha256"]
    assert review["status"] == "passed" and review["page_count"] == audit["pdf"]["pages"]
    assert [item["page"] for item in review["pages"]] == list(range(1, review["page_count"] + 1))
    for page in review["pages"]:
        assert page["status"] == "visually_reviewed"
        assert sha256(ROOT / page["image"]) == page["sha256"]
        compiled_page = build["page_images"][page["page"] - 1]
        assert compiled_page["page"] == page["page"]
        assert compiled_page["path"] == page["image"]
        assert compiled_page["sha256"] == page["sha256"]
    assert sha256(ROOT / "报告/elegantpaper.cls") == "29370455a0cdf0af94039a6483ba1433f5b067666771b0f96c4cca0159830173"

    paths = [ROOT / "报告/main.tex", ROOT / "报告/main.pdf", ROOT / "报告/references.bib",
             ROOT / "报告/elegantpaper.cls", ROOT / "docs/B交接给A.md"]
    paths += list((ROOT / "报告/sections").glob("*.tex"))
    paths += [p for p in (ROOT / "报告/images/B").rglob('*') if p.suffix in (".pdf", ".svg", ".png", ".json")]
    paths += list((ROOT / "scripts/B-figures").glob("*.py"))
    paths += [ROOT / "scripts/B-figures/README.md", ROOT/'scripts/B-figures/requirements-publication.txt',
              ROOT/'scripts/B-audit/check_control_panel.py',ROOT/'scripts/B-audit/check_research_panel.py',
              ROOT/'scripts/B-audit/check_research_condition_plan.py',
              ROOT/'docs/10-25 审稿视角复查与修订.md']
    paths += [p for p in (ROOT/'设计/controls-v2').iterdir() if p.is_file()]
    paths += [p for p in (ROOT/'设计/research-v3').iterdir() if p.is_file()]
    paths += [p for p in (ROOT/'设计/linker-v4').iterdir() if p.is_file()]
    paths += [ROOT/'scripts/B-audit/check_linker_controls.py', ROOT/'docs/10-25 Claude独立审稿复查与措辞修订.md']
    paths += [p for p in (ROOT/'报告/质粒图-SnapGene').iterdir() if p.is_file()]
    paths += [p for p in (ROOT / "docs/协作").glob("报告*") if p.is_file()]
    paths += [ROOT/'docs/协作/流程图ImageGen提示词-v3.json']
    files = {str(p.relative_to(ROOT)): sha256(p) for p in sorted(set(paths))}
    assert_unchanged(sequence)
    assert_unchanged(methods)
    payload = {
        "status": "ready", "stage": "research-v3 report", "design_version": "B-v1 plus controls-v2 and research-v3",
        "ready_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "scope": "Controlled replication plus C-terminal research extension, ten native SnapGene maps and ten main figures; historical P9 acceptance applies only to B-v1",
        "inputs": {sequence["path"]: sequence["sha256"], methods["path"]: methods["sha256"]},
        "files": files,
        "validation": {
            "compilation_record": "docs/协作/报告编译核验.json",
            "consistency_record": "docs/协作/报告一致性核验.json",
            "visual_record": "docs/协作/报告逐页核验.json",
            "review_record": "docs/协作/报告验收.md",
            "pdf_sha256": pdf_hash, "pdf_pages": review["page_count"],
            "all_pages_visually_reviewed": True,
            "compiler": "Tectonic 0.17.0", "pdf_tools": "Poppler 26.10.0",
            "actual_finalize_command": shlex.join([sys.executable, str(__file__)]),
        },
        "necessary_A_changes": ["Construct and empty-vector identity", "Calculated fusion mass and WB interpretation", "Matched nutrient/vehicle controls", "Localization interpretation and limits", "C-terminal genotype-tag-stimulus comparison and conditional follow-ups", "Honest planned contributions and AI declaration"],
        "remaining_independent_A": ["Exact siRNA targets/qPCR oligos and pilot-validated reagent/dose settings", "Authors' actual contributions, identity and submission confirmation"],
        "unresolved_report_issues": [], "user_decisions_required_for_B_report": [],
        "limitations": ["No wet-lab data or measured sequencing coverage", "Not the final complete course submission", "No author review or individual contribution is certified", "Historical B-v1 P9 does not constitute independent acceptance of this expansion; new work is locally reviewed", "Ignored build intermediates and page rasters remain local evidence; regenerate for a new PDF"],
    }
    for relative, expected in files.items():
        assert sha256(ROOT / relative) == expected, relative
    for relative, expected in {**build["inputs"], **record_hashes}.items():
        assert sha256(ROOT / relative) == expected, f"Input changed during finalization: {relative}"
    assert sha256(ROOT / "报告/main.pdf") == pdf_hash
    assert_unchanged(sequence)
    assert_unchanged(methods)
    ready_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "ready", "stage": "research-v3 report", "files": len(files), "pdf_pages": review["page_count"], "pdf_sha256": pdf_hash}))


if __name__ == "__main__":
    main()
