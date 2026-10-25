#!/usr/bin/env python3
"""Check accepted inputs against rendered figures and extracted final PDF text.

This is a report consistency audit, not a replacement for sequence/method
validation or manual page inspection. It writes only its own compact report.
"""

import csv
import json
import re
import subprocess
import shlex
import sys
from datetime import datetime, timezone

from handoff import ROOT, assert_unchanged, sha256, verify_gate


def table(relative):
    with (ROOT / relative).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    output = ROOT / "docs/协作/报告一致性核验.json"
    output.write_text(json.dumps({"status": "running", "started_at_utc": datetime.now(timezone.utc).isoformat()}) + "\n")
    _, sequence = verify_gate("docs/协作/sequence-ready.json")
    _, methods = verify_gate("docs/协作/methods-ready.json")
    checks = {}
    build_path = ROOT / "docs/协作/报告编译核验.json"
    build_hash = sha256(build_path)
    build = json.loads(build_path.read_text())
    assert build["status"] == "passed"
    assert build["sequence_gate"]["sha256"] == sequence["sha256"]
    assert build["methods_gate"]["sha256"] == methods["sha256"]
    assert build["script_sha256"] == sha256(ROOT / "scripts/B-figures/compile_report.py")
    for relative, expected in build["inputs"].items():
        assert sha256(ROOT / relative) == expected, f"Changed compilation input: {relative}"
    for entry in (build["pdf"], build["pdf_text"], build["compiler_log"], *build["page_images"]):
        assert sha256(ROOT / entry["path"]) == entry["sha256"], entry["path"]
    checks["compilation_inputs_PDF_text_log_and_page_images_match_bound_build"] = True
    map_data = json.loads((ROOT / "报告/images/B/plasmid-map.json").read_text())
    pipeline = json.loads((ROOT / "报告/images/B/pipeline-figures.json").read_text())
    for entry in (map_data, pipeline):
        assert entry["sequence_gate"]["sha256"] == sequence["sha256"]
        for relative, expected in entry["outputs"].items():
            assert sha256(ROOT / relative) == expected, f"Changed figure: {relative}"
        for relative, expected in entry["scripts"].items():
            assert sha256(ROOT / relative) == expected, f"Changed figure source: {relative}"
    assert pipeline["methods_gate"]["sha256"] == methods["sha256"]
    checks["figure_input_and_output_hashes"] = True
    publication = json.loads((ROOT / '报告/images/B/publication-figures.json').read_text())
    for relative, expected in {**publication['sources'], **publication['outputs']}.items():
        assert sha256(ROOT / relative) == expected, relative
    control_manifest = ROOT / '设计/controls-v2/manifest.json'
    controls = json.loads(control_manifest.read_text())
    control_audit = json.loads((control_manifest.parent/'independent-verification.json').read_text())
    assert control_audit['status'] == 'passed'
    assert control_audit['manifest_sha256'] == sha256(control_manifest)
    assert control_audit['audit_script_sha256'] == sha256(ROOT/'scripts/B-audit/check_control_panel.py')
    for relative, expected in {**controls['source_inputs'], **controls['outputs']}.items():
        assert sha256(ROOT / relative) == expected, relative
    checks['native_maps_publication_figures_and_mutants_match_verified_inputs'] = True
    research_path = ROOT / '设计/research-v3/manifest.json'
    research = json.loads(research_path.read_text())
    research_audit = json.loads((research_path.parent/'independent-verification.json').read_text())
    assert research_audit['status'] == 'passed' and all(research_audit['checks'].values())
    assert research_audit['manifest_sha256'] == sha256(research_path)
    assert research_audit['audit_script_sha256'] == sha256(ROOT/'scripts/B-audit/check_research_panel.py')
    assert research_audit['translation_audit_sha256'] == sha256(ROOT/'scripts/B-audit/check_control_panel.py')
    for relative, expected in {**research['inputs'], **research['outputs']}.items():
        assert sha256(ROOT / relative) == expected, relative
    display = json.loads((ROOT/'报告/质粒图-SnapGene/display-inputs.json').read_text())
    assert len(display['records']) == 10
    for row in display['records']:
        assert sha256(ROOT/row['source']) == row['source_sha256']
        assert sha256(ROOT/row['display_file']) == row['display_sha256']
    checks['five_research_constructs_and_ten_native_map_inputs_match_audits'] = True
    condition_audit = json.loads((ROOT/'docs/协作/报告研究条件核验.json').read_text())
    assert condition_audit['status'] == 'passed' and all(condition_audit['checks'].values())
    for relative, expected in condition_audit['inputs'].items():
        assert sha256(ROOT/relative) == expected, relative
    checks['research_control_inventory_workload_and_contrast_algebra_are_current'] = True
    imagegen = json.loads((ROOT/'报告/images/B/imagegen-flow-v3.json').read_text())
    for relative, expected in {**imagegen['sources'], **imagegen['outputs']}.items():
        assert sha256(ROOT/relative) == expected,relative
    import xml.etree.ElementTree as ET
    for figure in imagegen['figures']:
        svg = ROOT/'报告/images/B'/(figure['name']+'.svg')
        xml = svg.read_text(); tree = ET.fromstring(xml)
        assert '<image' not in xml and 'base64' not in xml and 'foreignObject' not in xml
        assert len(list(tree.iter('{http://www.w3.org/2000/svg}text')))==figure['editable_text']
        assert len(list(tree.iter('{http://www.w3.org/2000/svg}path')))>=figure['traced_paths']
    checks['ImageGen_vector_paths_and_editable_labels_match_bound_sources'] = True

    features = table("设计/construct_features.tsv")
    shown = map_data["features_displayed"]
    assert len(shown) == len(features), "Final design has unexpected multipart or omitted features"
    for row in features:
        matches = [item for item in shown if item["label"] == row["name"]]
        assert len(matches) == 1, row["name"]
        item = matches[0]
        assert item["start_1based"] == int(row["start_1based"])
        assert item["end_inclusive_1based"] == int(row["end_1based"])
        assert item["strand"] == int(row["strand"]), row["name"]
    checks["all_feature_intervals_and_displayed_strands_match_accepted_table"] = True
    metrics = json.loads((ROOT / "设计/construct_metrics.json").read_text())
    assert map_data["length_bp"] == metrics["plasmid_bp"]
    for enzyme in map_data["restriction_enzymes"]:
        expected_after = metrics["construct_top_cut_after_1based"][enzyme["enzyme"]]
        assert [x["first_base_after_cut_1based"] for x in enzyme["cuts"]] == [x + 1 for x in expected_after]
    checks["cleavage_labels_distinguished_from_recognition_starts"] = True

    fragments = table("设计/expected_fragments.tsv")
    for panel in pipeline["gel"]["panels"]:
        for lane in panel["lanes"]:
            expected = [int(x["size_bp"]) for x in fragments
                        if x["assay"] == lane["assay"] and x["sample"] == lane["sample"] and int(x["size_bp"]) > 0]
            assert lane["bands_bp"] == expected, lane
    expected_marker = [int(x["size_bp"]) for x in fragments if x["assay"] == "marker"]
    assert pipeline["gel"]["marker_bp"] == expected_marker
    checks["all_plotted_lane_sizes_and_marker_sizes_match_accepted_table"] = True

    pdf = ROOT / "报告/main.pdf"
    text_path = ROOT / "报告/build/main.txt"
    content = text_path.read_text()
    collapsed = re.sub(r"\s+", "", content)
    for row in table("设计/primers.tsv"):
        assert row["sequence_5to3"] in collapsed, f"Cloning primer absent from PDF: {row['name']}"
    for row in table("设计/sanger_primers.tsv"):
        assert row["sequence_5to3"] in collapsed, f"Sanger primer absent from PDF: {row['primer_id']}"
    for phrase in ("Author Contributions", "AI Use Declaration", "References"):
        assert phrase in content, phrase
    assert "pEGFP-N1-TFEB-Bv1" in collapsed
    assert "80.42" in content
    assert "??" not in content
    checks["PDF_contains_full_cloning_and_Sanger_oligos_and_required_sections"] = True
    log = (ROOT / "报告/build/tectonic-output.txt").read_text()
    serious = [line for line in log.splitlines()
               if re.search(r"Overfull|Undefined control sequence|(?:citation|reference).*undefined|undefined.*(?:citation|reference)|Missing character|^error:", line, re.I)]
    assert not serious, serious
    checks["no_serious_compile_messages_in_captured_output"] = True
    info = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)", info, re.M).group(1))
    assert pages > 0
    assert pages == build["pdf"]["pages"]
    assert_unchanged(sequence)
    assert_unchanged(methods)
    assert sha256(build_path) == build_hash, "Compilation record changed during audit"
    for relative, expected in build["inputs"].items():
        assert sha256(ROOT / relative) == expected, relative
    for entry in (build["pdf"], build["pdf_text"], build["compiler_log"], *build["page_images"]):
        assert sha256(ROOT / entry["path"]) == entry["sha256"], entry["path"]
    result = {
        "status": "passed", "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "actual_command": shlex.join([sys.executable, str(__file__)]),
        "scope": "report consistency only; manual page review is separate",
        "checks": checks, "sequence_gate": sequence, "methods_gate": methods,
        "compilation_record": {"path": "docs/协作/报告编译核验.json", "sha256": build_hash},
        "pdf": {"path": "报告/main.pdf", "sha256": sha256(pdf), "pages": pages},
        "pdf_text": {"path": "报告/build/main.txt", "sha256": sha256(text_path)},
        "script_sha256": sha256(ROOT / "scripts/B-figures/check_report_consistency.py"),
        "limitations": ["Does not revalidate all prose claims or calculations",
                        "Does not establish experimental success or human author review",
                        "Visual legibility and complete pagination require final page images"],
    }
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"status": "passed", "checks": len(checks), "pages": pages, "pdf_sha256": sha256(pdf)}))


if __name__ == "__main__":
    main()
