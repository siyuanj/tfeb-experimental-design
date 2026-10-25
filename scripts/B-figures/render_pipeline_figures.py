#!/usr/bin/env python3
"""Render P6-accepted cloning workflow and predicted DNA-band schematics."""

import csv
import json
import shlex
import sys
from datetime import datetime, timezone

from figure_layouts import predicted_gel, save_figure, workflow
from handoff import ROOT, assert_unchanged, sha256, verify_gate


def main():
    sequence, seq_snapshot = verify_gate("docs/协作/sequence-ready.json")
    methods, meth_snapshot = verify_gate("docs/协作/methods-ready.json")
    metrics = json.loads((ROOT / "设计/construct_metrics.json").read_text())
    spec = json.loads((ROOT / "设计/construct_spec.json").read_text())
    parameters = json.loads((ROOT / "设计/methods_parameters.json").read_text())
    with (ROOT / "设计/expected_fragments.tsv").open() as handle:
        fragments = list(csv.DictReader(handle, delimiter="\t"))
    for filename in ("设计/expected_fragments.tsv", "设计/methods_parameters.json"):
        if filename not in meth_snapshot["files"]:
            raise ValueError(f"Method gate does not bind figure input: {filename}")
    if parameters["design_version"] != sequence["design_version"]:
        raise ValueError("Sequence and method versions differ")

    def sizes(assay, sample):
        matches = [int(row["size_bp"]) for row in fragments
                   if row["assay"] == assay and row["sample"] == sample]
        if not matches or any(x < 0 for x in matches):
            raise ValueError(f"Missing or invalid expected fragments: {assay}/{sample}")
        return [x for x in matches if x > 0]

    def lane(label, assay, sample, note=None):
        bands = sizes(assay, sample)
        if note is None:
            note = " + ".join(f"{x:,}" for x in bands) + " bp" if bands else "no specific band"
        return {"label": label, "assay": assay, "sample": sample, "bands_bp": bands, "note": note}

    flank_pair = " / ".join(parameters["colony_PCR"]["pairs"]["colony_flanking"]["primers"])
    orientation_pair = " / ".join(parameters["colony_PCR"]["pairs"]["colony_orientation"]["primers"])
    enzyme_names = " / ".join(entry["name"] for entry in spec["enzymes"])
    panels = [
        {"title": "A  Insert PCR and digestion", "lanes": [
            lane("PCR", "cloning_PCR", "synthetic_TFEB_template"),
            lane("NTC", "cloning_PCR", "no_template_control"),
            lane("I", "preparative_digest", "PCR_product", f"{metrics['digested_insert_nominal_nt']:,} bp; terminal pieces omitted")
        ]},
        {"title": "B  Flanking colony PCR", "lanes": [
            lane("E", "colony_flanking", "empty_vector"),
            lane("C", "colony_flanking", "correct_construct"),
            lane("R", "colony_flanking", "reverse_TFEB_model"),
            lane("S", "colony_flanking", "retained_stop_model"),
            lane("NTC", "colony_flanking", "no_template_control"),
        ]},
        {"title": "C  Orientation PCR", "lanes": [
            lane("E", "colony_orientation", "empty_vector"),
            lane("C", "colony_orientation", "correct_construct"),
            lane("R", "colony_orientation", "reverse_TFEB_model"),
            lane("S", "colony_orientation", "retained_stop_model"),
            lane("NTC", "colony_orientation", "no_template_control"),
        ]},
        {"title": f"D  {enzyme_names} diagnostic digest", "lanes": [
            lane("E", "diagnostic_digest", "empty_vector"),
            lane("C", "diagnostic_digest", "correct_construct"),
            lane("R", "diagnostic_digest", "reverse_TFEB_model"),
            lane("S", "diagnostic_digest", "retained_stop_model"),
        ]},
    ]
    marker = sizes("marker", "NEB N3232")
    tick_labels = {10002: "10,002", 5001: "5,001", 2000: "2,000", 1000: "1,000", 500: "500/517"}
    if not set(tick_labels).issubset(marker) or 517 not in marker:
        raise ValueError("Marker changed; re-review grouping and tick-label choices")
    # Both manufacturer marker bands remain in the data, but intentionally
    # overlap visually. The single tick label does not imply their resolution.
    gel = predicted_gel(panels, marker, lower_bp=400, upper_bp=12000,
                        marker_tick_labels=tick_labels)
    flow_nodes = [
        ("1  Defined native TFEB template", f"{metrics['TFEB_reference_CDS_with_stop_nt']:,} nt synthetic dsDNA; {spec['tfeb']['transcript']}"),
        ("2  High-fidelity PCR", f"Q5 M0492; {metrics['PCR_product_bp']:,} bp tailed product"),
        ("3  Cleanup and separate double digests", f"{enzyme_names}; purify backbone and insert"),
        ("4  Directional T4 ligation", f"{metrics['digested_backbone_nominal_nt']:,} nt backbone + {metrics['digested_insert_nominal_nt']:,} nt insert; {parameters['ligation']['insert_to_vector_molar_ratio']}:1 insert:vector"),
        ("5  Transform and select", f"NEB 5-alpha; kanamycin {parameters['transformation']['kanamycin_ug_per_mL']} ug/mL"),
        ("6  Screen and prepare candidates", "Flanking/orientation PCR; miniprep; diagnostic digest"),
        ("7  Sanger acceptance before cell assays", f"{parameters['sanger']['primer_count']} primers; full target and both junctions"),
    ]
    controls = [
        "Planned procurement only. Check supplied sequence and stock identity before use.",
        "No-template control; target band purification if amplification is not specific.",
        "Remove enzymes and residual circular vector; discard short terminal/MCS fragments.",
        "Vector-only ligation estimates background. DNA amounts use measured concentrations.",
        "No-DNA control; separate pUC19 competence control uses ampicillin.",
        "Size and orientation screens do not exclude point mutations or a retained stop.",
        "Accept unambiguous concordant bidirectional reads. Rescreen or resequence failures.",
    ]
    flow = workflow(flow_nodes, controls, "PROPOSED CLONING WORKFLOW - SCHEMATIC")
    for snapshot in (seq_snapshot, meth_snapshot):
        assert_unchanged(snapshot)
    directory = ROOT / "报告/images/B"
    save_figure(gel, directory / "predicted-bands")
    save_figure(flow, directory / "cloning-workflow")
    outputs = {str(path.relative_to(ROOT)): sha256(path)
               for name in ("predicted-bands", "cloning-workflow")
               for ext in ("pdf", "svg", "png")
               for path in [directory / (name + "." + ext)]}
    manifest = {
        "status": "rendered; visual review recorded separately",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "actual_command": shlex.join([sys.executable, str(__file__)]),
        "sequence_gate": seq_snapshot, "methods_gate": meth_snapshot,
        "gel": {"panels": panels, "marker_product": "NEB N3232", "marker_bp": marker,
                "axis_limits_bp": [400, 12000], "tick_labels": tick_labels,
                "limitations": ["Predicted positions only; all band strokes have equal width",
                                "500/517 marker pair is not resolved; 3-bp differences are not expected to resolve",
                                "No specific band is distinct from a guaranteed blank lane",
                                "Uncut circular plasmids are not assigned linear migration positions",
                                "Reverse/stop constructs are hypothetical challenge models, not observed clones"]},
        "workflow": {"nodes": flow_nodes, "control_notes": controls},
        "scripts": {str(path.relative_to(ROOT)): sha256(path) for path in [
            ROOT / "scripts/B-figures/render_pipeline_figures.py", ROOT / "scripts/B-figures/figure_layouts.py",
            ROOT / "scripts/B-figures/handoff.py", ROOT / "scripts/B-figures/publication_style.py",
            ROOT / "scripts/B-figures/scientific_schematics.py"]},
        "outputs": outputs,
    }
    for snapshot in (seq_snapshot, meth_snapshot):
        assert_unchanged(snapshot)
    (directory / "pipeline-figures.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"outputs": outputs}, ensure_ascii=False))


if __name__ == "__main__":
    main()
