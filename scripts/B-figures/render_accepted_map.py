#!/usr/bin/env python3
"""Render the P5-accepted annotated plasmid without modifying its sequence."""

import json
import platform
import shlex
import sys
from datetime import datetime, timezone
from importlib.metadata import version

from figure_layouts import save_figure
from handoff import ROOT, assert_unchanged, sha256, verify_gate
from render_plasmid_map import draw_map, prepare_features, read_input, restriction_sites


def main():
    gate, snapshot = verify_gate("docs/协作/sequence-ready.json")
    source = ROOT / "设计/construct.gb"
    if "设计/construct.gb" not in snapshot["files"]:
        raise ValueError("Accepted gate does not bind the annotated construct")
    spec = json.loads((ROOT / "设计/construct_spec.json").read_text())
    record, source_format, source_hash = read_input(source)
    graphics, shown, omitted = prepare_features(record, {"source"})
    palette = {"TFEB": "#149A8A", "EGFP": "#799E28", "linker GDPPVAT": "#C78B2D",
               "TFEB-EGFP fusion ORF": "#326EB5", "KanR/NeoR": "#9962A9"}
    for graphic, feature in zip(graphics, shown):
        if feature["label"] in palette:
            feature["color"] = palette[feature["label"]]
            graphic.color = feature["color"]
    enzymes = restriction_sites(record, [item["name"] for item in spec["enzymes"]])
    fig = draw_map(record, graphics, shown, enzymes, spec["construct_name"])
    assert_unchanged(snapshot)
    prefix = ROOT / "报告/images/B/plasmid-map"
    save_figure(fig, prefix)
    outputs = {str(prefix.with_suffix("." + ext).relative_to(ROOT)):
               sha256(prefix.with_suffix("." + ext)) for ext in ("pdf", "svg", "png")}
    manifest = {
        "status": "rendered; visual review recorded separately", "kind": "predicted annotated sequence map",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "actual_command": shlex.join([sys.executable, str(__file__)]),
        "design_version": gate["design_version"], "sequence_gate": snapshot,
        "source": {"path": "设计/construct.gb", "format": source_format, "sha256": source_hash},
        "length_bp": len(record), "features_displayed": shown, "features_omitted": omitted,
        "restriction_enzymes": enzymes,
        "software": {"python": platform.python_version(), **{
            name: version(name) for name in ("biopython", "dna_features_viewer", "matplotlib")}},
        "scripts": {str(p.relative_to(ROOT)): sha256(p) for p in [
            ROOT / "scripts/B-figures/render_accepted_map.py", ROOT / "scripts/B-figures/render_plasmid_map.py",
            ROOT / "scripts/B-figures/handoff.py", ROOT / "scripts/B-figures/figure_layouts.py"]},
        "conventions": ["Feature coordinates are 0-based half-open, with 1-based inclusive equivalents",
                        "Enzyme labels are the first 1-based nucleotide after reference-strand cleavage",
                        "An explicit no-direction note suppresses the origin arrow after GenBank round-trip",
                        "PolyA annotations mark signal motifs, not the entire processing region",
                        "Overlapping fusion ORF and component annotations are intentional"],
        "outputs": outputs,
    }
    assert_unchanged(snapshot)
    prefix.with_suffix(".json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"length_bp": len(record), "displayed_parts": len(shown), "outputs": outputs}, ensure_ascii=False))


if __name__ == "__main__":
    main()
