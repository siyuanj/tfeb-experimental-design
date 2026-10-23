#!/usr/bin/env python3
"""Exercise boundary cases in multipart feature drawing before real rendering."""

import json
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from render_plasmid_map import MapError, prepare_features, restriction_sites


def main():
    # A positive-strand origin-spanning CDS has its terminal part at the start.
    record = SeqRecord(Seq("A" * 1000), id="synthetic_geometry_check")
    record.features = [SeqFeature(CompoundLocation([
        SimpleLocation(900, 1000, strand=1),
        SimpleLocation(0, 60, strand=1),
    ]), type="CDS", qualifiers={"label": ["origin-spanning"]})]
    graphics, shown, omitted = prepare_features(record, {"source"})
    assert [(g.start, g.end) for g in graphics] == [(900, 1000), (0, 60)]
    assert sum(g.end - g.start for g in graphics) == 160
    assert [g.strand for g in graphics] == [0, 1]
    assert not omitted and len(shown) == 2

    # Reverse-strand join order is extraction order, not genomic sorting.
    record.features = [SeqFeature(CompoundLocation([
        SimpleLocation(0, 60, strand=-1),
        SimpleLocation(900, 1000, strand=-1),
    ]), type="CDS")]
    graphics, _, _ = prepare_features(record, set())
    assert [g.strand for g in graphics] == [0, -1]
    assert [(g.start, g.end) for g in graphics] == [(0, 60), (900, 1000)]

    # Non-contiguous internal annotations must preserve the actual gap.
    record.features = [SeqFeature(CompoundLocation([
        SimpleLocation(200, 250, strand=1),
        SimpleLocation(400, 450, strand=1),
    ]), type="misc_feature")]
    graphics, _, _ = prepare_features(record, set())
    assert sum(g.end - g.start for g in graphics) == 100

    record.features[0].location = CompoundLocation([
        SimpleLocation(200, 250, strand=1),
        SimpleLocation(400, 450, strand=-1),
    ])
    try:
        prepare_features(record, set())
    except MapError:
        pass
    else:
        raise AssertionError("Mixed-strand join was not rejected")

    # EcoRI recognition motif starts at base 2; cleavage labels base 3.
    record = SeqRecord(Seq("AGAATTCTTTTTTTT"), id="cut_coordinate_check")
    cuts = restriction_sites(record, ["EcoRI"])[0]["cuts"]
    assert len(cuts) == 1
    assert cuts[0]["first_base_after_cut_1based"] == 3
    assert cuts[0]["cut_boundary_0based"] == 2
    record.features = [SeqFeature(SimpleLocation(2, 8, strand=1), type="rep_origin",
                       qualifiers={"note": ["origin interval, no directional arrow assigned"]})]
    graphics, shown, _ = prepare_features(record, set())
    assert graphics[0].strand == 0 and shown[0]["parsed_location_strand"] == 1
    assert shown[0]["direction_from_note"] is True
    result = {
        "status": "passed",
        "scope": "synthetic drawing geometry and reference-strand cut labels only",
        "checks": ["positive origin join", "reverse origin join", "internal gap preserved",
                   "mixed-strand join rejected", "recognition start distinguished from cut",
                   "explicit non-directional origin note respected after GenBank round-trip"],
        "does_not_validate": "TFEB design, biological function, or experimental cleavage",
    }
    output = Path(__file__).resolve().parents[2] / "docs/协作/报告绘图边界核验.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
