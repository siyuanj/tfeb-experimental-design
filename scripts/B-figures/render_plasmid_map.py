#!/usr/bin/env python3
"""Render annotated circular DNA, preserving multipart feature geometry.

Adapted from the course renderer (SHA-256
e79564a42038eb9440236b2f09672de885f2ded1ba280d0b8a2ae3a31895b6c7).
The original course script is unchanged. Coordinates are never inferred from
labels; joined locations are displayed as their individual annotated parts.
"""

import argparse
import hashlib
import io
import json
import math
import os
import platform
import shutil
import sys
import tempfile
import textwrap
from copy import deepcopy
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Restriction, SeqIO
from Bio.SeqFeature import CompoundLocation, ExactPosition
from dna_features_viewer import CircularGraphicRecord, GraphicFeature


class MapError(ValueError):
    """An input cannot be rendered without ambiguity or data loss."""


PALETTE = {
    "CDS": "#00A087",
    "promoter": "#3C5488",
    "rep_origin": "#B09C85",
    "terminator": "#F39B7F",
    "polyA_signal": "#F39B7F",
    "mRNA": "#91D1C2",
    "primer_bind": "#888888",
    "misc_feature": "#8491B4",
    "misc_RNA": "#4DBBD5",
}
FALLBACK_COLOR = "#8491B4"
CUT_COLOR = "#E64B35"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def feature_name(feature):
    for key in ("label", "gene", "name", "product", "locus_tag", "note"):
        values = feature.qualifiers.get(key)
        if values:
            value = values if isinstance(values, str) else " | ".join(map(str, values))
            return " ".join(value.split())
    return feature.type


def read_input(path):
    raw = path.read_bytes()
    suffix = path.suffix.lower()
    if suffix == ".dna":
        record = SeqIO.read(io.BytesIO(raw), "snapgene")
        file_format = "snapgene"
    elif suffix in {".gb", ".gbk", ".genbank"}:
        record = SeqIO.read(io.StringIO(raw.decode("utf-8-sig")), "genbank")
        file_format = "genbank"
    else:
        raise MapError("Input must be a SnapGene .dna or GenBank .gb/.gbk/.genbank file.")
    if str(record.annotations.get("topology", "")).lower() != "circular":
        raise MapError("The input must explicitly declare circular topology; no topology is inferred.")
    if not record.seq.defined or len(record) == 0:
        raise MapError("A nonempty, fully defined DNA sequence is required.")
    sequence = str(record.seq).upper()
    if set(sequence) - set("ACGTRYSWKMBDHVN"):
        raise MapError("The record contains unsupported DNA symbols or sequence gaps.")
    return record, file_format, hashlib.sha256(raw).hexdigest()


def prepare_features(record, omitted_types):
    graphic_features, shown, omitted = [], [], []
    for index, feature in enumerate(record.features):
        item = {
            "feature_index_0based": index,
            "type": feature.type,
            "label": feature_name(feature),
            "input_location": str(feature.location),
            "qualifiers": feature.qualifiers,
        }
        if feature.type in omitted_types:
            omitted.append({**item, "reason": "feature type explicitly omitted from the drawing"})
            continue
        if feature.location is None:
            omitted.append({**item, "reason": "input annotation has no location"})
            continue
        location = feature.location
        if isinstance(location, CompoundLocation) and location.operator != "join":
            raise MapError(f"Feature {index} uses unsupported operator {location.operator}.")
        parts = list(location.parts)
        if len({part.strand for part in parts}) > 1:
            raise MapError(f"Feature {index} joins mixed strands; resolve explicitly.")
        for part_index, part in enumerate(parts):
            if type(part.start) is not ExactPosition or type(part.end) is not ExactPosition:
                raise MapError(f"Feature {index} ({item['label']}) has fuzzy/unknown coordinates.")
            if part.ref is not None or part.ref_db is not None:
                raise MapError(f"Feature {index} ({item['label']}) refers to a remote sequence.")
            start, end = int(part.start), int(part.end)
            if not 0 <= start < end <= len(record):
                raise MapError(f"Feature {index} has an empty, reversed, or out-of-range part.")
            parsed_strand = part.strand
            notes = " ".join(feature.qualifiers.get("note", []))
            # GenBank round-tripping serializes strand=0 as an unsigned interval,
            # which readers interpret as +1. Preserve the author's explicit
            # non-directional annotation instead of inventing an origin arrow.
            non_directional = "no directional arrow assigned" in notes
            strand = 0 if non_directional else parsed_strand
            if strand not in (-1, 0, 1, None):
                raise MapError(f"Feature {index} has an unsupported strand value.")
            color = PALETTE.get(feature.type, FALLBACK_COLOR)
            # GenBank join order is biological extraction order. Only the last
            # part carries the terminal arrow; gaps are never filled in.
            arrow_strand = (strand or 0) if part_index == len(parts) - 1 else 0
            graphic = GraphicFeature(
                start=start, end=end, strand=arrow_strand, label=None,
                color=color, linecolor="white", linewidth=0.8,
            )
            graphic_features.append(graphic)
            label = item["label"]
            if len(parts) > 1:
                label += f" [part {part_index + 1}/{len(parts)}]"
            shown.append({
                **item, "label": label, "part_index_0based": part_index,
                "part_count": len(parts), "start_0based": start,
                "end_exclusive_0based": end, "start_1based": start + 1,
                "end_inclusive_1based": end, "strand": strand,
                "parsed_location_strand": parsed_strand,
                "direction_from_note": non_directional,
                "arrow_strand": arrow_strand, "color": color,
            })
    if not graphic_features:
        raise MapError("No located annotations remain to draw; nothing has been exported.")
    return graphic_features, shown, omitted


def restriction_sites(record, names):
    names = list(dict.fromkeys(names))
    if names and set(str(record.seq).upper()) - set("ACGT"):
        raise MapError("Restriction-site reporting requires unambiguous A/C/G/T sequence.")
    results = []
    for name in names:
        enzyme = getattr(Restriction, name, None)
        if enzyme is None or enzyme not in Restriction.AllEnzymes:
            raise MapError(f"Unknown Bio.Restriction enzyme: {name}")
        cuts = []
        for raw_position in enzyme.search(record.seq, linear=False):
            before = (int(raw_position) - 1) % len(record) + 1
            cuts.append({
                "search_result_1based": int(raw_position),
                "first_base_after_cut_1based": before,
                "last_base_before_cut_1based": (before - 2) % len(record) + 1,
                "cut_boundary_0based": before - 1,
            })
        results.append({
            "enzyme": name, "recognition_sequence": str(enzyme.site),
            "cuts": cuts, "cut_count": len(cuts),
        })
    return results


def polar(position, length, radius):
    angle = math.pi / 2 - 2 * math.pi * position / length
    return radius * math.cos(angle), -1 + radius * math.sin(angle)


def spread_labels(labels, minimum, maximum):
    """Keep vertical order and reserve enough space for wrapped labels."""
    labels.sort(key=lambda item: item["preferred_y"])
    if not labels:
        return
    total_height = sum(item["height"] for item in labels)
    if total_height > maximum - minimum:
        raise MapError("Too many labels for a readable map; omit selected feature types explicitly.")
    cursor = minimum
    for item in labels:
        half = item["height"] / 2
        item["label_y"] = max(item["preferred_y"], cursor + half)
        cursor = item["label_y"] + half
    cursor = maximum
    for item in reversed(labels):
        half = item["height"] / 2
        item["label_y"] = min(item["label_y"], cursor - half)
        cursor = item["label_y"] - half


def draw_map(record, graphic_features, shown, enzymes, title):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "pdf.fonttype": 42,
        "svg.fonttype": "none", "figure.facecolor": "white",
        "axes.facecolor": "white", "text.color": "#222222",
    })
    fig, ax = plt.subplots(figsize=(10, 7.2))
    circular = CircularGraphicRecord(
        sequence_length=len(record), features=graphic_features,
        feature_level_height=0.14,
    )
    _, (levels, _) = circular.plot(ax=ax, with_ruler=False)
    outer_radius = 1 + max(levels.values(), default=0) * 0.14 + 0.13
    labels = []
    for graphic, item in zip(graphic_features, shown):
        position = (graphic.start + graphic.end) / 2
        anchor = polar(position, len(record), 1 + levels[graphic] * 0.14)
        text = textwrap.fill(textwrap.shorten(item["label"], width=60, placeholder="…"), width=24)
        item["display_label"] = text
        labels.append({
            "anchor": anchor, "preferred_y": polar(position, len(record), outer_radius + 0.25)[1],
            "text": text, "color": item["color"], "enzyme": False,
            "height": 0.13 * (text.count("\n") + 1) + 0.07,
        })
    for result in enzymes:
        for cut in result["cuts"]:
            position = cut["cut_boundary_0based"]
            start, end = polar(position, len(record), outer_radius), polar(position, len(record), outer_radius + 0.10)
            inner = polar(position, len(record), 0.97)
            ax.plot([inner[0], start[0]], [inner[1], start[1]], color=CUT_COLOR, lw=0.6, alpha=0.45, linestyle="--")
            ax.plot([start[0], end[0]], [start[1], end[1]], color=CUT_COLOR, lw=1.5)
            labels.append({
                "anchor": end, "preferred_y": polar(position, len(record), outer_radius + 0.25)[1],
                "text": f"{result['enzyme']}  {cut['first_base_after_cut_1based']:,}",
                "color": CUT_COLOR, "enzyme": True, "height": 0.24,
            })
    vertical_extent = max(outer_radius + 0.45, max(
        (sum(item["height"] for item in labels if (item["anchor"][0] >= 0) == side) / 2 + 0.15)
        for side in (False, True)
    ))
    for side in (-1, 1):
        group = [item for item in labels if (1 if item["anchor"][0] >= 0 else -1) == side]
        spread_labels(group, -1 - vertical_extent, -1 + vertical_extent)
        text_x = side * (outer_radius + 0.64)
        for item in group:
            x, y = item["anchor"]
            elbow_x = side * (outer_radius + 0.35)
            ax.plot([x, elbow_x, text_x - side * 0.06], [y, item["label_y"], item["label_y"]],
                    color=item["color"], lw=0.7, alpha=0.8, zorder=0)
            ax.text(text_x, item["label_y"], item["text"], ha="left" if side > 0 else "right",
                    va="center", fontsize=11.5, color=item["color"],
                    fontweight="bold" if item["enzyme"] else "normal", linespacing=1.15)
    # Inner ticks use 0-based boundary offsets, independently of enzyme labels.
    tick_step = max(1, math.ceil(len(record) / 8 / 100) * 100)
    for position in range(0, len(record), tick_step):
        start, end = polar(position, len(record), 0.89), polar(position, len(record), 0.94)
        ax.plot([start[0], end[0]], [start[1], end[1]], color="#BBBBBB", lw=0.65)
        x, y = polar(position, len(record), 0.75)
        ax.text(x, y, str(position), fontsize=8.5, color="#777777", ha="center", va="center")
    ax.text(0, -0.82, textwrap.fill(title, width=12, break_long_words=True), ha="center", va="center",
            fontsize=12, fontweight="bold", linespacing=1.3)
    ax.text(0, -1.20, f"{len(record):,} bp", ha="center", va="center", fontsize=12, color="#666666")
    ax.set_xlim(-outer_radius - 2.55, outer_radius + 2.55)
    ax.set_ylim(-1 - vertical_extent - 0.25, -1 + vertical_extent + 0.25)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.subplots_adjust(left=0.025, right=0.975, top=0.97, bottom=0.09)
    return fig


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Annotated circular .dna or GenBank file; never modified")
    parser.add_argument("--output-prefix", required=True, type=Path, help="Output path without extension; existing files are never overwritten")
    parser.add_argument("--title", help="Center title; defaults to the input filename stem")
    parser.add_argument("--enzymes", nargs="+", default=[], metavar="ENZYME", help="Bio.Restriction names, e.g. EcoRI XhoI")
    parser.add_argument("--omit-type", action="append", default=[], metavar="TYPE", help="Omit this annotation type from the drawing; source is always omitted")
    parser.add_argument("--export-genbank", action="store_true", help="Also export parsed sequence and annotations as standard .gb (not native SnapGene)")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution (72–1200; default 300)")
    return parser.parse_args()


def main():
    args = parse_args()
    if not 72 <= args.dpi <= 1200:
        raise MapError("--dpi must be between 72 and 1200.")
    source = args.input.expanduser().resolve(strict=True)
    prefix = args.output_prefix.expanduser().absolute()
    extensions = ["png", "svg", "pdf"] + (["gb"] if args.export_genbank else []) + ["json"]
    outputs = {extension: Path(f"{prefix}.{extension}") for extension in extensions}
    for output in outputs.values():
        if output.exists() or output.is_symlink():
            raise MapError(f"Output already exists and will not be overwritten: {output}")
    record, source_format, source_hash = read_input(source)
    features, shown, omitted = prepare_features(record, {"source", *args.omit_type})
    enzymes = restriction_sites(record, args.enzymes)
    title = args.title.strip() if args.title is not None else source.stem
    if not title:
        raise MapError("The map title cannot be empty.")
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {"path": str(source), "sha256": source_hash, "format": source_format,
                   "record_id": record.id, "record_name": record.name},
        "map_title": title, "sequence_length_bp": len(record), "topology": "circular",
        "sequence_sha256": hashlib.sha256(str(record.seq).upper().encode("ascii")).hexdigest(),
        "software": {"python": platform.python_version(), **{
            name: version(name) for name in ("biopython", "dna_features_viewer", "matplotlib")
        }},
        "renderer_script_sha256": sha256(Path(__file__).resolve()),
        "coordinate_conventions": {
            "features": "0-based half-open [start,end); equivalent 1-based inclusive coordinates are also recorded",
            "enzyme_labels": "Bio.Restriction.search(linear=False): 1-based first base after the cut on the reference strand; not recognition-site start",
            "cut_boundary": "0-based boundary before the labeled base, normalized modulo sequence length",
            "inner_scale": "0-based boundary offsets; zero is the reference origin",
        },
        "features_displayed": shown, "features_omitted": omitted, "restriction_enzymes": enzymes,
        "limitations": [
            "Restriction cuts are sequence-based predictions; methylation, reaction conditions, and experimental cleavage are not tested.",
            "Joined locations are split into their actual parts; only the biological terminal part has an arrowhead. Fuzzy coordinates and mixed-strand joins are rejected.",
            "The diagram is not a native SnapGene map and does not validate construct function or reading frames.",
        ],
        "genbank_export": "Parsed sequence and feature annotations only; proprietary SnapGene state is not preserved." if args.export_genbank else None,
        "outputs": {},
    }
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".plasmid-map-", dir=prefix.parent) as staging_name:
        staging = Path(staging_name)
        fig = draw_map(record, features, shown, enzymes, title)
        try:
            for extension in ("png", "svg", "pdf"):
                fig.savefig(staging / f"map.{extension}", dpi=args.dpi, facecolor="white", bbox_inches="tight")
        finally:
            plt.close(fig)
        if args.export_genbank:
            exported = deepcopy(record)
            exported.annotations.setdefault("molecule_type", "DNA")
            SeqIO.write(exported, staging / "map.gb", "genbank")
        if sha256(source) != source_hash:
            raise MapError("The source changed during rendering; no outputs have been published.")
        for extension in extensions:
            if extension != "json":
                staged = staging / f"map.{extension}"
                manifest["outputs"][extension] = {"path": str(outputs[extension]), "sha256": sha256(staged), "bytes": staged.stat().st_size}
        (staging / "map.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        # Reserve every destination exclusively before publishing; JSON is copied last.
        reserved = []
        try:
            for output in outputs.values():
                descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
                os.close(descriptor)
                reserved.append(output)
            for extension in extensions:
                shutil.copyfile(staging / f"map.{extension}", outputs[extension])
        except BaseException:
            # Only paths exclusively created by this invocation may be removed.
            for output in reserved:
                output.unlink(missing_ok=True)
            raise
    print(json.dumps({"outputs": {key: str(value) for key, value in outputs.items()},
                      "length_bp": len(record), "features_displayed": len(shown),
                      "features_omitted": len(omitted), "source_sha256": source_hash}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (MapError, OSError, ValueError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
