#!/usr/bin/env python3
"""Independently review the B-v1 P6 outputs without writing any files.

This preserves the 48 checks used in the interactive review on 2026-10-23.
It imports neither build_methods nor sequence_checks. Primer sites are located
with string.find on both orientations, read coordinates are expanded from each
primer's one-based 3-prime coordinate, restriction cuts come from Biopython,
and molarity is recalculated with Decimal and explicit unit conversion.

Run with the project-authorized interpreter and -B. JSON is written to stdout;
the caller may redirect it to an approved audit file. Any failed check returns
exit status 1. Missing or malformed inputs raise an exception rather than
producing a successful report. The expected values bind this review to B-v1.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import sys
from collections import Counter
from datetime import datetime
from decimal import Decimal
from importlib.metadata import version
from pathlib import Path
from zoneinfo import ZoneInfo

from Bio import SeqIO
from Bio.Restriction import BamHI, HindIII
from Bio.Seq import Seq


ROOT = Path(__file__).resolve().parents[2]
METHOD_FILES = [
    "设计/sanger_primers.tsv",
    "设计/sanger_coverage.tsv",
    "设计/expected_fragments.tsv",
    "设计/methods_parameters.json",
]


def digest_bytes(data: bytes) -> str:
    """Hash exact file bytes, without sequence normalization."""
    return hashlib.sha256(data).hexdigest()


def hits(dna: str, motif: str) -> list[int]:
    """Return one-based circular matches using literal string searches."""
    result = []
    extended = dna + dna[:len(motif) - 1]
    position = extended.find(motif)
    while 0 <= position < len(dna):
        result.append(position + 1)
        position = extended.find(motif, position + 1)
    return result


def orientations(dna: str, oligo: str) -> list[tuple[int, int]]:
    """Locate an oligo and its reverse complement independently."""
    reverse = str(Seq(oligo).reverse_complement())
    return ([(start, 1) for start in hits(dna, oligo)]
            + [(start, -1) for start in hits(dna, reverse)])


def products(dna: str, first: str, second: str) -> list[int]:
    """Count bases between facing primers' 5-prime endpoints on a circle.

    Full-length exact matches are required. Products longer than 4,000 bases
    are outside the screening model used by the P6 tables. This does not model
    mismatched annealing, amplification efficiency, or primer interactions.
    """
    lengths = []
    for left, right in ((first, second), (second, first)):
        for forward_start, forward_direction in orientations(dna, left):
            for reverse_start, reverse_direction in orientations(dna, right):
                if forward_direction != 1 or reverse_direction != -1:
                    continue
                reverse_5prime = reverse_start + len(right) - 1
                length = (reverse_5prime - forward_start) % len(dna) + 1
                if len(left) + len(right) <= length <= min(4000, len(dna)):
                    lengths.append(length)
    return sorted(lengths)


def restriction_fragments(dna: str, circular: bool) -> list[int]:
    """Use Biopython cut positions, not the generator's literal cut offsets."""
    cuts = sorted(HindIII.search(Seq(dna), linear=not circular)
                  + BamHI.search(Seq(dna), linear=not circular))
    if circular:
        gaps = [right - left for left, right in
                zip(cuts, cuts[1:] + [cuts[0] + len(dna)])]
    else:
        boundaries = [1] + cuts + [len(dna) + 1]
        gaps = [right - left for left, right in
                zip(boundaries, boundaries[1:])]
    return sorted(gaps, reverse=True)


def decimal(value: object) -> Decimal:
    """Retain decimal quantities without binary floating-point comparison."""
    return Decimal(str(value))


def main() -> int:
    manifest_path = "docs/协作/sequence-ready.json"
    manifest_bytes = (ROOT / manifest_path).read_bytes()
    manifest = json.loads(manifest_bytes)
    files = list(manifest["files"]) + [manifest_path] + METHOD_FILES
    raw = {path: (ROOT / path).read_bytes() for path in files}
    # Preserve the same manifest bytes that were parsed above.
    raw[manifest_path] = manifest_bytes
    hashes = {path: digest_bytes(data) for path, data in raw.items()}
    checks = {
        "P5_status_ready": manifest["status"] == "ready",
        "all_P5_manifest_hashes_match": all(
            hashes[path] == expected
            for path, expected in manifest["files"].items()),
    }

    def table(path: str) -> list[dict]:
        return list(csv.DictReader(io.StringIO(raw[path].decode()), delimiter="\t"))

    def sequence(path: str, file_format: str = "fasta") -> str:
        return str(SeqIO.read(io.StringIO(raw[path].decode()), file_format).seq).upper()

    construct = sequence("设计/construct.fasta")
    vector = sequence("资料/序列/pEGFP-N1_U55762.1.gb", "genbank")
    amplicon = sequence("设计/PCR_product.fasta")
    primers = table("设计/sanger_primers.tsv")
    coverage = table("设计/sanger_coverage.tsv")
    fragments = table("设计/expected_fragments.tsv")
    parameters = json.loads(raw["设计/methods_parameters.json"])
    target = set(range(parameters["sanger"]["target_start_1based"],
                       parameters["sanger"]["target_end_1based"] + 1))

    read_sets = {}
    details = []
    for primer in primers:
        identifier = primer["primer_id"]
        oligo = primer["sequence_5to3"]
        start = int(primer["binding_start_1based"])
        end = int(primer["binding_end_1based"])
        step = 1 if primer["direction"] == "+" else -1
        located = orientations(construct, oligo)
        bound = construct[start - 1:end]
        if step == -1:
            bound = str(Seq(bound).reverse_complement())
        three_prime = end if step == 1 else start
        usable = {three_prime + step * offset for offset in range(51, 651)}
        read_sets[identifier] = usable
        row = next(row for row in coverage if row["primer_id"] == identifier)
        overlap = usable & target
        checks[f"{identifier}_binding_unique_and_coordinates"] = (
            located == [(start, step)] and bound == oligo
            and end - start + 1 == len(oligo) == int(primer["length_nt"]))
        checks[f"{identifier}_coverage_coordinates"] = (
            len(usable) == 600
            and min(usable) == int(primer["usable_start_1based"])
            == int(row["usable_start_1based"])
            and max(usable) == int(primer["usable_end_1based"])
            == int(row["usable_end_1based"])
            and len(overlap) == int(row["target_overlap_nt"])
            and min(overlap) == int(row["target_overlap_start_1based"])
            and max(overlap) == int(row["target_overlap_end_1based"]))
        details.append({"primer": identifier,
                        "binding": [start, end, primer["direction"]],
                        "usable": [min(usable), max(usable)],
                        "target_overlap": len(overlap)})

    depth = {direction: Counter(
        position for primer in primers if primer["direction"] == direction
        for position in read_sets[primer["primer_id"]]) for direction in "+-"}
    checks["full_2198nt_target_bidirectional"] = (
        len(target) == 2198
        and all(depth[direction][position] >= 1
                for direction in "+-" for position in target))
    checks["primer_binding_bases_covered_from_opposite_direction"] = all(
        all(depth["-" if primer["direction"] == "+" else "+"][position] >= 1
            for position in range(int(primer["binding_start_1based"]),
                                  int(primer["binding_end_1based"]) + 1))
        for primer in primers)
    overlaps = {}
    for direction in "+-":
        ordered = [primer for primer in primers if primer["direction"] == direction]
        overlaps[direction] = [len(read_sets[left["primer_id"]]
                                   & read_sets[right["primer_id"]])
                               for left, right in zip(ordered, ordered[1:])]
    checks["successive_reads_overlap_at_least_50nt"] = all(
        size >= 50 for sizes in overlaps.values() for size in sizes)

    reference = SeqIO.read(io.StringIO(raw["资料/序列/TFEB_reference.gb"].decode()),
                           "genbank")
    cds = str(next(feature for feature in reference.features
                   if feature.type == "CDS").extract(reference.seq)).upper()
    tfeb = cds[:-3]
    begin = construct.find(tfeb)
    end = begin + len(tfeb)
    models = {
        "correct_construct": construct,
        "empty_vector": vector,
        "reverse_TFEB_model": (construct[:begin]
                               + str(Seq(tfeb).reverse_complement()) + construct[end:]),
        "retained_stop_model": construct[:end] + cds[-3:] + construct[end:],
    }
    by_name = {primer["primer_id"]: primer["sequence_5to3"] for primer in primers}
    predictions = {}
    for assay, data in parameters["colony_PCR"]["pairs"].items():
        first, second = [by_name[name] for name in data["primers"]]
        for sample, dna in models.items():
            actual = products(dna, first, second)
            reported = sorted(int(row["size_bp"]) for row in fragments
                              if row["assay"] == assay and row["sample"] == sample
                              and int(row["size_bp"]) > 0)
            checks[f"{assay}_{sample}"] = actual == reported
            predictions[f"{assay}/{sample}"] = actual
    for sample, dna in models.items():
        checks[f"digest_{sample}"] = restriction_fragments(dna, True) == sorted(
            [int(row["size_bp"]) for row in fragments
             if row["assay"] == "diagnostic_digest" and row["sample"] == sample],
            reverse=True)
    checks["preparative_PCR_digest"] = restriction_fragments(amplicon, False) == [1439, 11, 7]
    checks["empty_vector_digest"] = restriction_fragments(vector, True) == [4695, 38]

    cloning = table("设计/primers.tsv")
    forward = next(primer for primer in cloning if primer["direction"] == "forward")
    reverse = next(primer for primer in cloning if primer["direction"] == "reverse")
    forward_core = forward["annealing_sequence_5to3"]
    reverse_core = reverse["annealing_sequence_5to3"]
    # Use both complete oligos around the intervening template sequence. The
    # generator instead adds separately extracted tails to the complete CDS.
    independent_pcr = (forward["sequence_5to3"]
                       + tfeb[len(forward_core):len(tfeb) - len(reverse_core)]
                       + str(Seq(reverse["sequence_5to3"]).reverse_complement()))
    checks["full_cloning_primers_and_PCR"] = (
        forward["sequence_5to3"].endswith(forward_core)
        and reverse["sequence_5to3"].endswith(reverse_core)
        and tfeb.startswith(forward_core)
        and tfeb.endswith(str(Seq(reverse_core).reverse_complement()))
        and independent_pcr == amplicon and len(amplicon) == 1457)
    protein = str(Seq(construct[630:2799]).translate())
    checks["fusion_translation_722aa_and_terminal_stop"] = (
        protein.endswith("*") and "*" not in protein[:-1]
        and len(protein[:-1]) == 722
        and protein[:-1] == sequence("设计/fusion_protein.fasta"))

    ligation = parameters["ligation"]
    # ng -> g -> mol -> fmol gives the factor 1,000,000 below.
    vector_fmol = (decimal(ligation["vector_ng"]) * Decimal(1000000)
                   / (decimal(ligation["backbone_bp"]) * 660))
    insert_fmol = (decimal(ligation["insert_ng"]) * Decimal(1000000)
                   / (decimal(ligation["insert_bp"]) * 660))
    checks["ligation_molarity_20_60_fmol_and_ratio3"] = (
        vector_fmol == 20 and insert_fmol == 60 and insert_fmol / vector_fmol == 3)
    checks["ligation_volume_and_concentration"] = (
        sum(decimal(ligation[key]) for key in
            ["vector_uL", "insert_uL", "buffer10X_uL", "ligase_uL", "water_uL"]) == 20
        and decimal(ligation["vector_uL"])
        * decimal(ligation["vector_stock_target_ng_per_uL"]) == decimal(ligation["vector_ng"])
        and decimal(ligation["insert_uL"])
        * decimal(ligation["insert_stock_target_ng_per_uL"]) == decimal(ligation["insert_ng"])
        and (decimal(ligation["vector_ng"]) + decimal(ligation["insert_ng"])) / 20
        == decimal(ligation["total_DNA_ng_per_uL"]))
    reaction = parameters["cloning_PCR"]
    checks["cloning_PCR_50uL"] = (
        decimal(reaction["mix_2X_uL"]) + 2 * decimal(reaction["primer_10uM_each_uL"])
        + decimal(reaction["template_uL"]) + decimal(reaction["water_uL"]) == 50)
    reaction = parameters["digest"]
    checks["digest_50uL"] = (
        2 * decimal(reaction["each_uL"]) + decimal(reaction["buffer10X_uL"])
        + decimal(reaction["DNA_uL"]) + decimal(reaction["water_uL"]) == 50)
    reaction = parameters["colony_PCR"]
    checks["colony_PCR_25uL"] = (
        decimal(reaction["mix2X_uL"]) + 2 * decimal(reaction["primer10uM_each_uL"])
        + decimal(reaction["water_uL"]) == 25)
    checks["parameters_all_bound_input_hashes_match"] = all(
        digest_bytes((ROOT / path).read_bytes()) == expected
        for path, expected in parameters["inputs"].items())
    checks["no_inputs_changed_during_readonly_check"] = all(
        digest_bytes((ROOT / path).read_bytes()) == expected
        for path, expected in hashes.items())

    failed = [name for name, passed in checks.items() if not passed]
    result = {
        "status": "passed" if not failed else "failed",
        "read_only_independent_review": True,
        "checked_at": datetime.now(ZoneInfo("America/New_York")).isoformat(),
        "independent_path": (
            "No project-module imports; string.find exact primer matches in both "
            "orientations; one-based primer 3-prime coordinate expansion; "
            "Biopython restriction cuts; Decimal unit conversion"),
        "passed": sum(checks.values()), "checks": len(checks),
        "failed": failed, "check_results": checks,
        "input_hashes": hashes,
        "script_sha256": digest_bytes(Path(__file__).read_bytes()),
        "software": {"python": sys.version.split()[0], "biopython": version("biopython")},
        "command": ("'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' "
                    "-B scripts/B-methods/independent_review.py"),
        "method_hashes": {path: hashes[path] for path in METHOD_FILES},
        "construct_bp": len(construct),
        "target_1based": [min(target), max(target)], "target_nt": len(target),
        "minimum_depth": {direction: min(depth[direction][position] for position in target)
                          for direction in "+-"},
        "adjacent_read_overlap_nt": overlaps, "primer_checks": details,
        "PCR_predictions_bp": predictions,
        "digest_bp": {name: restriction_fragments(dna, True) for name, dna in models.items()},
        "molarity_fmol": {"vector": str(vector_fmol), "insert": str(insert_fmol)},
        "ligation_ng": {"vector": ligation["vector_ng"], "insert": ligation["insert_ng"]},
        "ligation_uL": {key: ligation[key] for key in
                        ["vector_uL", "insert_uL", "water_uL", "total_uL"]},
        "limitations": [
            "Exact full-length oligo matches only; no genome-wide off-target or mismatch PCR model",
            "Predicted offsets 51-650, not actual Sanger reads",
            "No wet-lab validation",
        ],
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
