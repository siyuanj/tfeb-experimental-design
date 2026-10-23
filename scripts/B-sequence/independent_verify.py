#!/usr/bin/env python3
"""Audit B-v1 with literal sequence operations; never import the builder.

Only the Python standard library is used. Coordinates independently audited from
NM_007162.3 and U55762.1 are pinned below. The sole output is
design/independent_verification.json (the directory name is Chinese on disk).
This verifies an in-silico design, not a physical template or experimental result.
"""

import csv
import hashlib
import io
import json
import math
import re
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "设计" / "independent_verification.json"
DNA_BASES = set("ACGT")
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
TFEB_INTERVAL = (303, 1733)
EGFP_INTERVAL = (679, 1398)
PINNED_SEQUENCES = {
    "NM_007162.3": "0702312ceb2a00b30cc060d1ba809bfeeef5105c132be49512818041cc26e70e",
    "NM_007162.2": "080f92f539ae5015adefc68de1865d8ead5371931c5eb14cfbc9a955c1fb1a0f",
    "NP_009093.1": "84f910f982077df035fe31d4dc4baa3731c043c6efb3682ed2fac9cd807a0ae9",
    "U55762.1": "63918e4dcae6eed7f21a2d430f8b24a7172acc9a174e307d933550034189efde",
}
# NCBI standard genetic code in T,C,A,G order for each codon position.
# Source: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
CODON_AA = (
    "FFLLSSSSYY**CC*W"
    "LLLLPPPPHHQQRRRR"
    "IIIMTTTTNNKKSSRR"
    "VVVVAAAADDEEGGGG"
)
CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
GENETIC_CODE = dict(zip(CODONS, CODON_AA, strict=True))
# Literal average neutral free-amino-acid masses in Da. This independent
# calculation uses the published IUPAC/Biopython average-mass convention;
# it does not call Bio.SeqUtils.molecular_weight or read builder output values.
# Source: https://github.com/biopython/biopython/blob/biopython-188/Bio/Data/IUPACData.py
FREE_AA_MASS = {
    "A": 89.0932, "C": 121.1582, "D": 133.1027, "E": 147.1293,
    "F": 165.1891, "G": 75.0666, "H": 155.1546, "I": 131.1729,
    "K": 146.1876, "L": 131.1729, "M": 149.2113, "N": 132.1179,
    "P": 115.1305, "Q": 146.1445, "R": 174.2010, "S": 105.0926,
    "T": 119.1192, "V": 117.1463, "W": 204.2252, "Y": 181.1885,
}
WATER_MASS = 18.0153
ENZYMES = {
    "HindIII": ("AAGCTT", 1, 5, "AGCT"),
    "BamHI": ("GGATCC", 1, 5, "GATC"),
}


class AuditFailure(Exception):
    """Keep a machine-readable reason for positive and negative checks."""

    def __init__(self, code, message):
        self.code = code
        super().__init__(message)


def need(condition, code, message):
    if not condition:
        raise AuditFailure(code, message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def reverse_complement(sequence):
    need(set(sequence) <= DNA_BASES, "ambiguous_dna", "Non-ACGT DNA")
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def translate(sequence):
    need(len(sequence) % 3 == 0, "reading_frame", "ORF length is not divisible by three")
    need(set(sequence) <= DNA_BASES, "ambiguous_dna", "Non-ACGT ORF")
    return "".join(GENETIC_CODE[sequence[i:i + 3]] for i in range(0, len(sequence), 3))


def sites(sequence, motif, circular=False):
    search = sequence + sequence[:len(motif) - 1] if circular else sequence
    return [i for i in range(len(sequence)) if search.startswith(motif, i)]


def parse_genbank(text):
    origins = list(re.finditer(r"^ORIGIN[^\n]*\n(.*?)^//\s*$", text, re.M | re.S))
    need(len(origins) == 1, "genbank_format", "Expected exactly one GenBank record")
    sequence = re.sub(r"[\s0-9]", "", origins[0].group(1)).upper()
    locus = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+(bp|aa)\b([^\n]*)", text, re.M)
    version = re.search(r"^VERSION\s+(\S+)", text, re.M)
    need(locus is not None and version is not None, "genbank_format", "Missing LOCUS/VERSION")
    need(len(sequence) == int(locus[1]), "genbank_length", "LOCUS length differs from ORIGIN")
    allowed = DNA_BASES if locus[2] == "bp" else AMINO_ACIDS
    need(set(sequence) <= allowed, "genbank_alphabet", "Unexpected sequence characters")
    return {"sequence": sequence, "version": version[1], "circular": "circular" in locus[3]}


def parse_fasta(text, alphabet):
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    need(lines and lines[0].startswith(">"), "fasta_format", "Missing FASTA header")
    need(sum(line.startswith(">") for line in lines) == 1, "fasta_format", "Expected one FASTA record")
    sequence = "".join(lines[1:]).upper()
    need(sequence and set(sequence) <= alphabet, "fasta_alphabet", "Unexpected FASTA characters")
    return sequence


def feature(text, label):
    header = text.split("\nORIGIN", 1)[0]
    blocks = re.split(r"(?m)(?=^ {5}\S)", header)
    matches = [block for block in blocks if f'/label="{label}"' in block]
    need(len(matches) == 1, "feature_annotation", f"Expected one feature: {label}")
    block = matches[0]
    first = re.match(r" {5}(\S+)\s+(.+)\n", block)
    need(first is not None, "feature_annotation", f"Malformed feature: {label}")
    location = first[2]
    complement = location.startswith("complement(")
    if complement:
        location = location[len("complement("):-1]
    interval = re.fullmatch(r"(\d+)\.\.(\d+)", location)
    need(interval is not None, "feature_annotation", f"Unsupported location: {label}")
    protein = re.search(r'/translation="([^"]+)"', block)
    return {
        "type": first[1], "start": int(interval[1]), "end": int(interval[2]),
        "strand": -1 if complement else 1,
        "translation": re.sub(r"\s", "", protein[1]) if protein else None,
    }


def evaluate_candidate(sequence, context):
    """Use biological/sequence invariants, not equality to a constructed string."""
    found = {name: sites(sequence, values[0], circular=True) for name, values in ENZYMES.items()}
    need(all(len(value) == 1 for value in found.values()), "restriction_site_multiplicity",
         "Each cloning enzyme must cut once; internal sites are forbidden")
    hind, bam = found["HindIII"][0], found["BamHI"][0]
    need(hind < bam, "site_order", "Expected HindIII before BamHI at the reference origin")
    start = hind + 6 + len(context["kozak"])
    need(sequence[hind + 6:start] == context["kozak"], "initiation_context", "Incorrect initiation context")
    egfp_positions = sites(sequence, context["egfp_cds"])
    need(len(egfp_positions) == 1, "egfp_integrity", "Exactly one complete reference EGFP CDS is required")
    egfp_start = egfp_positions[0]
    end = egfp_start + len(context["egfp_cds"])
    need(start < bam < egfp_start, "fusion_orientation", "Incorrect TFEB-to-EGFP module order")
    insert_end = bam - len(context["adjustment"])
    insert = sequence[start:insert_end]
    need(insert != reverse_complement(context["tfeb_no_stop"]), "insert_orientation",
         "TFEB is reverse-complemented")
    protein_with_stop = translate(sequence[start:end])
    need(protein_with_stop.endswith("*"), "terminal_stop", "Missing final stop codon")
    need("*" not in protein_with_stop[:-1], "premature_stop", "Internal stop in fusion ORF")
    need(insert == context["tfeb_no_stop"], "tfeb_identity", "TFEB sequence or stop-removal boundary differs")
    linker = sequence[insert_end:egfp_start]
    need(linker == context["linker_nt"], "linker_identity", "Incorrect linker or frame adjustment")
    need(protein_with_stop[:-1] == context["expected_protein"], "fusion_translation",
         "Complete translated fusion differs from TFEB + linker + EGFP")
    return {"orf_start_1based": start + 1, "orf_end_1based": end,
            "tfeb_end_1based": insert_end, "egfp_start_1based": egfp_start + 1,
            "linker_start_1based": insert_end + 1, "linker_end_1based": egfp_start,
            "orf_nt_including_stop": end - start, "protein": protein_with_stop[:-1],
            "recognition_starts_1based": {key: value[0] + 1 for key, value in found.items()}}


def audit(report):
    def read(relative):
        path = ROOT / relative
        data = path.read_bytes()
        report["inputs"].append({"path": str(path.relative_to(ROOT)), "bytes": len(data),
                                  "file_sha256": digest(data)})
        return data.decode("utf-8")

    def passed(name, detail):
        report["checks"].append({"name": name, "passed": True, "detail": detail})

    spec = json.loads(read("设计/construct_spec.json"))
    report["design_version"] = spec["design_version"]
    need(spec["design_version"] == "B-v1", "unsupported_spec", "This audit is pinned to the reviewed B-v1 design")
    need(tuple(spec["tfeb"]["cds_1based"]) == TFEB_INTERVAL and
         tuple(spec["vector"]["egfp_cds_1based"]) == EGFP_INTERVAL,
         "reference_coordinates", "Spec differs from independently audited source coordinates")
    sources = {}
    for item in spec["reference_inputs"]:
        record = parse_genbank(read(item["path"]))
        accession = item["accession_version"]
        need(accession in PINNED_SEQUENCES and record["version"] == accession,
             "reference_accession", f"Unexpected accession: {accession}")
        sequence_hash = digest(record["sequence"].encode())
        need(sequence_hash == PINNED_SEQUENCES[accession] == item["sequence_sha256"],
             "reference_hash", f"Reference sequence hash differs: {accession}")
        need(report["inputs"][-1]["file_sha256"] == item["file_sha256"],
             "reference_file_hash", f"Source file changed since spec creation: {accession}")
        sources[accession] = record
    need(len(spec["reference_inputs"]) == len(PINNED_SEQUENCES) and
         set(sources) == set(PINNED_SEQUENCES), "reference_set", "Missing or duplicated required reference")
    vector = sources["U55762.1"]["sequence"]
    tfeb_cds = sources["NM_007162.3"]["sequence"][302:1733]
    tfeb_no_stop = tfeb_cds[:-3]
    tfeb_protein = sources["NP_009093.1"]["sequence"]
    egfp_cds = vector[678:1398]
    need(sources["U55762.1"]["circular"] and len(vector) == 4733,
         "vector_topology", "Expected the 4,733 bp circular source")
    need(sources["NM_007162.2"]["sequence"] == sources["NM_007162.3"]["sequence"] + "A" * 10,
         "transcript_versions", "Historical transcript differs beyond terminal poly(A)")
    need(translate(tfeb_cds) == tfeb_protein + "*" and tfeb_cds[-3:] == "TGA",
         "tfeb_reference_translation", "TFEB reference does not translate to NP_009093.1")
    need(translate(egfp_cds).endswith("*") and "*" not in translate(egfp_cds)[:-1]
         and egfp_cds[-3:] == "TAA", "egfp_reference_translation", "Invalid EGFP reference ORF")
    need(parse_fasta(read(spec["template"]["sequence_path"]), DNA_BASES) == tfeb_cds,
         "template_identity", "Planned synthetic template differs from the native reference CDS")
    passed("pinned_sources_and_independent_reference_translation",
           {"TFEB_CDS_nt": len(tfeb_cds), "TFEB_aa": len(tfeb_protein),
            "EGFP_aa": len(translate(egfp_cds)) - 1})

    need({item["name"] for item in spec["enzymes"]} == set(ENZYMES),
         "enzyme_set", "Unexpected cloning enzymes")
    for item in spec["enzymes"]:
        expected_enzyme = ENZYMES[item["name"]]
        observed_enzyme = (item["site"], item["top_cut_offset"], item["bottom_cut_offset"], item["overhang"])
        need(observed_enzyme == expected_enzyme and item["overhang_type"] == "5prime",
             "enzyme_geometry", f"Incorrect cut geometry: {item['name']}")
        need(not sites(tfeb_no_stop, item["site"]), "internal_template_site", "TFEB contains a cloning site")
    vector_sites = {name: sites(vector, values[0], circular=True) for name, values in ENZYMES.items()}
    need(all(len(value) == 1 for value in vector_sites.values()), "vector_sites", "Nonunique vector cloning site")
    h, b = vector_sites["HindIII"][0], vector_sites["BamHI"][0]
    need(h < b < 678, "mcs_order", "Unexpected source MCS site order")
    rows = list(csv.DictReader(io.StringIO(read("设计/primers.tsv")), delimiter="\t"))
    need(len(rows) == 2 and {row["direction"] for row in rows} == {"forward", "reverse"},
         "primer_rows", "Expected one forward and one reverse cloning primer")
    primers = {row["direction"]: row for row in rows}
    definitions = {item["direction"]: item for item in spec["primers"]}
    for direction, row in primers.items():
        definition = definitions[direction]
        anneal_length = definition["annealing_nt"]
        anneal = (tfeb_no_stop[:anneal_length] if direction == "forward" else
                  reverse_complement(tfeb_no_stop[-anneal_length:]))
        tail = definition["protective_bases"] + definition["restriction_site"] + definition["adjustment"]
        need(row["name"] == definition["name"] and row["sequence_5to3"] == tail + anneal
             and row["annealing_sequence_5to3"] == anneal,
             "primer_sequence", f"Primer or reverse-complement error: {direction}")
        for field in ("protective_bases", "restriction_site", "adjustment"):
            need(row[field] == definition[field], "primer_components", f"Incorrect primer component: {field}")
        expected_interval = (1, anneal_length, 1) if direction == "forward" else (
            len(tfeb_no_stop) - anneal_length + 1, len(tfeb_no_stop), -1)
        actual_interval = tuple(int(row[key]) for key in (
            "template_CDS_start_1based", "template_CDS_end_1based", "strand"))
        need(actual_interval == expected_interval, "primer_binding_coordinates", f"Incorrect binding interval: {direction}")
        need(int(row["full_length_nt"]) == len(tail + anneal) and int(row["annealing_length_nt"]) == len(anneal),
             "primer_lengths", f"Incorrect primer lengths: {direction}")
        for sequence, field in ((tail + anneal, "full_gc_percent"), (anneal, "annealing_gc_percent")):
            gc = 100 * sum(sequence.count(base) for base in "GC") / len(sequence)
            need(abs(float(row[field]) - gc) <= 0.00051, "primer_gc", f"GC mismatch: {direction}/{field}")
        terminal = anneal[-12:]
        need(len(sites(tfeb_cds, terminal)) + len(sites(reverse_complement(tfeb_cds), terminal)) == 1,
             "primer_terminal_specificity", f"Nonunique terminal 12mer in defined template: {direction}")
    passed("primers_sequences_tails_reverse_complement_and_binding",
           {direction: {"full_nt": len(row["sequence_5to3"]), "annealing_nt": len(row["annealing_sequence_5to3"])}
            for direction, row in primers.items()})

    kozak = definitions["forward"]["adjustment"]
    adjustment = reverse_complement(definitions["reverse"]["adjustment"])
    need(kozak == spec["fusion"]["kozak_added"] == "ACC" and
         adjustment == spec["fusion"]["frame_adjustment_coding_strand"] == "GG",
         "fusion_adjustments", "Incorrect initiation or reading-frame additions")
    linker_nt = adjustment + vector[b:678]
    linker_aa = translate(linker_nt)
    need(linker_nt == spec["fusion"]["expected_linker_nt"] and
         linker_aa == spec["fusion"]["expected_linker_aa"], "linker_design", "Spec linker differs from actual residual vector")
    direct_reconstruction = vector[:h + 6] + kozak + tfeb_no_stop + adjustment + vector[b:]
    f, r = primers["forward"], primers["reverse"]
    f_len, r_len = len(f["annealing_sequence_5to3"]), len(r["annealing_sequence_5to3"])
    pcr = f["sequence_5to3"] + tfeb_no_stop[f_len:-r_len] + reverse_complement(r["sequence_5to3"])
    pcr_sites = {name: sites(pcr, values[0]) for name, values in ENZYMES.items()}
    need(all(len(value) == 1 for value in pcr_sites.values()), "pcr_sites", "PCR must contain one site per cloning enzyme")
    ph, pb = pcr_sites["HindIII"][0], pcr_sites["BamHI"][0]
    insert_top = pcr[ph + 1:pb + 1]
    ligation_reconstruction = vector[:h + 1] + insert_top + vector[b + 1:]
    need(ligation_reconstruction == direct_reconstruction, "independent_reconstructions",
         "Primer-derived cut/ligation differs from direct reference-coordinate reconstruction")
    for name, offset in (("HindIII", ph), ("BamHI", pb)):
        need(pcr[offset + 1:offset + 5] == ENZYMES[name][3], "cohesive_end", f"Wrong 5-prime overhang: {name}")
    passed("explicit_PCR_top_bottom_cut_and_ligation",
           {"PCR_nt": len(pcr), "insert_nominal_nt": len(insert_top),
            "insert_ds_span_nt": (pb + 5) - (ph + 1),
            "backbone_nominal_nt": len(vector) - (b - h),
            "backbone_ds_span_nt": len(vector) - (b - h) + 4,
            "vector_top_cut_after_1based": [h + 1, b + 1],
            "PCR_top_cut_after_1based": [ph + 1, pb + 1],
            "junctions": {"upstream": "AAGCTT" + kozak + tfeb_no_stop[:12],
                          "downstream": tfeb_no_stop[-12:] + linker_nt + egfp_cds[:12]}})

    context = {"kozak": kozak, "adjustment": adjustment, "tfeb_no_stop": tfeb_no_stop,
               "egfp_cds": egfp_cds, "linker_nt": linker_nt,
               "expected_protein": tfeb_protein + linker_aa + translate(egfp_cds)[:-1]}
    fasta = parse_fasta(read("设计/construct.fasta"), DNA_BASES)
    genbank_text = read("设计/construct.gb")
    gb = parse_genbank(genbank_text)
    need(gb["circular"], "construct_topology", "Final GenBank must declare circular topology")
    need(fasta == gb["sequence"] == direct_reconstruction, "construct_full_sequence",
         "Final DNA differs from independently reconstructed, reference-origin sequence")
    evaluated = evaluate_candidate(fasta, context)
    protein = evaluated.pop("protein")
    need(parse_fasta(read("设计/fusion_protein.fasta"), AMINO_ACIDS) == protein,
         "protein_file", "Protein FASTA differs from independent full-ORF translation")
    fusion_feature = feature(genbank_text, "TFEB-EGFP fusion ORF")
    need((fusion_feature["start"], fusion_feature["end"], fusion_feature["strand"]) ==
         (evaluated["orf_start_1based"], evaluated["orf_end_1based"], 1)
         and fusion_feature["translation"] == protein,
         "fusion_annotation", "GenBank fusion CDS coordinates or translation are incorrect")
    delta = len(fasta) - len(vector)
    for label, first, last, strand, kind, _source in spec["source_feature_coordinates_1based"]:
        if last <= h + 1:
            mapped_first, mapped_last = first, last
        else:
            need(first > b + 1, "source_feature_overlap", f"Source feature crosses replaced MCS: {label}")
            mapped_first, mapped_last = first + delta, last + delta
        observed = feature(genbank_text, label)
        need(observed["type"] == kind and (observed["start"], observed["end"]) == (mapped_first, mapped_last),
             "backbone_feature_coordinates", f"Incorrect transformed feature: {label}")
        need(strand == 0 or observed["strand"] == strand,
             "backbone_feature_strand", f"Incorrect feature strand: {label}")
        need(fasta[mapped_first - 1:mapped_last] == vector[first - 1:last],
             "backbone_feature_sequence", f"Source feature sequence was altered: {label}")
    mass = math.fsum(FREE_AA_MASS[aa] for aa in protein) - (len(protein) - 1) * WATER_MASS
    metrics = json.loads(read("设计/construct_metrics.json"))
    expected_metrics = {
        "plasmid_bp": len(fasta), "empty_vector_bp": len(vector),
        "TFEB_reference_CDS_with_stop_nt": len(tfeb_cds), "TFEB_coding_no_stop_nt": len(tfeb_no_stop),
        "PCR_product_bp": len(pcr), "removed_vector_nominal_nt": b - h,
        "digested_backbone_nominal_nt": len(vector) - (b - h), "digested_insert_nominal_nt": len(insert_top),
        "TFEB_aa": len(tfeb_protein), "linker_aa": len(linker_aa), "EGFP_aa": len(translate(egfp_cds)) - 1,
        "fusion_aa": len(protein), "fusion_orf_nt_including_stop": evaluated["orf_nt_including_stop"],
        "fusion_start_1based": evaluated["orf_start_1based"], "fusion_end_1based": evaluated["orf_end_1based"],
        "TFEB_end_1based": evaluated["tfeb_end_1based"], "EGFP_start_1based": evaluated["egfp_start_1based"],
        "linker_sequence": linker_nt, "linker_peptide": linker_aa,
        "vector_top_cut_after_1based": [h + 1, b + 1], "PCR_top_cut_after_1based": [ph + 1, pb + 1],
    }
    for key, expected_value in expected_metrics.items():
        need(metrics[key] == expected_value, "metric_mismatch", f"Main-line metric differs: {key}")
    for name, nominal in (("insert", len(insert_top)), ("backbone", len(vector) - (b - h))):
        need(metrics[name]["ds_span_bp"] == nominal + 4 and metrics[name]["watson_nt"] == nominal
             and metrics[name]["crick_nt"] == nominal, "duplex_span_metric", f"Incorrect nominal/ds-span distinction: {name}")
        left, right = ("agct", "gatc") if name == "insert" else ("gatc", "agct")
        need(metrics[name]["five_prime_end"] == ["5'", left] and
             metrics[name]["three_prime_end"] == ["5'", right],
             "cohesive_end_metric", f"Reported fragment overhangs are incorrect: {name}")
    for enzyme, motif in (("HindIII", "AAGCTT"), ("BamHI", "GGATCC"),
                          ("XhoI", "CTCGAG"), ("KpnI", "GGTACC")):
        cut_offset = 5 if enzyme == "KpnI" else 1
        cut_positions = [position + cut_offset for position in sites(fasta, motif, circular=True)]
        need(metrics["construct_top_cut_after_1based"][enzyme] == cut_positions,
             "restriction_metric", f"Incorrect final restriction map positions: {enzyme}")
    need(abs(metrics["fusion_molecular_weight_Da"] - mass) < 0.01,
         "molecular_weight", "Main-line molecular weight differs by at least 0.01 Da")
    report["verified_construct"] = {**evaluated, "plasmid_bp": len(fasta), "fusion_aa": len(protein),
                                    "fusion_molecular_weight_Da": round(mass, 6),
                                    "mass_difference_Da": mass - metrics["fusion_molecular_weight_Da"],
                                    "DNA_sha256": digest(fasta.encode()), "protein_sha256": digest(protein.encode()),
                                    "linker_nt": linker_nt, "linker_aa": linker_aa}
    passed("complete_construct_ORF_protein_mass_backbone_annotations_and_metrics", report["verified_construct"])

    start = evaluated["orf_start_1based"] - 1
    end = start + len(tfeb_no_stop)
    mutations = [
        ("reverse_TFEB_insert", fasta[:start] + reverse_complement(tfeb_no_stop) + fasta[end:], "insert_orientation"),
        ("one_base_frameshift", fasta[:start + 90] + "A" + fasta[start + 90:], "reading_frame"),
        ("restore_native_TFEB_stop", fasta[:end] + "TGA" + fasta[end:], "premature_stop"),
        ("introduce_internal_HindIII", fasta[:start + 300] + "AAGCTT" + fasta[start + 306:], "restriction_site_multiplicity"),
    ]
    for name, mutated, expected_code in mutations:
        try:
            evaluate_candidate(mutated, context)
        except AuditFailure as error:
            rejected = error.code == expected_code
            report["negative_controls"].append({"name": name, "rejected": True,
                                                "expected_reason": expected_code, "actual_reason": error.code,
                                                "reason_matches": rejected, "detail": str(error),
                                                "mutant_sequence_sha256": digest(mutated.encode())})
            need(rejected, "negative_control_wrong_reason", f"{name} failed for an unexpected reason: {error.code}")
        else:
            report["negative_controls"].append({"name": name, "rejected": False,
                                                "expected_reason": expected_code})
            raise AuditFailure("negative_control_accepted", f"Semantic validation accepted {name}")
    passed("four_semantic_negative_controls", "All rejected by their intended invariant; no exact-sequence equality check used")
    for entry in report["inputs"]:
        need(digest((ROOT / entry["path"]).read_bytes()) == entry["file_sha256"],
             "concurrent_input_change", f"Input changed during audit: {entry['path']}")
    passed("input_snapshot_unchanged", "Every recorded input hash was rechecked after validation")


def main():
    report = {
        "status": "running", "success": False,
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "actual_command": shlex.join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]),
        "python_version": sys.version.split()[0],
        "script_path": str(Path(__file__).resolve().relative_to(ROOT)),
        "script_sha256": digest(Path(__file__).read_bytes()),
        "implementation": "Python standard library only; no builder imports, pydna, Bio, or primer3",
        "independence": "Pinned source coordinates and sequence hashes; literal cut geometry; independent codon lookup and mass summation",
        "mass_convention": "Average neutral unmodified linear protein; free amino acid masses minus (n-1)*18.0153 Da; same convention as main-line average MW",
        "inputs": [], "checks": [], "negative_controls": [],
        "limitations": [
            "No wet-lab validation or verification of acquired physical stocks",
            "No whole-transcriptome off-target analysis; terminal 12mers checked only in the defined synthetic template",
            "Primer3 thermodynamics are not independently recomputed by this standard-library audit",
            "The average-mass convention is shared with the main calculation; algorithm and constants are independently embedded",
            "Reference-origin output coordinates are required; arbitrary circular rotations are not accepted",
            "Backbone feature coordinate mapping is checked against the approved spec, not independently reannotated from functional assays",
        ],
        "self_review_scope": [
            "One-based inclusive versus zero-based half-open boundaries",
            "Reverse-primer annealing orientation and frame-adjustment reverse complement",
            "Five-prime sticky ends, equal strand lengths and four-base larger duplex spans",
            "Full-ORF terminal stop, annotation mapping and mutation rejection before exact equality",
            "Only the designated verification JSON is written; failed checks return a nonzero exit code",
        ],
    }
    try:
        need(len(sys.argv) == 1, "arguments", "Run without arguments")
        audit(report)
        report.update(status="pass", success=True)
        exit_code = 0
    except Exception as error:
        report.update(status="fail", error_code=getattr(error, "code", type(error).__name__), error=str(error))
        exit_code = 1
    OUTPUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "output": str(OUTPUT),
                      "passed_checks": len(report["checks"]),
                      "negative_controls": len(report["negative_controls"]),
                      "error": report.get("error")}, ensure_ascii=False))
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
