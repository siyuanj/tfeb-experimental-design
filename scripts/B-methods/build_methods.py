#!/usr/bin/env python3
"""Generate P6 diagnostics only after verifying the P5 ready manifest."""
from __future__ import annotations

import csv
import importlib.metadata
import json
import math
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Restriction
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import primer3

from sequence_checks import (binding_sites, choose_primer, digest_literal,
                             occurrences, pcr_products, rc, sha256)

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "设计"
C = ROOT / "docs/协作"


def load_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def write_tsv(path, records):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def fasta(path):
    return str(SeqIO.read(path, "fasta").seq).upper()


def input_gate():
    ready_path = C / "sequence-ready.json"
    if not ready_path.exists():
        raise RuntimeError("P5 sequence-ready.json is absent; no P6 output generated")
    ready = load_json(ready_path)
    if ready.get("status") != "ready":
        raise RuntimeError("P5 has not passed its readiness gate")
    # The manifest schema is explicit; unknown schemas must not silently pass.
    hashes = ready.get("files", ready.get("artifacts", ready.get("outputs")))
    if isinstance(hashes, list):
        hashes = {item["path"]: item.get("sha256", item.get("file_sha256")) for item in hashes}
    if not isinstance(hashes, dict) or not hashes:
        raise RuntimeError("Readiness manifest has no supported file hash map")
    checked = {}
    for relative, expected in hashes.items():
        if isinstance(expected, dict):
            expected = expected.get("sha256", expected.get("file_sha256"))
        actual = sha256(ROOT / relative)
        if actual != expected:
            raise RuntimeError(f"P5 manifest hash mismatch: {relative}")
        checked[relative] = actual
    for required in ("设计/construct_spec.json", "设计/construct.fasta", "设计/construct.gb",
                     "设计/primers.tsv", "设计/verification_summary.json",
                     "设计/PCR_product.fasta", "设计/fusion_protein.fasta"):
        if required not in checked:
            raise RuntimeError(f"Required input is not bound by P5 manifest: {required}")
    checked["docs/协作/sequence-ready.json"] = sha256(ready_path)
    return ready, checked


def main():
    old_ready = C / "methods-ready.json"
    if old_ready.exists():
        prior_hash = sha256(old_ready)
        write_json(old_ready, {"status":"revalidation_required",
                              "reason":"P6 generator invoked;recheck all regenerated outputs before ready",
                              "previous_ready_sha256":prior_hash})
    ready, input_hashes = input_gate()
    spec = load_json(D / "construct_spec.json")
    for item in spec["reference_inputs"]:
        actual = sha256(ROOT / item["path"])
        if actual != item["file_sha256"]:
            raise RuntimeError(f"Source hash changed: {item['path']}")
        input_hashes[item["path"]] = actual
    sequence = fasta(D / "construct.fasta")
    gb = SeqIO.read(D / "construct.gb", "genbank")
    vector_record = SeqIO.read(ROOT / spec["vector"]["path"], "genbank")
    vector = str(vector_record.seq).upper()
    reference = SeqIO.read(ROOT / "资料/序列/TFEB_reference.gb", "genbank")
    cds_feature = next(f for f in reference.features if f.type == "CDS")
    cds = str(cds_feature.extract(reference.seq)).upper()
    tfeb = cds[:-3]
    tfeb_start = sequence.index(tfeb)
    tfeb_end = tfeb_start + len(tfeb)
    egfp_feature = next(f for f in vector_record.features
                       if f.type == "CDS" and f.qualifiers.get("gene") == ["egfp"])
    egfp = str(egfp_feature.extract(vector_record.seq)).upper()
    egfp_start = sequence.index(egfp)
    egfp_end = egfp_start + len(egfp)
    orf = sequence[tfeb_start:egfp_end]
    linker = sequence[tfeb_end:egfp_start]
    checks = {}
    checks["gb_fasta_same_sequence_and_circular"] = str(gb.seq).upper() == sequence and gb.annotations.get("topology") == "circular"
    checks["reference_translation_and_stop"] = str(Seq(cds).translate()) == cds_feature.qualifiers["translation"][0] + "*"
    checks["TFEB_and_EGFP_exactly_once"] = len(occurrences(sequence, tfeb)) == len(occurrences(sequence, egfp)) == 1
    protein = str(Seq(orf).translate())
    def fusion_is_valid(s):
        if s.count(tfeb) != 1 or s.count(egfp) != 1:
            return False
        coding = s[s.index(tfeb):s.index(egfp) + len(egfp)]
        if len(coding) % 3:
            return False
        translated = str(Seq(coding).translate())
        expected = cds_feature.qualifiers["translation"][0] + spec["fusion"]["expected_linker_aa"] + egfp_feature.qualifiers["translation"][0] + "*"
        return translated.endswith("*") and "*" not in translated[:-1] and translated == expected
    checks["no_internal_stop_and_expected_fusion"] = fusion_is_valid(sequence)
    checks["linker_matches_spec"] = linker == spec["fusion"]["expected_linker_nt"]
    checks["fusion_output_matches_independent_translation"] = protein[:-1] == fasta(D / "fusion_protein.fasta")
    with (D / "primers.tsv").open() as stream:
        primers = list(csv.DictReader(stream, delimiter="\t"))
    fprimer = next(p for p in primers if p["direction"] == "forward")
    rprimer = next(p for p in primers if p["direction"] == "reverse")
    fa, ra = fprimer["annealing_sequence_5to3"], rprimer["annealing_sequence_5to3"]
    checks["full_primers_match_declared_segments"] = all(
        p["sequence_5to3"] == p["protective_bases"] + p["restriction_site"] + p["adjustment"] + p["annealing_sequence_5to3"]
        and p["sequence_5to3"].endswith(p["annealing_sequence_5to3"]) for p in primers)
    checks["forward_primer_exact_start"] = fa == tfeb[:len(fa)]
    checks["reverse_primer_rc_exact_end_no_stop"] = rc(ra) == tfeb[-len(ra):]
    ftail = fprimer["sequence_5to3"][:-len(fa)]
    rtail = rprimer["sequence_5to3"][:-len(ra)]
    pcr = ftail + tfeb + rc(rtail)
    checks["tailed_PCR_independent_reconstruction"] = pcr == fasta(D / "PCR_product.fasta")
    enzymes = [(e["name"], e["site"], e["top_cut_offset"]) for e in spec["enzymes"]]
    # This protocol is intentionally tied to the accepted pair, not arbitrary enzymes.
    checks["selected_pair_matches_method_sources"] = [e[0] for e in enzymes] == ["HindIII", "BamHI"]
    vd, pd, cd = [digest_literal(s, enzymes, circular) for s, circular in
                  ((vector, True), (pcr, False), (sequence, True))]
    vh, vb = [e["cut_top_0based_boundary"] for e in vd["cuts"]]
    ph, pb = [e["cut_top_0based_boundary"] for e in pd["cuts"]]
    reconstructed = vector[:vh] + pcr[ph:pb] + vector[vb:]
    checks["independent_literal_construct_reconstruction"] = reconstructed == sequence
    checks["no_internal_HindIII_or_BamHI_in_TFEB"] = all(not occurrences(tfeb, e[1]) for e in enzymes)
    checks["literal_cuts_match_Biopython"] = all(
        [x["cut_top_0based_boundary"] + 1 for x in cd["cuts"] if x["enzyme"] == name]
        == getattr(Restriction, name).search(Seq(sequence), linear=False)
        for name, _, _ in enzymes)
    # Recognition-site and frame checks must detect intentional defective models.
    reverse_model = sequence[:tfeb_start] + rc(tfeb) + sequence[tfeb_end:]
    stop_model = sequence[:tfeb_end] + cds[-3:] + sequence[tfeb_end:]
    frame_model = sequence[:tfeb_end] + "A" + sequence[tfeb_end:]
    checks["reverse_model_rejected_by_fusion_validator"] = not fusion_is_valid(reverse_model)
    checks["native_stop_retention_model_rejected"] = not fusion_is_valid(stop_model)
    checks["frameshift_model_rejected_by_fusion_validator"] = not fusion_is_valid(frame_model)
    if not all(checks.values()):
        write_json(C / "方法独立核验.json", {"status": "failed", "checks": checks, "inputs": input_hashes})
        raise RuntimeError("Independent sequence checks failed")

    # Include 20 upstream vector bases, the full junction, and the terminal stop.
    target_start = next(e["site_start_1based"] - 21 for e in cd["cuts"] if e["enzyme"] == "HindIII")
    target_end = egfp_end
    sanger = []
    for strand, prefix in ((1, "S-F"), (-1, "S-R")):
        needed = target_start if strand == 1 else target_end - 1
        for i in range(10):
            center = needed - 90 if strand == 1 else needed + 100
            window = (needed - 200, needed - 60) if strand == 1 else (needed + 60, needed + 200)
            p = choose_primer(sequence, window, strand, center, f"{prefix}{i+1}", must_cover_0=needed)
            sanger.append(p)
            if strand == 1:
                if p["usable_end_1based"] >= target_end:
                    break
                needed = p["usable_end_1based"] - 50
            else:
                if p["usable_start_1based"] <= target_start + 1:
                    break
                needed = p["usable_start_1based"] + 48
        else:
            raise RuntimeError("Primer walking did not finish in ten reads per strand")
    coverage = []
    for p in sanger:
        a = max(target_start + 1, p["usable_start_1based"])
        b = min(target_end, p["usable_end_1based"])
        coverage.append({"primer_id": p["primer_id"], "direction": p["direction"],
                         "binding_start_1based": p["binding_start_1based"],
                         "binding_end_1based": p["binding_end_1based"],
                         "usable_start_1based": p["usable_start_1based"],
                         "usable_end_1based": p["usable_end_1based"],
                         "target_overlap_start_1based": a, "target_overlap_end_1based": b,
                         "target_overlap_nt": max(0, b - a + 1),
                         "planning_assumption": "extended offsets51-650;600nt;not observed reads"})
    strand_depth = {strand: [sum(p["direction"] == strand and p["usable_start_1based"] <= x <= p["usable_end_1based"]
                                for p in sanger) for x in range(target_start + 1, target_end + 1)]
                    for strand in ("+", "-")}
    checks["full_target_bidirectional_coverage"] = all(min(depth) >= 1 for depth in strand_depth.values())
    checks["all_Sanger_primers_unique_exact_binding"] = all(len(binding_sites(sequence, p["sequence_5to3"])) == 1 for p in sanger)
    if not all(checks.values()):
        raise RuntimeError("Sanger coverage did not pass")
    by_name = {p["primer_id"]: p for p in sanger}
    # Reuse a vector pair and an insert/vector pair without extra oligo orders.
    internal_forward = next(p for p in sanger if p["direction"] == "+"
                            and tfeb_start + 1 <= p["binding_start_1based"]
                            and p["binding_end_1based"] <= tfeb_end)
    assays = {"colony_flanking": (by_name["S-F1"], by_name["S-R1"]),
              "colony_orientation": (internal_forward, by_name["S-R1"])}
    fragments = []
    def add(assay, sample, reagents, lengths, meaning, limitation="", definition="nominal bp"):
        for i, size in enumerate(lengths or [0], 1):
            fragments.append({"assay": assay, "sample": sample, "reagents": reagents,
                              "fragment_id": i, "size_bp": size,
                              "length_definition": definition,
                              "visibility": "no specific band predicted" if size == 0 else ("below routine1%gel useful range" if size < 100 else "predicted;not observed"),
                              "interpretation": meaning, "limitation": limitation})
    add("marker", "NEB N3232", "1kb ladder", [10002,8001,6001,5001,4001,3001,2000,1500,1000,517,500],
        "manufacturer listed sizes", "500 and517 co-migrate on routine1%agarose")
    add("cloning_PCR", "synthetic_TFEB_template", "TFEB-HindIII-F + TFEB-BamHI-R", [len(pcr)], "tailed full-length PCR target")
    add("cloning_PCR", "no_template_control", "TFEB-HindIII-F + TFEB-BamHI-R", [], "no specific product expected", "contamination and primer dimers are empirical possibilities")
    add("preparative_digest", "empty_vector", "HindIII-HF + BamHI-HF", vd["fragments_bp"], "recover larger backbone", "small MCS may be invisible", "top-strand cut-to-cut nt")
    add("preparative_digest", "PCR_product", "HindIII-HF + BamHI-HF", pd["fragments_bp"], "recover central insert", "7/11nt terminal pieces not visible;strand lengths differ only for discarded terminals", "top-strand cut-to-cut nt")
    for sample, sample_seq in (("correct_construct", sequence), ("empty_vector", vector),
                               ("reverse_TFEB_model", reverse_model), ("retained_stop_model", stop_model)):
        add("diagnostic_digest", sample, "HindIII-HF + BamHI-HF",
            digest_literal(sample_seq, enzymes)["fragments_bp"], "gross insert/backbone size check",
            "reverse model can match correct;3bp stop difference not resolvable;Sanger required", "top-strand cut-to-cut nt")
        for assay, pair in assays.items():
            products = pcr_products(sample_seq, *(p["sequence_5to3"] for p in pair))
            add(assay, sample, " + ".join(p["primer_id"] for p in pair), [p["length_bp"] for p in products],
                "insert-size screen" if assay == "colony_flanking" else "orientation screen",
                "no product requires positive PCR control;3bp difference not resolvable")
    for assay, pair in assays.items():
        add(assay, "no_template_control", " + ".join(p["primer_id"] for p in pair), [], "no specific product expected")
    # Prespecified controls must actually discriminate the modeled failures.
    flanking_rows = [r for r in fragments if r["assay"] == "colony_flanking"]
    sizes = {r["sample"]: r["size_bp"] for r in flanking_rows}
    checks["flanking_distinguishes_empty_vector"] = sizes["correct_construct"] != sizes["empty_vector"] > 0
    orientation_rows = {r["sample"]: r["size_bp"] for r in fragments if r["assay"] == "colony_orientation"}
    checks["orientation_pair_discriminates_reverse_model"] = orientation_rows["correct_construct"] > 0 and orientation_rows["reverse_TFEB_model"] == orientation_rows["empty_vector"] == 0
    # For this accepted pair both termini are four-base 5-prime overhangs.
    backbone_bp, insert_bp = max(vd["fragments_bp"]), pb - ph
    checks["equal_strand_lengths_for_molarity"] = all(e["bottom_cut_offset"] - e["top_cut_offset"] == 4 for e in spec["enzymes"])
    vector_ng, insert_ng = 20 * backbone_bp * .00066, 60 * insert_bp * .00066
    total_ng_per_uL = (vector_ng + insert_ng) / 20
    # 20 and10 ng/uL stocks are design dilution targets, not observed measurements.
    vector_uL, insert_uL = vector_ng / 20, insert_ng / 10
    water_uL = 20 - 2 - 1 - vector_uL - insert_uL
    checks["ligation_concentration_and_volume"] = 1 <= total_ng_per_uL <= 10 and water_uL > 0
    pcr_extension = math.ceil(len(pcr) / 1000 * 30 / 5) * 5
    colony_longest = max(r["size_bp"] for r in fragments if r["assay"].startswith("colony"))
    colony_extension = math.ceil(colony_longest / 1000 * 60 / 5) * 5
    pair_details = {}
    for name, pair in assays.items():
        tm = [primer3.calc_tm(p["sequence_5to3"], mv_conc=50, dv_conc=1.8,
                             dntp_conc=.8, dna_conc=200) for p in pair]
        pair_details[name] = {"primers": [p["primer_id"] for p in pair],
                              "model_tm_C": tm, "initial_Ta_C": 58,
                              "optimization_gradient_C": [55,58,61],
                              "heterodimer_tm_C": primer3.calc_heterodimer(*(p["sequence_5to3"] for p in pair),
                                  mv_conc=50, dv_conc=1.8, dntp_conc=.8, dna_conc=200).tm}
    checks["colony_primer_pairs_thermodynamic_screen"] = all(abs(v["model_tm_C"][0] - v["model_tm_C"][1]) <= 5 and v["heterodimer_tm_C"] < 45 for v in pair_details.values())
    parameters = {
        "status": "designed_not_experimentally_validated", "design_version": spec["design_version"],
        "inputs": input_hashes, "coordinate_convention": "1-based inclusive tables;0-based half-open implementation",
        "source_ids": [f"M{i:02}" for i in [1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]],
        "cloning_PCR": {"product": "NEB M0492", "total_uL":50,"mix_2X_uL":25,
                        "primer_10uM_each_uL":2.5,"template_target_ng":1,"template_target_ng_per_uL":1,
                        "template_uL":1,"water_uL":19,"initial_denaturation_C_s":[98,30],"cycles":30,
                        "denaturation_C_s":[98,10],"annealing_C_s":[spec["primer_thermodynamics"]["initial_PCR_Ta_C"],20],
                        "extension_C_s":[72,pcr_extension],"final_extension_C_s":[72,120],"gradient_C":[62,65,68],
                        "template_type":spec["template"]["type"],"wet_lab_optimized":False},
        "digest": {"enzymes":["HindIII-HF R3104","BamHI-HF R3136"],"each_stock_U_per_uL":20,
                   "each_uL":1,"DNA_target_ng":500,"total_uL":50,"buffer10X_uL":5,"C":37,"min":60,
                   "DNA_target_ng_per_uL":50,"DNA_uL":10,"water_uL":33,"glycerol_percent":2,
                   "cleanup_required":True,"heat_inactivation_used":False},
        "digest_gel_loading": {"loading_dye":"NEB B7024 Purple Gel Loading Dye, with SDS",
                               "stock_X":6,"final_X":1,"dye_to_sample_volume_ratio":[1,5],
                               "scope":"Direct preparative and diagnostic gel loading of unpurified BamHI-HF digests only; not added to ligation",
                               "SYBR_Safe_or_GelRed":"post-electrophoresis staining only","source_ids":["M19","M20"]},
        "ligation": {"length_for_molarity":"mean strand length;equals top-strand nominal length for selected equal5prime overhangs",
                     "backbone_bp":backbone_bp,"insert_bp":insert_bp,"vector_fmol":20,"insert_fmol":60,
                     "insert_to_vector_molar_ratio":3,"vector_ng":vector_ng,"insert_ng":insert_ng,
                     "vector_stock_target_ng_per_uL":20,"insert_stock_target_ng_per_uL":10,
                     "vector_uL":vector_uL,"insert_uL":insert_uL,"buffer10X_uL":2,"ligase_uL":1,
                     "ligase_U":400,"water_uL":water_uL,"total_uL":20,"total_DNA_ng_per_uL":total_ng_per_uL,
                     "C_range":[20,25],"min":10,"heat_inactivation_C_min":[65,10]},
        "transformation": {"cells":"NEB5-alpha C2987 DH5alpha derivative","cells_uL":50,"ligation_uL":2,
                           "ice_before_min":30,"heat_shock_C_s":[42,30],"ice_after_min":5,"SOC_uL":950,
                           "outgrowth_C_min_rpm":[37,60,250],"plate_uL":100,"dilutions":[1,10,100],
                           "kanamycin_ug_per_mL":30,"pUC19_control_antibiotic":"ampicillin100ug/mL"},
        "colony_PCR": {"product":"OneTaq M0482","colonies_initial_screen":8,"total_uL":25,
                       "mix2X_uL":12.5,"primer10uM_each_uL":.5,"water_uL":11.5,
                       "template":"trace colony material;patch saved","initial_lysis_C_s":[94,300],"cycles":30,
                       "denaturation_C_s":[94,30],"annealing_C_s":[58,30],"extension_C_s":[68,colony_extension],
                       "final_extension_C_s":[68,300],"pairs":pair_details,
                       "Tm_model":{"mv_mM":50,"Mg_mM":1.8,"total_dNTP_mM":.8,"primer_nM":200,"salt":"SantaLucia;monovalent assumption"}},
        "sanger": {"service":"Eurofins Standard Read purified plasmid","primer_count":len(sanger),
                   "target_start_1based":target_start+1,"target_end_1based":target_end,
                   "target_nt":target_end-target_start,"usable_offsets_from_3prime":[51,650],
                   "usable_nt_per_read":600,"minimum_forward_depth":min(strand_depth["+"]),
                   "minimum_reverse_depth":min(strand_depth["-"]),"quality_threshold":"Q20 both directions;no unresolved ambiguity",
                   "Tm_model":{"mv_mM":50,"Mg_mM":0,"dNTP_mM":0,"primer_nM":50,"model":"SantaLucia","actual_service_chemistry":False},
                   "design_constraints":{"length_nt":[18,25],"Tm_C":[50,62],"GC_percent":[35,60],"hairpin_or_homodimer_Tm_max_C":40,
                                         "homopolymer_max":3,"unique_full_match_in_final_plasmid":True}},
        "limitations":["All outputs are predictions;no wet-lab data","Exact-match specificity only within supplied constructs",
                       "Reverse model is a hypothetical failure challenge,not a product of compatible directional ligation",
                       "Sanger windows require validation from actual quality-trimmed traces;whole backbone is not fully covered"]}
    checks["PCR_volume_50uL"] = 25+2.5+2.5+1+19 == 50
    checks["colony_PCR_volume_25uL"] = 12.5+.5+.5+11.5 == 25
    checks["digest_volume_50uL"] = 5+1+1+10+33 == 50
    checks["digestion_fragments_sum_to_input"] = sum(vd["fragments_bp"]) == len(vector) and sum(pd["fragments_bp"]) == len(pcr) and sum(cd["fragments_bp"]) == len(sequence)
    if not all(checks.values()):
        write_json(C / "方法独立核验.json", {"status":"failed","checks":checks,"inputs":input_hashes})
        raise RuntimeError("Method validation failed")
    _, current_hashes = input_gate()
    if current_hashes != input_hashes:
        raise RuntimeError("P5 inputs changed during method calculation;retry on a stable snapshot")
    for item in spec["reference_inputs"]:
        if sha256(ROOT / item["path"]) != item["file_sha256"]:
            raise RuntimeError("Source changed during method calculation")
    write_tsv(D / "expected_fragments.tsv", fragments)
    write_tsv(D / "sanger_primers.tsv", sanger)
    write_tsv(D / "sanger_coverage.tsv", coverage)
    write_json(D / "methods_parameters.json", parameters)
    outputs = {str(p.relative_to(ROOT)):sha256(p) for p in [D / n for n in
               ("expected_fragments.tsv","sanger_primers.tsv","sanger_coverage.tsv","methods_parameters.json")]}
    result = {"status":"passed","checked_at":datetime.now(ZoneInfo("America/New_York")).isoformat(),
              "design_version":spec["design_version"],"checks":checks,"inputs":input_hashes,"outputs":outputs,
              "independent_path":"literal motif search,manual top-cut offsets,reference feature extraction;no B-sequence imports or pydna assembly",
              "metrics":{"construct_bp":len(sequence),"PCR_bp":len(pcr),"digest_backbone_bp":backbone_bp,
                         "digest_insert_bp":insert_bp,"fusion_aa":len(protein)-1,
                         "fusion_MW_Da":ProteinAnalysis(protein[:-1]).molecular_weight(),
                         "sanger_target_1based":[target_start+1,target_end],"sanger_target_nt":target_end-target_start,
                         "Sanger_primers":len(sanger)},
              "software":{"python":sys.version.split()[0],"biopython":importlib.metadata.version("biopython"),
                          "primer3-py":importlib.metadata.version("primer3-py")},
              "scripts":{str(p.relative_to(ROOT)):sha256(p) for p in Path(__file__).parent.glob("*.py")},
              "command":"'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/B-methods/build_methods.py"}
    write_json(C / "方法独立核验.json", result)
    print(json.dumps({"status":result["status"],"checks":len(checks),"metrics":result["metrics"]},ensure_ascii=False,indent=2))


if __name__ == "__main__":
    main()
