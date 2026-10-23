#!/usr/bin/env python3
"""Freeze reviewed primary references and explicit design choices for B-v1."""
import hashlib
import json
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[2]

def require(condition, message):
    if not condition:
        raise ValueError(message)

def main():
    source = ROOT / "资料/序列"
    out = ROOT / "设计"
    out.mkdir(exist_ok=True)
    current = SeqIO.read(source / "TFEB_reference.gb", "genbank")
    old = SeqIO.read(source / "TFEB_NM_007162.2.gb", "genbank")
    protein = SeqIO.read(source / "TFEB_NP_009093.1.gp", "genbank")
    vector = SeqIO.read(source / "pEGFP-N1_U55762.1.gb", "genbank")
    f = next(f for f in current.features if f.type == "CDS" and f.qualifiers.get("gene") == ["TFEB"])
    old_f = next(f for f in old.features if f.type == "CDS" and f.qualifiers.get("gene") == ["TFEB"])
    cds = f.extract(current.seq).upper()
    require(current.id == "NM_007162.3" and protein.id == "NP_009093.1" and vector.id == "U55762.1", "Unexpected accession")
    require(cds.translate(cds=True) == protein.seq and cds == old_f.extract(old.seq), "CDS/protein mismatch")
    require(old.seq == current.seq + "A" * 10, "Unexpected transcript difference")
    require(len(cds) == 1431 and len(vector) == 4733, "Unexpected source length")
    for filename, sequence, identifier, description in [
        ("TFEB_cds.fasta", cds, current.id + "_CDS", "303..1733; includes terminal TGA"),
        ("TFEB_cds_no_stop.fasta", cds[:-3], current.id + "_coding", "1428 nt; no cloning tails"),
        ("TFEB_protein.fasta", protein.seq, protein.id, "476 aa; native TFEB"),
        ("pEGFP-N1_U55762.1.fasta", vector.seq.upper(), vector.id, "4733 bp; circular empty vector")]:
        SeqIO.write(SeqRecord(sequence, id=identifier, description=description), source / filename, "fasta")
    inputs = json.loads((source / "source_manifest.json").read_text())
    spec = {
        "design_version": "B-v1", "construct_name": "pEGFP-N1-TFEB-Bv1",
        "reference_inputs": inputs,
        "tfeb": {"transcript": "NM_007162.3", "protein": "NP_009093.1", "cds_1based": [303,1733], "strand": 1,
                 "reference_cds_nt_including_stop": 1431, "coding_nt_without_stop": 1428, "protein_aa": 476,
                 "native_stop": "TGA", "native_start_retained": True, "native_stop_retained": False},
        "template": {"type": "planned synthetic linear dsDNA", "name": "TFEB native CDS NM_007162.3",
                     "sequence_path": "资料/序列/TFEB_cds.fasta", "nt": 1431,
                     "procurement_plan": "Plan custom native CDS dsDNA synthesis without codon optimization, e.g. IDT gBlocks HiFi (sequence-verified clonal intermediate). Not acquired; vendor feasibility/QC and physical template identity remain to be checked, and final clones require complete Sanger verification.",
                     "rationale": "A defined intron-free isoform-1 template avoids uncertain HeLa transcript isoforms; PCR/restriction/ligation remain the construction route."},
        "vector": {"accession": "U55762.1", "path": "资料/序列/pEGFP-N1_U55762.1.gb", "length_bp":4733,
                   "topology":"circular", "physical_stock": "planned empty pEGFP-N1 stock; verify identity before use",
                   "egfp_cds_1based":[679,1398], "egfp_aa":239, "egfp_stop":"TAA"},
        "fusion": {"orientation":"TFEB-EGFP", "tag_location":"C-terminal to TFEB", "upstream_initiation_context":"AAGCTTACCATGG",
                   "kozak_added":"ACC", "kozak_rationale":"A at -3 and G at +4; shorter context avoids the predicted 66 C hairpin of GCCACC + native TFEB start.",
                   "frame_adjustment_coding_strand":"GG", "frame_adjustment_reverse_primer":"CC",
                   "expected_linker_nt":"GGGGATCCACCGGTCGCCACC", "expected_linker_aa":"GDPPVAT",
                   "egfp_ATG_retained_as_internal_methionine":True, "terminal_stop":"TAA",
                   "not_identical_to_addgene_38119":True},
        "enzymes": [
            {"name":"HindIII", "product":"HindIII-HF R3104", "site":"AAGCTT", "top_cut_offset":1, "bottom_cut_offset":5, "overhang":"AGCT", "overhang_type":"5prime"},
            {"name":"BamHI", "product":"BamHI-HF R3136", "site":"GGATCC", "top_cut_offset":1, "bottom_cut_offset":5, "overhang":"GATC", "overhang_type":"5prime"}],
        "primers": [
            {"name":"TFEB-HindIII-F", "direction":"forward", "protective_bases":"ACTGAC", "restriction_site":"AAGCTT", "adjustment":"ACC", "adjustment_role":"translation initiation context", "annealing_nt":22},
            {"name":"TFEB-BamHI-R", "direction":"reverse", "protective_bases":"TATGAC", "restriction_site":"GGATCC", "adjustment":"CC", "adjustment_role":"reverse complement of coding-strand GG frame adjustment", "annealing_nt":22}],
        "primer_thermodynamics": {"mv_conc_mM":50.0,"dv_conc_mM":2.0,"dntp_conc_mM_total":0.8,"dna_conc_nM":500.0,
                                  "structure_temperature_C":65.0,"tm_method":"santalucia","salt_corrections_method":"santalucia",
                                  "scope":"Primer3 explicit model; 50 mM monovalent salt is a modeling assumption, not Q5 proprietary buffer composition",
                                  "acceptance":{"annealing_tm_difference_max_C":2.0,"full_oligo_hairpin_tm_max_C":55.0,"dimer_tm_max_C":55.0,"terminal_12mer_unique_in_defined_template":True},
                                  "initial_PCR_Ta_C":65, "optimization_gradient_C":[62,65,68], "wet_lab_optimized":False},
        "source_feature_coordinates_1based": [
            ["CMV promoter",1,589,1,"promoter","PT3027-5 p2"],
            ["SV40 polyA signal 1",1552,1557,1,"polyA_signal","PT3027-5 p2; signal motif only"],
            ["SV40 polyA signal 2",1581,1586,1,"polyA_signal","PT3027-5 p2; signal motif only"],
            ["f1 origin",1649,2104,-1,"rep_origin","PT3027-5 p2; packages EGFP noncoding strand"],
            ["KanR promoter -35",2166,2171,1,"promoter","PT3027-5 p2"],
            ["KanR promoter -10",2189,2194,1,"promoter","PT3027-5 p2"],
            ["SV40 early enhancer 1",2278,2349,1,"enhancer","PT3027-5 p2"],
            ["SV40 early enhancer 2",2350,2421,1,"enhancer","PT3027-5 p2"],
            ["SV40 origin",2445,2580,0,"rep_origin","PT3027-5 p2; no directional arrow assigned"],
            ["SV40 early promoter element",2501,2507,1,"promoter","PT3027-5 p2; not whole promoter boundary"],
            ["KanR/NeoR",2629,3423,1,"CDS","U55762.1 CDS; Kan for bacteria, G418 for mammalian cells"],
            ["HSV TK polyA signal 1",3659,3664,1,"polyA_signal","PT3027-5 p2; signal motif only"],
            ["HSV TK polyA signal 2",3672,3677,1,"polyA_signal","PT3027-5 p2; signal motif only"],
            ["pUC origin",4008,4651,0,"rep_origin","PT3027-5 p2; origin interval, no directional arrow assigned"]],
        "coordinate_convention":"Tables: 1-based inclusive. Python intervals: 0-based half-open. Digest nominal size: top-strand cut-to-cut nt, distinct from pydna ds span.",
        "limitations":["No wet-lab validation or acquired stocks", "No transcriptome-wide specificity claim", "Tag effects and expression need cellular validation", "Reference backbone was compiled by manufacturer, not a sequenced physical stock"]
    }
    (out / "construct_spec.json").write_text(json.dumps(spec, ensure_ascii=False, indent=2) + "\n")
    summary={"status":"passed", "accessions":[current.id,old.id,protein.id,vector.id],"cds_translation_matches_NP":True,
             "historic_current_CDS_identical":True,"historic_only_difference":"10 terminal A bases", "cds_sha256":hashlib.sha256(str(cds).encode()).hexdigest(),
             "inputs":inputs}
    (out / "reference_verification.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print("Reference checks passed; explicit construct specification written.")

if __name__ == "__main__":
    main()
