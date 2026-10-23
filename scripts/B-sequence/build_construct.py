#!/usr/bin/env python3
"""Simulate tailed PCR, physical restriction ends, ligation and translation."""
import csv
import hashlib
import importlib.metadata as metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import primer3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction, molecular_weight
from Bio.Restriction import HindIII, BamHI, XhoI, KpnI
from pydna.amplify import pcr
from pydna.dseqrecord import Dseqrecord

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "设计"

def require(condition, message):
    if not condition:
        raise ValueError(message)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_tsv(name, rows):
    with (OUT / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

def fasta(name, sequence, identifier, description):
    SeqIO.write(SeqRecord(Seq(str(sequence)), id=identifier, description=description), OUT / name, "fasta")

def tm_args(spec):
    p = spec["primer_thermodynamics"]
    return dict(mv_conc=p["mv_conc_mM"], dv_conc=p["dv_conc_mM"],
                dntp_conc=p["dntp_conc_mM_total"], dna_conc=p["dna_conc_nM"])

def all_indices(text, pattern):
    return [i for i in range(len(text)-len(pattern)+1) if text.startswith(pattern, i)]

def main():
    spec = json.loads((OUT / "construct_spec.json").read_text())
    for item in spec["reference_inputs"]:
        require(digest(ROOT / item["path"]) == item["file_sha256"], "Source file changed")
    cds = SeqIO.read(ROOT / spec["template"]["sequence_path"], "fasta").seq.upper()
    transcript = SeqIO.read(ROOT / "资料/序列/TFEB_reference.gb", "genbank")
    reference_cds = next(f for f in transcript.features if f.type == "CDS" and f.qualifiers.get("gene") == ["TFEB"]).extract(transcript.seq).upper()
    require(cds == reference_cds, "Template FASTA differs from primary reference")
    tfeb = cds[:-3]
    vector_record = SeqIO.read(ROOT / spec["vector"]["path"], "genbank")
    vector_seq = vector_record.seq.upper()
    native_protein = SeqIO.read(ROOT / "资料/序列/TFEB_protein.fasta", "fasta").seq
    require(native_protein == cds.translate(cds=True), "Derived TFEB protein disagrees with reference")
    egfp_cds = vector_seq[678:1398]
    egfp_protein = egfp_cds.translate(cds=True)
    require(cds[-3:] == "TGA" and len(tfeb) == 1428, "TFEB stop/length differs")
    primers = []
    thermodynamics = []
    args = tm_args(spec)
    structure_args = dict(args, temp_c=spec["primer_thermodynamics"]["structure_temperature_C"])
    for row in spec["primers"]:
        anneal = str(tfeb[:row["annealing_nt"]] if row["direction"] == "forward" else tfeb[-row["annealing_nt"]:].reverse_complement())
        tail = row["protective_bases"] + row["restriction_site"] + row["adjustment"]
        full = tail + anneal
        require(len(row["protective_bases"]) == 6, "Unexpected clamp length")
        anneal_tm = primer3.calc_tm(anneal, tm_method="santalucia", salt_corrections_method="santalucia", **args)
        structures = {"hairpin": primer3.calc_hairpin(full, **structure_args).todict(),
                      "homodimer": primer3.calc_homodimer(full, **structure_args).todict()}
        template = str(cds)
        for segment in (anneal, anneal[-12:]):
            count = len(all_indices(template, segment)) + len(all_indices(str(cds.reverse_complement()), segment))
            require(count == 1, f"Nonunique primer sequence in defined template: {row['name']}")
        require(structures["hairpin"]["tm"] < 55 and structures["homodimer"]["tm"] < 55, "Unacceptable primer structure")
        primers.append(dict(name=row["name"],direction=row["direction"],sequence_5to3=full,
                            protective_bases=row["protective_bases"],restriction_site=row["restriction_site"],
                            adjustment=row["adjustment"],adjustment_role=row["adjustment_role"],annealing_sequence_5to3=anneal,
                            full_length_nt=len(full),annealing_length_nt=len(anneal),full_gc_percent=round(gc_fraction(full)*100,3),
                            annealing_gc_percent=round(gc_fraction(anneal)*100,3),annealing_tm_C=round(anneal_tm,3),
                            full_oligo_tm_C_not_for_initial_annealing=round(primer3.calc_tm(full,**args),3),
                            template_CDS_start_1based=1 if row["direction"] == "forward" else len(tfeb)-len(anneal)+1,
                            template_CDS_end_1based=len(anneal) if row["direction"] == "forward" else len(tfeb),
                            strand=1 if row["direction"] == "forward" else -1))
        thermodynamics.append({"primer":row["name"],"full_oligo":structures,"unique_full_annealing_site_in_template":True,"unique_terminal_12mer_in_template":True})
    fwd, rev = [r["sequence_5to3"] for r in primers]
    require(abs(primers[0]["annealing_tm_C"]-primers[1]["annealing_tm_C"]) <= 2, "Tm mismatch")
    hetero = primer3.calc_heterodimer(fwd, rev, **structure_args).todict()
    ends = [primer3.calc_end_stability(a,b,**structure_args).todict() for a,b in [(fwd,rev),(rev,fwd)]]
    require(hetero["tm"] < 55 and all(e["tm"] < 55 for e in ends), "Primer pair dimer risk")
    # PCR must be performed against the actual defined template, including its stop.
    amplicon = pcr(fwd, rev, Dseqrecord(cds), limit=20)
    expected_pcr = fwd[:-22] + str(tfeb) + str(Seq(rev[:-22]).reverse_complement())
    require(str(amplicon.seq).upper() == expected_pcr, "pydna PCR disagrees with literal primer-tail reconstruction")
    vector = Dseqrecord(vector_seq, circular=True)
    for enzyme in (HindIII,BamHI):
        require(not enzyme.search(tfeb), f"Internal {enzyme} site")
        require(len(enzyme.search(vector_seq, linear=False)) == 1, f"Nonunique vector {enzyme}")
    vector_parts = vector.cut(HindIII,BamHI)
    insert_parts = amplicon.cut(HindIII,BamHI)
    require(len(vector_parts) == 2 and len(insert_parts) == 3, "Unexpected digest products")
    backbone = max(vector_parts, key=len)
    insert = max(insert_parts, key=len)
    # Nonmatching cohesive ends must reject vector self-circularization.
    try:
        backbone.looped()
    except (ValueError, TypeError):
        self_ligation_rejected = True
    else:
        self_ligation_rejected = False
    require(self_ligation_rejected, "Vector ends unexpectedly self-compatible")
    assembled = (backbone + insert).looped()
    anchor = str(vector_seq[:60])
    shift = str(assembled.seq).upper().find(anchor)
    require(shift >= 0, "Cannot normalize construct to original vector origin")
    assembled = assembled.shifted(shift)
    seq = Seq(str(assembled.seq).upper())
    cuts_vector = [e.search(vector_seq,linear=False)[0]-1 for e in (HindIII,BamHI)]
    cuts_insert = [e.search(amplicon.seq)[0]-1 for e in (HindIII,BamHI)]
    nominal_backbone = len(vector_seq) - (cuts_vector[1]-cuts_vector[0])
    nominal_insert = cuts_insert[1]-cuts_insert[0]
    require(len(seq) == nominal_backbone+nominal_insert, "Ligation length mismatch")
    tfeb_start = str(seq).find(str(tfeb))
    require(tfeb_start >= 0 and str(seq).count(str(tfeb)) == 1, "TFEB copy/identity mismatch")
    egfp_start = str(seq).find(str(egfp_cds))
    tfeb_end = tfeb_start + len(tfeb)
    linker = seq[tfeb_end:egfp_start]
    require(str(linker) == spec["fusion"]["expected_linker_nt"], "Unexpected fusion junction")
    require(str(linker.translate()) == spec["fusion"]["expected_linker_aa"], "Unexpected linker translation")
    orf_end = egfp_start + len(egfp_cds)
    orf = seq[tfeb_start:orf_end]
    fusion = orf.translate(cds=True)
    require(fusion == native_protein+linker.translate()+egfp_protein, "Unexpected fusion residues")
    require(str(seq[tfeb_start-3:tfeb_start]) == "ACC" and str(seq[tfeb_start+3]) == "G", "Initiation context mismatch")
    delta = len(seq)-len(vector_seq)
    features = []
    feature_rows = []
    def add_feature(label,start,end,strand,kind,source,note=""):
        require(0 <= start < end <= len(seq), "Feature outside construct")
        q={"label":[label],"note":[source + ("; "+note if note else "")]}
        if kind == "CDS":
            q["translation"]=[str(seq[start:end].translate(cds=True))]
            q["codon_start"]=["1"]
        features.append(SeqFeature(FeatureLocation(start,end,strand=strand),type=kind,qualifiers=q))
        feature_rows.append(dict(name=label,start_1based=start+1,end_1based=end,strand=strand,type=kind,source=source,note=note))
    for label,start,end,strand,kind,source in spec["source_feature_coordinates_1based"]:
        if end <= cuts_vector[0]:
            dest_start,dest_end=start-1,end
        elif start-1 >= cuts_vector[1]:
            dest_start,dest_end=start-1+delta,end+delta
        else:
            raise ValueError(f"Source feature intersects removed segment: {label}")
        require(seq[dest_start:dest_end] == vector_seq[start-1:end], "Backbone feature sequence changed")
        add_feature(label,dest_start,dest_end,strand,kind,source)
    add_feature("TFEB-EGFP fusion ORF",tfeb_start,orf_end,1,"CDS","B-v1; NM_007162.3 + defined linker + U55762.1 EGFP","includes terminal TAA")
    add_feature("TFEB",tfeb_start,tfeb_end,1,"misc_feature","NM_007162.3 CDS 303..1730","native stop omitted; native start retained")
    add_feature("linker GDPPVAT",tfeb_end,egfp_start,1,"misc_feature","B-v1 GG frame bases + residual BamHI-to-EGFP vector segment")
    add_feature("EGFP",egfp_start,orf_end,1,"misc_feature","U55762.1 CDS 679..1398","internal ATG retained; terminal TAA")
    add_feature("initiation context ACC",tfeb_start-3,tfeb_start,1,"misc_feature","B-v1 forward-primer tail","-3 A and +4 G")
    add_feature("retained upstream MCS",590,cuts_vector[0],1,"misc_feature","U55762.1 MCS fragment","not a contiguous original MCS after cloning")
    add_feature("retained downstream MCS",egfp_start-19,egfp_start-7,1,"misc_feature","U55762.1 BamHI-to-AgeI segment","overlaps translated linker")
    for enzyme in (HindIII,BamHI):
        site=str(seq).find(enzyme.site)
        add_feature(str(enzyme)+" site",site,site+len(enzyme.site),1,"misc_feature","sequence-exact restriction recognition site")
    record=SeqRecord(seq,id="TFEB_Bv1",name="TFEB_Bv1",description=spec["construct_name"]+"; in-silico design, not experimentally validated")
    record.annotations={"molecule_type":"DNA","topology":"circular","date":"23-OCT-2026"}
    record.features=features
    SeqIO.write(record,OUT/"construct.gb","genbank")
    fasta("construct.fasta",seq,spec["construct_name"],f"circular; {len(seq)} bp; predicted design")
    fasta("fusion_protein.fasta",fusion,"TFEB-EGFP_Bv1",f"{len(fusion)} aa; native TFEB + GDPPVAT + EGFP")
    fasta("PCR_product.fasta",amplicon.seq,"TFEB_Bv1_PCR",f"{len(amplicon)} bp; includes both complete primer tails")
    fasta("fusion_orf.fasta",orf,"TFEB-EGFP_Bv1_ORF","includes terminal TAA")
    write_tsv("primers.tsv",primers)
    write_tsv("construct_features.tsv",sorted(feature_rows,key=lambda r:(r['start_1based'],r['end_1based'])))
    junction_rows=[
        dict(junction="5prime vector-insert",start_1based=tfeb_start-9+1,end_1based=tfeb_start+18,
             dna_5to3=str(seq[tfeb_start-9:tfeb_start+18]),translation="upstream AAGCTTACC untranslated | "+str(seq[tfeb_start:tfeb_start+18].translate()),interpretation="HindIII; ACC context; native TFEB ATG"),
        dict(junction="TFEB-linker-EGFP",start_1based=tfeb_end-18+1,end_1based=egfp_start+18,
             dna_5to3=str(seq[tfeb_end-18:egfp_start+18]),translation=str(seq[tfeb_end-18:egfp_start+18].translate()),interpretation="last 6 TFEB aa | GDPPVAT | first 6 EGFP aa"),
        dict(junction="terminal stop",start_1based=orf_end-18+1,end_1based=orf_end,
             dna_5to3=str(seq[orf_end-18:orf_end]),translation=str(seq[orf_end-18:orf_end].translate()),interpretation="EGFP-derived final TAA; no TFEB stop")]
    write_tsv("junctions.tsv",junction_rows)
    def end_info(fragment):
        return {"ds_span_bp":len(fragment),"watson_nt":len(fragment.seq.watson),"crick_nt":len(fragment.seq.crick),
                "five_prime_end":fragment.seq.five_prime_end(),"three_prime_end":fragment.seq.three_prime_end()}
    metrics={"plasmid_bp":len(seq),"empty_vector_bp":len(vector_seq),"TFEB_reference_CDS_with_stop_nt":len(cds),
             "TFEB_coding_no_stop_nt":len(tfeb),"PCR_product_bp":len(amplicon),"removed_vector_nominal_nt":cuts_vector[1]-cuts_vector[0],
             "digested_backbone_nominal_nt":nominal_backbone,"digested_insert_nominal_nt":nominal_insert,
             "backbone":end_info(backbone),"insert":end_info(insert),"all_vector_fragments":[end_info(f) for f in vector_parts],
             "all_PCR_fragments":[end_info(f) for f in insert_parts],"TFEB_aa":len(native_protein),"linker_aa":len(linker)//3,
             "EGFP_aa":len(egfp_protein),"fusion_aa":len(fusion),"fusion_molecular_weight_Da":molecular_weight(fusion,seq_type="protein"),
             "fusion_orf_nt_including_stop":len(orf),"fusion_start_1based":tfeb_start+1,"fusion_end_1based":orf_end,
             "TFEB_end_1based":tfeb_end,"EGFP_start_1based":egfp_start+1,"linker_sequence":str(linker),"linker_peptide":str(linker.translate()),
             "vector_top_cut_after_1based":cuts_vector,"PCR_top_cut_after_1based":cuts_insert,
             "construct_top_cut_after_1based":{str(e):[p-1 for p in e.search(seq,linear=False)] for e in (HindIII,BamHI,XhoI,KpnI)}}
    (OUT/"construct_metrics.json").write_text(json.dumps(metrics,indent=2)+"\n")
    thermo={"parameters":spec["primer_thermodynamics"],"primers":thermodynamics,"heterodimer":hetero,"directional_3prime_end_stability":ends,
            "rejected_candidate":{"forward_context":"GCCACC","full_oligo_hairpin_Tm_C":primer3.calc_hairpin("ACTGACAAGCTTGCCACC"+str(tfeb[:22]),**structure_args).tm,
                                  "reason":"66 C short hairpin above prespecified 55 C screening cutoff; select ACC context with -3 A and +4 G"},
            "specificity_scope":"Exact matches and terminal 12mer against the specified synthetic template, both strands; not genome/transcriptome validation"}
    (OUT/"primer_checks.json").write_text(json.dumps(thermo,indent=2)+"\n")
    outputs=["construct.gb","construct.fasta","fusion_protein.fasta","fusion_orf.fasta","PCR_product.fasta","primers.tsv","construct_features.tsv","junctions.tsv","construct_metrics.json","primer_checks.json"]
    summary={"status":"passed_primary_path", "design_version":spec["design_version"],"checked_at_utc":datetime.now(timezone.utc).isoformat(),
             "checks":{"source_hashes":True,"reference_translation":True,"primer_direction_and_template_match":True,"primer_thermodynamics_screen":True,
                       "no_internal_cloning_sites":True,"pydna_PCR_matches_literal_tails":True,"restriction_fragment_count":True,"physical_cohesive_end_ligation":True,
                       "vector_self_ligation_rejected":self_ligation_rejected,"full_ORF_and_terminal_stop":True,"TFEB_and_EGFP_unchanged":True,"feature_sequences_preserved":True},
             "metrics":metrics,"software":{"python":platform.python_version(),**{p:metadata.version(p) for p in ["pydna","biopython","primer3-py"]}},
             "inputs":{"设计/construct_spec.json":digest(OUT/"construct_spec.json"),**{r["path"]:digest(ROOT/r["path"]) for r in spec["reference_inputs"]}},
             "scripts":{str(p.relative_to(ROOT)):digest(p) for p in [Path(__file__).resolve(),Path(__file__).with_name("prepare_spec.py")]},
             "outputs":{"设计/"+p:digest(OUT/p) for p in outputs},"independent_path":"pending; do not publish sequence-ready until independent verification passes",
             "command":"'/Volumes/Work/课程/10-18 分子生物学实验/.venv-plasmid/bin/python' scripts/B-sequence/build_construct.py",
             "limitations":spec["limitations"]}
    (OUT/"verification_summary.json").write_text(json.dumps(summary,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({"status":summary["status"],"metrics":metrics},indent=2))

if __name__ == "__main__":
    main()
