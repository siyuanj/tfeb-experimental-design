#!/usr/bin/env python3
"""Merge primary sequence/design and methods evidence into one source ledger."""
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def main():
    rows=[]
    def add(identifier,title,url,version,use,location,path="",classification="source fact",notes=""):
        file=ROOT/path if path else None
        rows.append(dict(source_id=identifier,title=title,url_or_doi=url,version_or_accessed=version,use=use,evidence_location=location,
                         raw_file_path=path or "not stored; use URL and evidence location",raw_file_sha256=hashlib.sha256(file.read_bytes()).hexdigest() if file and file.exists() else "not applicable",
                         classification=classification,notes=notes))
    for i,r in enumerate(json.loads((ROOT/"资料/序列/source_manifest.json").read_text()),1):
        add(f"S{i:02}",r["accession_version"],r["url"],r["accession_version"]+"; retrieved 2026-10-23","TFEB/version cross-check or complete empty-vector sequence","GenBank FEATURES and ORIGIN",r["path"])
    add("S05","Clontech pEGFP-N1 Vector Information","https://engibody.com/UPLOADEDFILES/VIAFCKEDITOR/files/AT1670-pEGFP-N1-Vector-Information.pdf","PT3027-5 / PR93634, 1999-03-19; accessed 2026-10-23","Feature coordinate mapping, marker identity, fusion logic","pp1-2 feature locations; p3 reference-sequence limitation","docs/核查/B/sources/Clontech_PT3027-5_PR93634.pdf",notes="Manufacturer original on third-party mirror; raw PDF local/ignored; not fully sequenced physical stock")
    add("S06","Roczniak-Ferguson et al. TFEB links mTORC1 signaling to lysosome homeostasis","https://doi.org/10.1126/scisignal.2002790","2012; accessed 2026-10-23","C-terminal EGFP precedent in HeLa nutrient response","Methods Cell culture and transfection; Figure3; local text lines122/157","docs/核查/B/sources/RoczniakFerguson2012_ScienceSignaling.txt",notes="Full text reused read-only from literature worktree and primary PMC methods checked; no empirical claim for B-v1")
    add("S07","Addgene pEGFP-N1-TFEB #38119","https://www.addgene.org/38119/","accessed 2026-10-23","C-terminal EGFP, deleted stop, historical reference and enzyme pair","Gene/Insert; Cloning Information",notes="Full finished plasmid sequence not downloaded; not equated to B-v1")
    add("S08","IDT gBlocks and gBlocks HiFi Gene Fragments","https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments","accessed 2026-10-23","Proposed native CDS synthetic template procurement","Product details: HiFi1000-3000bp, sequence-verified clonal intermediate",classification="manufacturer fact + project design choice",notes="No order or physical stock; sequence-specific synthesis feasibility not checked")
    add("S09","Primer3-py thermodynamic API","https://libnano.github.io/primer3-py/api/bindings.html","2.4.0; accessed 2026-10-23","Explicit Tm and secondary-structure model inputs","calc_tm, calc_hairpin, calc_homodimer, calc_heterodimer, calc_end_stability",notes="50mM monovalent is a modeling assumption, not Q5 buffer composition")
    add("S10","BamHI-HF NEB R3136","https://www.neb.com/en-us/products/r3136-bamhi-hf","accessed 2026-10-23","Selected 3prime cloning enzyme","Properties; reaction conditions; usage notes",notes="37C, rCutSmart; cleanup before ligation; full details in methods evidence")
    add("S11","NEB Cleavage Close to End of DNA Fragments","https://www.neb.com/en-us/tools-and-resources/usage-guidelines/cleavage-close-to-the-end-of-dna-fragments","accessed 2026-10-23","6nt protective clamps","BamHI-HF and HindIII-HF rows; experimental conditions; legend",classification="manufacturer fact + design choice",notes="3-5bp+++ means50-100%; not guarantee for actual amplicon")
    add("S12","Park et al. A splicing variant of TFEB negatively regulates TFEB-autophagy pathway","https://doi.org/10.1038/s41598-021-00613-y","2021; accessed 2026-10-23","N-terminal EGFP comparator","Methods Antibodies, reagents, and plasmids",notes="pEGFP-C1 precedent, not direct tag-orientation superiority test")
    add("S13","Martina and Puertollano Rag GTPases recruit TFEB and MITF to lysosomes","https://doi.org/10.1083/jcb.201209135","2013; accessed 2026-10-23","N-terminal regulatory region and tag-effect limitation","Figure6 and N-terminal deletion/Rag interaction results",notes="Does not show all N-terminal tags inactivate TFEB")
    add("S14","Course requirement summary","local:docs/10-23 课程要求与待办.md","2026-10-23","Required B deliverables","week2 course slides pp2,3,5-7; summarized without redistributing lecture","docs/10-23 课程要求与待办.md",notes="Source course PDF remains private in old course repository")
    for r in csv.DictReader((ROOT/"docs/协作/方法来源.tsv").open(),delimiter="\t"):
        add(r["source_id"],r["title"],r["url"],r["accessed_date"],r["project_use"],r["evidence_location"],classification=r["classification"],notes=r["verified_fact"])
    path=ROOT/"资料/B来源索引.tsv"
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    print(f"Wrote {len(rows)} source entries; local-only PDFs stay ignored.")

if __name__ == "__main__":
    main()
