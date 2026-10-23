#!/usr/bin/env python3
"""Publish P5 only when primary and independent evidence still match files."""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]

def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()

def check(path,value):
    if sha(path)!=value:
        raise ValueError(f"Stale evidence: {path}")

def main():
    primary_path=ROOT/"设计/verification_summary.json"
    primary=json.loads(primary_path.read_text())
    independent=json.loads((ROOT/"设计/independent_verification.json").read_text())
    if primary["status"] not in ("passed_primary_path","passed") or not all(primary["checks"].values()):
        raise ValueError("Primary validation failed")
    if independent["status"]!="pass" or not independent["success"]:
        raise ValueError("Independent validation failed")
    if len(independent["negative_controls"])<4 or not all(c["rejected"] and c["reason_matches"] for c in independent["negative_controls"]):
        raise ValueError("Negative controls failed")
    for group in (primary["inputs"],primary["outputs"],primary["scripts"]):
        for path,value in group.items(): check(path,value)
    for item in independent["inputs"]: check(item["path"],item["file_sha256"])
    check(independent["script_path"],independent["script_sha256"])
    primary["status"]="passed"
    primary["independent_path"]={"path":"设计/independent_verification.json","sha256":sha("设计/independent_verification.json"),
                                 "status":"pass","check_groups":len(independent["checks"]),"negative_controls_passed":len(independent["negative_controls"])}
    primary_path.write_text(json.dumps(primary,indent=2,ensure_ascii=False)+"\n")
    files={**primary["inputs"],**primary["outputs"],**primary["scripts"],
           **{i["path"]:i["file_sha256"] for i in independent["inputs"]}}
    for path in ["设计/verification_summary.json","设计/independent_verification.json","设计/reference_verification.json",
                 "设计/environment.json","scripts/B-sequence/independent_verify.py","资料/序列/TFEB_protein.fasta",
                 "资料/序列/TFEB_cds_no_stop.fasta","资料/序列/pEGFP-N1_U55762.1.fasta"]:
        files[path]=sha(path)
    ready={"status":"ready","stage":"P5","design_version":primary["design_version"],"published_at_utc":datetime.now(timezone.utc).isoformat(),
           "construct":"pEGFP-N1-TFEB-Bv1","files":files,"metrics":primary["metrics"],
           "validation":"pydna primary path + independent standard-library reconstruction/translation and four semantic negative controls",
           "limitations":primary["limitations"],
           "downstream_rule":"Any change to sequence, primers, spec or bound files requires withdrawal and revalidation of affected methods/report ready markers."}
    (ROOT/"docs/协作/sequence-ready.json").write_text(json.dumps(ready,indent=2,ensure_ascii=False)+"\n")
    print(f"P5 ready published; {len(files)} bound files; independent checks and four negative controls passed.")

if __name__=="__main__":
    main()
