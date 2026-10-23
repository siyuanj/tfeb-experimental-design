#!/usr/bin/env python3
"""Export display values from the accepted design, without recomputing biology."""

import json
from handoff import ROOT, assert_unchanged, sha256, verify_gate


def main():
    _, snapshot = verify_gate("docs/协作/sequence-ready.json")
    metrics = json.loads((ROOT / "设计/construct_metrics.json").read_text())
    spec = json.loads((ROOT / "设计/construct_spec.json").read_text())
    values = {
        "ConstructName": spec["construct_name"],
        "ConstructBp": f"{metrics['plasmid_bp']:,}",
        "EmptyVectorBp": f"{metrics['empty_vector_bp']:,}",
        "TFEBReferenceNt": f"{metrics['TFEB_reference_CDS_with_stop_nt']:,}",
        "TFEBCodingNt": f"{metrics['TFEB_coding_no_stop_nt']:,}",
        "PCRBp": f"{metrics['PCR_product_bp']:,}",
        "BackboneBp": f"{metrics['digested_backbone_nominal_nt']:,}",
        "InsertBp": f"{metrics['digested_insert_nominal_nt']:,}",
        "RemovedMCSBp": str(metrics['removed_vector_nominal_nt']),
        "FusionAa": str(metrics['fusion_aa']),
        "TFEBAa": str(metrics['TFEB_aa']),
        "EGFPAa": str(metrics['EGFP_aa']),
        "FusionMassKDa": f"{metrics['fusion_molecular_weight_Da'] / 1000:.2f}",
        "FusionORFNt": f"{metrics['fusion_orf_nt_including_stop']:,}",
        "LinkerPeptide": metrics['linker_peptide'],
    }
    contents = "% Generated from P5-accepted inputs; rerun write_construct_values.py.\n"
    contents += "\n".join("\\newcommand{\\" + key + "}{" + value + "}" for key, value in values.items()) + "\n"
    assert_unchanged(snapshot)
    output = ROOT / "报告/sections/construct-values.tex"
    output.write_text(contents)
    summary = {"sequence_gate": snapshot, "values": values,
               "output": {str(output.relative_to(ROOT)): sha256(output)}}
    (ROOT / "docs/协作/报告数值输入.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(values, ensure_ascii=False))


if __name__ == "__main__":
    main()
