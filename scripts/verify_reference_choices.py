#!/usr/bin/env python3
"""Verify fixed TFEB/vector references; do not assemble a final construct."""

import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[1]
SEQUENCES = ROOT / "资料" / "序列"


def require(condition, message):
    """Fail explicitly, including when Python optimization is enabled."""
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def read_record(name, accession):
    path = SEQUENCES / name
    record = SeqIO.read(path, "genbank")
    require(record.id == accession, f"Unexpected accession in {name}")
    return record, {
        "path": str(path.relative_to(ROOT)),
        "accession_version": record.id,
        "file_sha256": sha256(path.read_bytes()),
        "sequence_sha256": sha256(str(record.seq).upper().encode()),
        "length": len(record.seq),
    }


def tfeb_cds(record):
    features = [f for f in record.features if f.type == "CDS"
                and f.qualifiers.get("gene") == ["TFEB"]]
    require(len(features) == 1, "Expected one TFEB CDS")
    feature = features[0]
    require(feature.location.strand == 1, "Unexpected TFEB strand")
    require(feature.qualifiers.get("codon_start", ["1"]) == ["1"],
            "Unexpected codon_start")
    require(feature.qualifiers["protein_id"] == ["NP_009093.1"],
            "Unexpected TFEB protein accession")
    coding = feature.extract(record.seq).upper()
    protein = coding.translate(table=1, cds=True)
    require(str(protein) == feature.qualifiers["translation"][0],
            "Translation does not match GenBank annotation")
    return feature, coding, protein


def write_fasta(name, sequence, identifier, description):
    SeqIO.write(SeqRecord(sequence, id=identifier, description=description),
                SEQUENCES / name, "fasta")


def main():
    require(not (ROOT / "docs/协作/sequence-ready.json").exists(),
            "Legacy reference writer is disabled after P5 handoff; use scripts/B-sequence/README.md")
    current, current_meta = read_record("TFEB_reference.gb", "NM_007162.3")
    historic, historic_meta = read_record("TFEB_NM_007162.2.gb", "NM_007162.2")
    protein_record, protein_meta = read_record("TFEB_NP_009093.1.gp", "NP_009093.1")
    vector, vector_meta = read_record("pEGFP-N1_U55762.1.gb", "U55762.1")
    feature, cds, protein = tfeb_cds(current)
    old_feature, old_cds, old_protein = tfeb_cds(historic)
    require(str(cds) == str(old_cds), "Historical and selected CDS differ")
    require(str(protein) == str(old_protein) == str(protein_record.seq),
            "Protein records disagree")
    require(str(historic.seq) == str(current.seq) + "A" * 10,
            "Transcript-version difference is not the expected terminal 10 A")
    require(vector.annotations.get("topology") == "circular", "Vector is not circular")
    require(len(vector) == 4733, "Unexpected pEGFP-N1 reference length")
    egfp_features = [f for f in vector.features if f.type == "CDS"
                     and f.qualifiers.get("gene") == ["egfp"]]
    require(len(egfp_features) == 1, "Expected one EGFP CDS")
    egfp = egfp_features[0]
    egfp_cds = egfp.extract(vector.seq).upper()
    egfp_protein = egfp_cds.translate(table=11, cds=True)
    require(str(egfp_protein) == egfp.qualifiers["translation"][0],
            "Vector EGFP translation disagrees with annotation")

    write_fasta("TFEB_cds.fasta", cds, "NM_007162.3_CDS",
                "human WT TFEB isoform 1; 303..1733; includes terminal TGA")
    write_fasta("TFEB_cds_no_stop.fasta", cds[:-3], "NM_007162.3_CDS_no_stop",
                "TFEB coding portion for C-terminal EGFP fusion; no cloning tails")
    write_fasta("TFEB_protein.fasta", protein, "NP_009093.1",
                "human WT TFEB isoform 1; no tag or linker")
    write_fasta("pEGFP-N1_U55762.1.fasta", vector.seq.upper(), "U55762.1",
                "pEGFP-N1 empty vector reference; circular")
    summary = {
        "scope": "Reference sequence selection only; no final fusion assembly or primers",
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "software": {"python": platform.python_version(), "biopython": Bio.__version__},
        "inputs": [current_meta, historic_meta, protein_meta, vector_meta],
        "tfeb": {
            "transcript": current.id, "protein": protein_record.id,
            "cds_start_1based": int(feature.location.start) + 1,
            "cds_end_1based_inclusive": int(feature.location.end), "strand": "+",
            "cds_length_including_stop_nt": len(cds),
            "coding_length_excluding_stop_nt": len(cds) - 3,
            "start_codon": str(cds[:3]), "native_stop": str(cds[-3:]),
            "protein_length_aa": len(protein),
            "cds_sha256": sha256(str(cds).encode()),
            "protein_sha256": sha256(str(protein).encode()),
            "historical_cds_identical": True,
            "historical_protein_identical": True,
            "historical_mrna_difference": "NM_007162.2 has 10 additional terminal A bases",
        },
        "vector": {
            "accession": vector.id, "length_bp": len(vector), "topology": "circular",
            "egfp_cds_start_1based": int(egfp.location.start) + 1,
            "egfp_cds_end_1based_inclusive": int(egfp.location.end),
            "egfp_length_aa": len(egfp_protein), "egfp_stop": str(egfp_cds[-3:]),
            "egfp_translation_matches_annotation": True,
            "annotation_limit": "Original GenBank annotates CDSs, not a complete vector feature map",
        },
        "not_verified": ["Physical template or vector stock", "Addgene #38119 full sequence",
                         "Final junctions and linker", "Primers and restriction cloning",
                         "Final fusion size or molecular weight"],
    }
    output = ROOT / "设计" / "reference_verification.json"
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
