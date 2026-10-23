#!/usr/bin/env python3
"""Retrieve explicit versioned primary records without modifying other checkouts."""
import hashlib
import json
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "资料/序列"
RECORDS = [
    ("NM_007162.3", "nuccore", "TFEB_reference.gb"),
    ("NM_007162.2", "nuccore", "TFEB_NM_007162.2.gb"),
    ("NP_009093.1", "protein", "TFEB_NP_009093.1.gp"),
    ("U55762.1", "nuccore", "pEGFP-N1_U55762.1.gb"),
]

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for accession, database, filename in RECORDS:
        url = ("https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?"
               f"id={accession}&db={database}&report=genbank&retmode=text")
        data = urllib.request.urlopen(url, timeout=60).read()
        if not data.startswith(b"LOCUS") or b"ORIGIN" not in data:
            raise ValueError(f"Invalid primary record {accession}")
        path = OUT / filename
        if path.exists() and path.read_bytes() != data:
            raise ValueError(f"Existing source differs: {path}; review before replacing")
        path.write_bytes(data)
        record = SeqIO.read(path, "genbank")
        if record.id != accession:
            raise ValueError(f"Wrong version returned for {accession}")
        previous = ROOT / ".worktrees/literature-review/资料/序列" / filename
        manifest.append({"accession_version": accession, "path": str(path.relative_to(ROOT)),
                         "url": url, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                         "length": len(record), "file_sha256": hashlib.sha256(data).hexdigest(),
                         "sequence_sha256": hashlib.sha256(str(record.seq).upper().encode()).hexdigest(),
                         "matches_prior_worktree_file": previous.exists() and previous.read_bytes() == data})
        time.sleep(0.4)
    (OUT / "source_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(manifest, indent=2, ensure_ascii=False))

if __name__ == "__main__":
    main()
