#!/usr/bin/env python3
"""Derive linker-only attribution controls (TFEB-GDPPVAT-stop, no EGFP) from frozen designs.

WT-linker comes from the accepted B-v1 WT fusion and 3A-linker from the research-v3
3A fusion. In both, bp 2080..2799 (EGFP including its TAA) is replaced by a single TAA
stop directly after the GDPPVAT linker, so the HindIII and BamHI sites are retained.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from handoff import ROOT, verify_gate, assert_unchanged

OUT = ROOT / '设计/linker-v4'
LINKER_START, LINKER_END = 2058, 2079   # 0-based half-open: bp 2059..2079 (GDPPVAT)
EGFP_END = 2799                         # EGFP including its TAA ends at bp 2799
STOP = 'TAA'
DELTA = (EGFP_END - LINKER_END) - len(STOP)   # 717 bp removed in total
SOURCES = {
    'WT': ROOT / '设计/construct.gb',
    '3A': ROOT / '设计/research-v3/TFEB-3A-EGFP.gb',
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def derive(record, variant):
    dna = str(record.seq)
    assert str(Seq(dna[LINKER_START:LINKER_END]).translate()) == 'GDPPVAT'
    assert dna[LINKER_END:LINKER_END + 3] == 'ATG' and dna[EGFP_END - 3:EGFP_END] == 'TAA'
    derived = deepcopy(record)
    new_dna = dna[:LINKER_END] + STOP + dna[EGFP_END:]
    features = []
    for feature in record.features:
        label = feature.qualifiers.get('label', [''])[0]
        start, end = int(feature.location.start), int(feature.location.end)
        if label in ('TFEB-EGFP fusion ORF', 'EGFP'):
            continue
        if start >= LINKER_END and end <= EGFP_END:
            continue
        if start >= EGFP_END:
            feature = deepcopy(feature)
            feature.location = FeatureLocation(start - DELTA, end - DELTA, strand=feature.location.strand)
        elif end > LINKER_END:
            raise ValueError('Unsupported annotation crossing the EGFP replacement: ' + label)
        features.append(feature)
    derived.features = features
    derived.seq = Seq(new_dna)
    orf_end = LINKER_END + len(STOP)
    protein = derived.seq[630:orf_end].translate()
    assert str(protein).endswith('*') and '*' not in str(protein[:-1])
    assert len(protein) - 1 == 483
    derived.features.append(SeqFeature(
        FeatureLocation(630, orf_end, strand=1), type='CDS',
        qualifiers={'label': ['TFEB-linker ORF'], 'translation': [str(protein[:-1])], 'codon_start': [1],
                    'note': ['TFEB 1-476 plus GDPPVAT linker; native TGA absent; TAA stop replaces EGFP']}))
    derived.features.append(SeqFeature(
        FeatureLocation(LINKER_END, orf_end, strand=1), type='misc_feature',
        qualifiers={'label': ['TAA stop'], 'note': ['Replaces EGFP (B-v1 bp 2080..2799)']}))
    derived.features.sort(key=lambda f: int(f.location.start))
    name = f'TFEB-{variant}-linker'
    derived.name = f'TFEB_{variant}_L'
    derived.id = derived.name
    derived.description = f'Proposed {name}-only; derived from frozen design; not an experimental clone'
    return name, derived, protein


def main():
    _, gate = verify_gate('docs/协作/sequence-ready.json')
    script = Path(__file__).resolve()
    hashes = {str(p.relative_to(ROOT)): sha(p) for p in [*SOURCES.values(), script]}
    OUT.mkdir(exist_ok=True)
    records = []
    for variant, source in SOURCES.items():
        name, record, protein = derive(SeqIO.read(source, 'genbank'), variant)
        SeqIO.write(record, OUT / (name + '.gb'), 'genbank')
        SeqIO.write(record, OUT / (name + '.fasta'), 'fasta')
        SeqIO.write(SeqRecord(protein[:-1], id=name, description='Designed translation'),
                    OUT / (name + '-protein.fasta'), 'fasta')
        records.append({'name': name, 'variant': variant, 'bp': len(record), 'protein_aa': len(protein) - 1,
                        'orf_start_1based': 631, 'orf_end_1based': LINKER_END + len(STOP),
                        'removed_bp': EGFP_END - LINKER_END, 'added_bp': len(STOP)})
    assert_unchanged(gate)
    assert hashes == {str(p.relative_to(ROOT)): sha(p) for p in [*SOURCES.values(), script]}
    files = [p for p in OUT.iterdir() if p.suffix in ('.gb', '.fasta')]
    manifest = {
        'status': 'designed; independent audit required', 'records': records, 'inputs': hashes,
        'outputs': {str(p.relative_to(ROOT)): sha(p) for p in sorted(files)},
        'design_choices': {'stop_codon': 'TAA, matching the EGFP stop it replaces',
                           'replacement': 'B-v1 bp 2080..2799 (EGFP including TAA) replaced by TAA; delta -717 bp',
                           'retained': 'HindIII (622) and BamHI (2061) sites and the GDPPVAT linker'},
        'limitations': ['No wet-lab validation',
                        'Whole-plasmid synthesis and sequencing; the B-v1 Sanger primers S-F5 and S-R2 bind inside EGFP and are absent']}
    (OUT / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    main()
