#!/usr/bin/env python3
"""Independently check the linker-only constructs from the frozen WT sequence.

The expected DNA is rebuilt here from the raw B-v1 sequence and explicit edits, not
from the builder's code path; translation uses the standard-library codon table.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime
from zoneinfo import ZoneInfo
import json
from Bio import SeqIO
from check_control_panel import translate

ROOT = Path(__file__).resolve().parents[2]


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def main():
    out = ROOT / '设计/linker-v4'
    manifest = json.loads((out / 'manifest.json').read_text())
    for rel, h in {**manifest['inputs'], **manifest['outputs']}.items():
        assert digest(ROOT / rel) == h, rel
    wt = str(SeqIO.read(ROOT / '设计/construct.gb', 'genbank').seq)
    native = str(SeqIO.read(ROOT / '资料/序列/TFEB_NP_009093.1.gp', 'genbank').seq)
    edits = {'WT': {}, '3A': {2026: 'G', 2027: 'C', 2029: 'G', 2030: 'C', 2035: 'G', 2036: 'C'}}
    assert wt[2058:2079] == 'GGGGATCCACCGGTCGCCACC' and translate(wt[2058:2079]) == 'GDPPVAT'
    assert wt[2079:2082] == 'ATG' and wt[2796:2799] == 'TAA'
    reports = []
    for variant in ('WT', '3A'):
        name = f'TFEB-{variant}-linker'
        record = SeqIO.read(out / (name + '.gb'), 'genbank')
        dna = str(record.seq)
        expected = list(wt)
        for pos, base in edits[variant].items():
            expected[pos - 1] = base
        expected = ''.join(expected)
        expected = expected[:2079] + 'TAA' + expected[2799:]
        assert dna == expected, name
        assert len(dna) == 5417 and record.annotations['topology'] == 'circular'
        aa = translate(dna[630:2082])
        assert aa.endswith('*') and '*' not in aa[:-1] and len(aa) - 1 == 483
        native_expected = list(native)
        for pos in (466, 467, 469) if variant == '3A' else ():
            assert native[pos - 1] == 'S'
            native_expected[pos - 1] = 'A'
        assert aa[:-1] == ''.join(native_expected) + 'GDPPVAT'
        assert dna.count('AAGCTT') == 1 and dna.index('AAGCTT') == 621
        assert dna.count('GGATCC') == 1 and dna.index('GGATCC') == 2060
        assert 'ATGGTGAGCAAGGGCGAGGAG' not in dna                       # EGFP start is gone
        assert wt[2079:2796] not in dna
        assert dna[:2079] == expected[:2079] and dna[2082:] == wt[2799:]
        hind_cut, bam_cut = 622, 2061                                    # last base before each top-strand cut
        insert, backbone = bam_cut - hind_cut, len(dna) - (bam_cut - hind_cut)
        assert (insert, backbone) == (1439, 3978)
        assert str(SeqIO.read(out / (name + '.fasta'), 'fasta').seq) == dna
        assert str(SeqIO.read(out / (name + '-protein.fasta'), 'fasta').seq) == aa[:-1]
        base_record = SeqIO.read(ROOT / '设计/construct.gb', 'genbank')
        for f in record.features:
            assert 0 <= int(f.location.start) < int(f.location.end) <= len(dna)
            label = f.qualifiers.get('label', [''])[0]
            if label in ('CMV promoter', 'KanR/NeoR', 'pUC origin', 'f1 origin'):
                old = next(g for g in base_record.features if g.qualifiers.get('label', [''])[0] == label)
                assert str(f.extract(record.seq)) == str(old.extract(base_record.seq)), label
            if label == 'TFEB-linker ORF':
                assert f.qualifiers['translation'] == [aa[:-1]]
        labels = [f.qualifiers.get('label', [''])[0] for f in record.features]
        assert 'EGFP' not in labels and 'TFEB-EGFP fusion ORF' not in labels
        reports.append({'name': name, 'bp': len(dna), 'protein_aa': len(aa) - 1,
                        'specific_edits_1based': edits[variant], 'HindIII_BamHI_insert_nt': insert,
                        'HindIII_BamHI_backbone_nt': backbone})
    checks = {k: True for k in (
        'all_input_output_hashes', 'exact_full_DNA_sequence_from_frozen_WT', 'native_residue_numbering_and_translation',
        'single_TAA_stop_no_internal_stop', 'EGFP_fully_removed', 'HindIII_BamHI_sites_retained_and_unique',
        'digest_geometry', 'FASTA_GenBank_protein_agree', 'annotation_bounds_and_preserved_backbone_features')}
    result = {'status': 'passed', 'checked_at': datetime.now(ZoneInfo('America/New_York')).isoformat(),
              'checks': checks, 'constructs': reports, 'manifest_sha256': digest(out / 'manifest.json'),
              'audit_script_sha256': digest(Path(__file__).resolve()),
              'translation_audit_sha256': digest(ROOT / 'scripts/B-audit/check_control_panel.py'),
              'limitations': ['Sequence design only; no synthesis, cloning or expression result']}
    (out / 'independent-verification.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': 'passed', 'constructs': reports}, ensure_ascii=False))


if __name__ == '__main__':
    main()
