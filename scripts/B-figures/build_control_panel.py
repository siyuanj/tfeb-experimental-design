#!/usr/bin/env python3
"""Derive point-mutant sequences and reproducible control-design figures.

No wet-lab data are generated. Existing B-v1 inputs are read-only. GCC is the
explicit alanine codon choice; this is not a codon-optimization claim.
"""
from copy import deepcopy
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
import hashlib
import json

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import HindIII, BamHI
from handoff import ROOT, verify_gate, assert_unchanged

OUT = ROOT / '设计/controls-v2'
FIG = ROOT / '报告/images/B'
from publication_style import FEATURE_COLORS as PALETTE
VARIANTS = [('WT', []), ('S211A', [211]), ('S142A', [142]),
            ('S142A-S211A', [142, 211])]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def label(feature):
    return feature.qualifiers.get('label', [''])[0]


def main():
    _, gate = verify_gate('docs/协作/sequence-ready.json')
    inputs = [ROOT/'设计/construct.gb', ROOT/'设计/PCR_product.fasta',
              ROOT/'资料/序列/pEGFP-N1_U55762.1.gb', ROOT/'设计/construct_spec.json',
              ROOT/'资料/序列/TFEB_NP_009093.1.gp', Path(__file__).resolve(), Path(__file__).with_name('publication_style.py')]
    before = {str(p.relative_to(ROOT)):sha(p) for p in inputs}
    wt = SeqIO.read(inputs[0], 'genbank')
    tfeb = next(f for f in wt.features if label(f) == 'TFEB')
    fusion = next(f for f in wt.features if f.type == 'CDS' and label(f) == 'TFEB-EGFP fusion ORF')
    start = int(tfeb.location.start)
    donor = SeqIO.read(inputs[1], 'fasta')
    donor_start = str(donor.seq).find(str(tfeb.extract(wt.seq)))
    assert donor_start >= 0 and str(donor.seq).count(str(tfeb.extract(wt.seq))) == 1
    native = str(SeqIO.read(inputs[4], 'genbank').seq)
    assert str(tfeb.extract(wt.seq).translate()) == native
    OUT.mkdir(parents=True, exist_ok=True)
    panel, mutation_rows = {}, []
    for name, sites in VARIANTS:
        record = deepcopy(wt)
        bases, donor_bases = list(str(wt.seq)), list(str(donor.seq))
        changes = []
        for aa in sites:
            offset = start + 3*(aa-1)
            old = str(wt.seq[offset:offset+3])
            assert str(Seq(old).translate()) == native[aa-1] == 'S'
            bases[offset:offset+3] = 'GCC'
            donor_offset = donor_start+3*(aa-1)
            donor_bases[donor_offset:donor_offset+3] = 'GCC'
            change = {'variant':name, 'substitution':f'S{aa}A', 'residue_1based':aa,
                      'WT_codon':old, 'mutant_codon':'GCC',
                      'plasmid_codon_start_1based':offset+1,
                      'plasmid_codon_end_1based':offset+3,
                      'changed_bases_1based':[offset+i+1 for i,(a,b) in enumerate(zip(old,'GCC')) if a != b]}
            changes.append(change)
            mutation_rows.append(change)
            record.features.append(SeqFeature(FeatureLocation(offset,offset+3,strand=1),
                type='misc_feature', qualifiers={'label':[f'S{aa}A'],
                'note':[f'Designed {old} to GCC; numbering NP_009093.1; no wet-lab validation']}))
        record.seq = Seq(''.join(bases))
        record.id = 'TFEB_'+name
        record.name = record.id
        record.description = f'Proposed pEGFP-N1 TFEB-{name}-EGFP; derived from B-v1; not experimentally validated'
        for feature in record.features:
            if feature.type == 'CDS' and label(feature) == 'TFEB-EGFP fusion ORF':
                feature.qualifiers['translation'] = [str(feature.extract(record.seq).translate(to_stop=True))]
        protein = str(fusion.extract(record.seq).translate())
        assert len(protein) == 723 and protein[-1] == '*' and '*' not in protein[:-1]
        expected = list(str(fusion.extract(wt.seq).translate()))
        for aa in sites:
            expected[aa-1] = 'A'
        assert protein == ''.join(expected)
        assert len(record) == len(wt)
        for enzyme in (HindIII,BamHI):
            assert enzyme.search(record.seq,linear=False) == enzyme.search(wt.seq,linear=False)
        stem = OUT/f'TFEB-{name}-EGFP'
        SeqIO.write(record, stem.with_suffix('.gb'), 'genbank')
        SeqIO.write(record, stem.with_suffix('.fasta'), 'fasta')
        SeqIO.write(SeqRecord(Seq(''.join(donor_bases)), id='donor_'+name,
                    description='Proposed synthetic donor with B-v1 cloning flanks'), OUT/f'donor-{name}.fasta','fasta')
        assert str(SeqIO.read(stem.with_suffix('.gb'),'genbank').seq) == str(record.seq)
        panel[name] = (record, changes)
    empty = SeqIO.read(inputs[2], 'genbank')
    spec = json.loads(inputs[3].read_text())
    # The original GenBank lacks several backbone annotations; use the existing
    # manufacturer-coordinate table rather than estimating them from a drawing.
    empty.features = [SeqFeature(FeatureLocation(a-1,b,strand=s),type=t,
                       qualifiers={'label':[name],'note':[source]})
                      for name,a,b,s,t,source in spec['source_feature_coordinates_1based']]
    empty.features.append(SeqFeature(FeatureLocation(678,1398,strand=1),type='CDS',
                                    qualifiers={'label':['EGFP']}))
    empty.description = 'U55762.1 sequence; selected existing manufacturer annotations added in a derived copy'
    SeqIO.write(empty, OUT/'pEGFP-N1-annotated.gb','genbank')
    maps = {}
    for name, record in [('empty', empty)] + [(name, pair[0]) for name,pair in panel.items()]:
        maps[name] = [{'name':label(f), 'start_1based':int(f.location.start)+1,
                       'end_1based':int(f.location.end), 'strand':f.location.strand,
                       'color':PALETTE[label(f)]} for f in record.features
                      if label(f) in PALETTE and label(f) != 'TFEB-EGFP fusion ORF']
    assert_unchanged(gate)
    assert before == {str(p.relative_to(ROOT)):sha(p) for p in inputs}
    outputs=list(OUT.glob('*.gb'))+list(OUT.glob('*.fasta'))
    manifest={'status':'computed; independent audit and visual inspection required',
              'created_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
              'design_version':'controls-v2','source_inputs':before,'sequence_gate':gate['sha256'],
              'mutations':mutation_rows,'displayed_features':maps,
              'outputs':{str(p.relative_to(ROOT)):sha(p) for p in outputs},
              'limitations':['No wet-lab clone, phenotype or synthesis QC',
                             'New siRNA targets and qPCR primers are not sequence-finalized',
                             'Point mutants cannot be distinguished by ordinary HindIII/BamHI size checks']}
    (OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'variants':list(panel),'mutations':mutation_rows,'output_files':len(outputs)},indent=2))


if __name__=='__main__':
    main()
