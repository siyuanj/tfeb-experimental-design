#!/usr/bin/env python3
"""Independently check the exact intended nucleotide edits and tag removal."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime
from zoneinfo import ZoneInfo
import json
from Bio import SeqIO
from check_control_panel import translate

ROOT=Path(__file__).resolve().parents[2]
def digest(p):return sha256(p.read_bytes()).hexdigest()

def main():
    out=ROOT/'设计/research-v3'; manifest=json.loads((out/'manifest.json').read_text())
    for rel,h in {**manifest['inputs'],**manifest['outputs']}.items():assert digest(ROOT/rel)==h,rel
    wt=str(SeqIO.read(ROOT/'设计/construct.gb','genbank').seq)
    native=str(SeqIO.read(ROOT/'资料/序列/TFEB_NP_009093.1.gp','genbank').seq)
    edits={'WT':{},'3A':{2026:'G',2027:'C',2029:'G',2030:'C',2035:'G',2036:'C'},
           '3T':{2027:'C',2030:'C',2036:'C'}}
    names=[('3A','EGFP'),('3T','EGFP'),('WT','untagged'),('3A','untagged'),('3T','untagged')]
    reports=[]
    for variant,tag in names:
        name=f'TFEB-{variant}-{tag}';record=SeqIO.read(out/(name+'.gb'),'genbank');dna=str(record.seq)
        expected=list(wt)
        for pos,base in edits[variant].items():expected[pos-1]=base
        expected=''.join(expected)
        if tag=='untagged':expected=expected[:2058]+'TGA'+expected[2799:]
        assert dna==expected,name
        assert len(dna)==(5396 if tag=='untagged' else 6134)
        assert record.annotations['topology']=='circular'
        native_expected=list(native)
        for pos in (466,467,469) if variant!='WT' else ():
            assert native[pos-1]=='S';native_expected[pos-1]=variant[-1]
        core=translate(dna[630:2058]);assert core==''.join(native_expected)
        end=2061 if tag=='untagged' else 2799;aa=translate(dna[630:end])
        assert aa[-1]=='*' and '*' not in aa[:-1]
        if tag=='EGFP':assert dna[2058:]==wt[2058:]
        else:
            assert dna[2058:2061]=='TGA' and dna[2061:]==wt[2799:]
            assert wt[2079:2796] not in dna
        assert dna.count('GGATCC')==(0 if tag=='untagged' else 1)
        assert dna.count('AAGCTT')==1 and dna.index('AAGCTT')==621
        assert str(SeqIO.read(out/(name+'.fasta'),'fasta').seq)==dna
        assert str(SeqIO.read(out/(name+'-protein.fasta'),'fasta').seq)==aa[:-1]
        for f in record.features:
            assert 0<=int(f.location.start)<int(f.location.end)<=len(dna)
            label=f.qualifiers.get('label',[''])[0]
            if label in ('CMV promoter','KanR/NeoR','pUC origin','f1 origin'):
                old=next(g for g in SeqIO.read(ROOT/'设计/construct.gb','genbank').features if g.qualifiers.get('label',[''])[0]==label)
                assert str(f.extract(record.seq))==str(old.extract(SeqIO.read(ROOT/'设计/construct.gb','genbank').seq))
            if label in ('untagged TFEB ORF','TFEB-EGFP fusion ORF'):assert f.qualifiers['translation']==[aa[:-1]]
        reports.append({'name':name,'bp':len(dna),'protein_aa':len(aa)-1,'specific_edits_1based':edits[variant],
                        'tag_removed':tag=='untagged','BamHI_site_count':dna.count('GGATCC')})
    checks={k:True for k in ('all_input_output_hashes','exact_full_DNA_sequence','native_residue_numbering_and_translation',
             'terminal_stop_no_internal_stop','tag_removal_and_backbone_identity','restriction_changes_explicit',
             'FASTA_GenBank_protein_agree','annotation_bounds_and_preserved_backbone_features')}
    result={'status':'passed','checked_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
            'checks':checks,'constructs':reports,'manifest_sha256':digest(out/'manifest.json'),
            'audit_script_sha256':digest(Path(__file__).resolve()),
            'translation_audit_sha256':digest(Path(__file__).with_name('check_control_panel.py')),
            'limitations':['Independent exact-edit audit, not wet-lab or functional validation','No exact siRNA/qPCR oligo validation']}
    (out/'independent-verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':'passed','checks':len(checks),'constructs':reports},indent=2))

if __name__=='__main__':main()
