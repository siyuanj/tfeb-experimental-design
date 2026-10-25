#!/usr/bin/env python3
"""Independently verify intended mutant identity, donor assembly and bindings."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime
from zoneinfo import ZoneInfo
import json
from Bio import SeqIO
from Bio.Data.CodonTable import standard_dna_table

ROOT=Path(__file__).resolve().parents[2]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def translate(dna):
    if len(dna)%3:
        raise ValueError('Frame is not divisible by three')
    table=standard_dna_table.forward_table
    return ''.join('*' if dna[i:i+3] in standard_dna_table.stop_codons else table[dna[i:i+3]]
                   for i in range(0,len(dna),3))


def cut(dna,site):
    if dna.count(site)!=1:
        raise ValueError('Expected one literal recognition site: '+site)
    return dna.index(site)+1


def main():
    out=ROOT/'设计/controls-v2'
    manifest=json.loads((out/'manifest.json').read_text())
    checks={}
    for name,expected in {**manifest['source_inputs'],**manifest['outputs']}.items():
        assert digest(ROOT/name)==expected,name
    checks['source_and_output_hashes']=True
    native=str(SeqIO.read(ROOT/'资料/序列/TFEB_NP_009093.1.gp','genbank').seq)
    wt=str(SeqIO.read(ROOT/'设计/construct.gb','genbank').seq)
    vector=str(SeqIO.read(ROOT/'资料/序列/pEGFP-N1_U55762.1.gb','genbank').seq)
    # Coordinates are independently read from the accepted construct table:
    # TFEB 631..2058; full ORF 631..2799, inclusive, one-based.
    assert len(native)==476 and native[141]==native[210]=='S'
    assert translate(wt[630:2058])==native
    changes={'WT':{},'S211A':{1261:'G'},'S142A':{1054:'G',1055:'C'},
             'S142A-S211A':{1054:'G',1055:'C',1261:'G'}}
    aa_changes={'WT':{},'S211A':{211:'A'},'S142A':{142:'A'},'S142A-S211A':{142:'A',211:'A'}}
    summary=[]
    for variant,expected in changes.items():
        record=SeqIO.read(out/f'TFEB-{variant}-EGFP.gb','genbank')
        dna=str(record.seq)
        assert len(dna)==6134 and record.annotations['topology']=='circular'
        actual={i+1:b for i,(a,b) in enumerate(zip(wt,dna)) if a!=b}
        assert actual==expected,(variant,actual)
        aa=translate(dna[630:2799])
        wt_aa=translate(wt[630:2799])
        assert len(aa)==723 and aa[-1]=='*' and '*' not in aa[:-1]
        observed={i+1:b for i,(a,b) in enumerate(zip(wt_aa,aa)) if a!=b}
        assert observed==aa_changes[variant],(variant,observed)
        assert dna[2058:2799]==wt[2058:2799]
        assert str(SeqIO.read(out/f'TFEB-{variant}-EGFP.fasta','fasta').seq)==dna
        for f in record.features:
            if f.type=='CDS' and 'TFEB-EGFP fusion ORF' in f.qualifiers.get('label',[]):
                assert f.qualifiers['translation'][0]==aa[:-1]
        donor=str(SeqIO.read(out/f'donor-{variant}.fasta','fasta').seq)
        vleft,vright=cut(vector,'AAGCTT'),cut(vector,'GGATCC')
        dleft,dright=cut(donor,'AAGCTT'),cut(donor,'GGATCC')
        assembled=vector[:vleft]+donor[dleft:dright]+vector[vright:]
        assert assembled==dna,variant
        sites={site:cut(dna,site) for site in ('AAGCTT','GGATCC')}
        assert sites=={'AAGCTT':622,'GGATCC':2061}
        summary.append({'variant':variant,'bp':len(dna),'protein_aa':len(aa)-1,
                        'changed_nt_1based':list(actual),'changed_aa_1based':list(observed),
                        'donor_reassembly_matches':True,'nominal_digest_nt':[1439,4695]})
    checks['exact_DNA_and_protein_substitutions_only']=True
    checks['donor_reassembly_from_original_vector']=True
    checks['frame_stop_linker_EGFP_and_GenBank_translation']=True
    checks['restriction_geometry_retained']=True
    empty=SeqIO.read(out/'pEGFP-N1-annotated.gb','genbank')
    assert str(empty.seq)==vector
    checks['empty_vector_sequence_unchanged']=True
    # Every plotted interval must match the corresponding exported annotation.
    map_inputs={'empty':empty,**{v:SeqIO.read(out/f'TFEB-{v}-EGFP.gb','genbank') for v in changes}}
    for name,features in manifest['displayed_features'].items():
        record=map_inputs[name]
        for drawn in features:
            match=[f for f in record.features if drawn['name'] in f.qualifiers.get('label',[])]
            assert len(match)==1
            f=match[0]
            assert int(f.location.start)+1==drawn['start_1based']
            assert int(f.location.end)==drawn['end_1based']
    checks['all_selected_map_intervals_match_GenBank']=True
    payload={'status':'passed','checked_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
             'checks':checks,'constructs':summary,'manifest_sha256':digest(out/'manifest.json'),
             'audit_script_sha256':digest(Path(__file__).resolve()),
             'limitations':['No experimental synthesis, sequencing or phenotypes',
                            'Only specified point mutants and donor joins tested; no siRNA/qPCR oligo validation']}
    (out/'independent-verification.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'status':payload['status'],'checks':len(checks),'constructs':summary},indent=2))


if __name__=='__main__':
    main()
