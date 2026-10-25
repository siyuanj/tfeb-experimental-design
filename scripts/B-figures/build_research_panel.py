#!/usr/bin/env python3
"""Derive an explicitly specified C-tail/tag factorial panel from frozen B-v1."""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from handoff import ROOT, verify_gate, assert_unchanged

OUT = ROOT/'设计/research-v3'
SITES = (466,467,469)
DESIGNS = [('3A','EGFP'),('3T','EGFP'),('WT','untagged'),('3A','untagged'),('3T','untagged')]

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    _, gate = verify_gate('docs/协作/sequence-ready.json')
    sources=[ROOT/'设计/construct.gb', ROOT/'资料/序列/TFEB_cds.fasta',Path(__file__).resolve()]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in sources}
    base=SeqIO.read(sources[0],'genbank')
    native=SeqIO.read(sources[1],'fasta').seq
    assert len(native)==1431 and str(native[-3:])=='TGA'
    assert base.seq[630:2058]==native[:-3]
    OUT.mkdir(exist_ok=True)
    records=[]
    for variant,tag in DESIGNS:
        record=deepcopy(base); dna=list(str(base.seq)); mutations=[]
        for site in SITES if variant!='WT' else ():
            offset=630+3*(site-1); old=str(base.seq[offset:offset+3])
            assert old=='AGC'
            new='GCC' if variant=='3A' else 'ACC'
            dna[offset:offset+3]=new
            mutations.append({'residue':site,'from_codon':old,'to_codon':new,'codon_start_1based':offset+1})
        dna=''.join(dna)
        if tag=='untagged':
            dna=dna[:2058]+'TGA'+dna[2799:]
            features=[]
            for f in record.features:
                name=f.qualifiers.get('label',[''])[0]
                if name=='TFEB-EGFP fusion ORF':continue
                start,end=int(f.location.start),int(f.location.end)
                if start>=2058 and end<=2799:continue
                if start>=2799:
                    f.location=FeatureLocation(start-738,end-738,strand=f.location.strand)
                elif end>2058:
                    raise ValueError('Unsupported annotation crossing replacement: '+name)
                features.append(f)
            record.features=features
            record.features.append(SeqFeature(FeatureLocation(630,2061,strand=1),type='CDS',
                qualifiers={'label':['untagged TFEB ORF']}))
        record.seq=Seq(dna)
        name=f'TFEB-{variant}-{tag}'
        record.name='TFEB_'+variant+'_'+('U' if tag=='untagged' else 'G')
        record.id=record.name
        record.description=f'Proposed {name}; derived B-v1 sequence; not an experimental clone'
        for site in SITES if variant!='WT' else ():
            offset=630+3*(site-1)
            record.features.append(SeqFeature(FeatureLocation(offset,offset+3,strand=1),type='misc_feature',
                qualifiers={'label':[f'S{site}'+variant[-1]],'note':['Designed substitution; NP_009093.1 numbering']}))
        if variant!='WT':
            record.features.append(SeqFeature(FeatureLocation(2025,2037,strand=1),type='misc_feature',
                qualifiers={'label':[variant],'note':['Cluster S466/S467/S469; intervening F468 is unchanged']}))
        end=2061 if tag=='untagged' else 2799
        protein=record.seq[630:end].translate()
        assert str(protein).endswith('*') and '*' not in str(protein[:-1])
        assert len(protein)-1==(476 if tag=='untagged' else 722)
        for f in record.features:
            if f.type=='CDS' and f.qualifiers.get('label',[''])[0] in ('untagged TFEB ORF','TFEB-EGFP fusion ORF'):
                f.qualifiers['translation']=[str(protein[:-1])]
        SeqIO.write(record,OUT/(name+'.gb'),'genbank')
        SeqIO.write(record,OUT/(name+'.fasta'),'fasta')
        SeqIO.write(SeqRecord(protein[:-1],id=name,description='Designed translation'),OUT/(name+'-protein.fasta'),'fasta')
        records.append({'name':name,'variant':variant,'tag':tag,'bp':len(record),'protein_aa':len(protein)-1,
                        'mutations':mutations,'orf_start_1based':631,'orf_end_1based':end})
    assert_unchanged(gate)
    assert hashes=={str(p.relative_to(ROOT)):sha(p) for p in sources}
    files=[p for p in OUT.iterdir() if p.suffix in ('.gb','.fasta')]
    manifest={'status':'designed; independent audit required','records':records,'inputs':hashes,
              'outputs':{str(p.relative_to(ROOT)):sha(p) for p in files},
              'design_choices':{'A_codon':'GCC','T_codon':'ACC','untagged_replacement':'B-v1 2059..2799 replaced with native TGA; delta -738 bp'},
              'limitations':['No wet-lab validation','3T may be phosphorylated; neither 3T nor 3A establishes a physiological phosphorylation state',
                             'Untagged design uses whole-plasmid synthesis and sequencing; original BamHI/Sanger plan does not transfer']}
    (OUT/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(records,indent=2))

if __name__=='__main__':main()
