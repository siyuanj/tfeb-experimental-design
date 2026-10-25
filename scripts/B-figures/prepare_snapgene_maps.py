#!/usr/bin/env python3
"""Make presentation-only GenBank copies using SnapGene's documented dialect.

Sequence and shown coordinates are preserved. Minor and overlapping annotations
are omitted only from these display copies; complete design records are retained.
"""
from copy import deepcopy
from pathlib import Path
import hashlib
import json
from Bio import SeqIO
from Bio.SeqFeature import Reference
from publication_style import FEATURE_COLORS, NPG, SOURCE

ROOT = Path(__file__).resolve().parents[2]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    out = ROOT/'报告/质粒图-SnapGene'
    out.mkdir(exist_ok=True)
    entries = []
    old = [(key, ROOT/'设计/controls-v2'/('pEGFP-N1-annotated.gb' if key == 'empty' else f'TFEB-{key}-EGFP.gb'))
           for key in ('empty','WT','S211A','S142A','S142A-S211A')]
    new = [(p.stem,p) for p in sorted((ROOT/'设计/research-v3').glob('*.gb'))]
    for key, source in old+new:
        record = SeqIO.read(source, 'genbank')
        display = deepcopy(record)
        title = ('pEGFP-N1' if key == 'empty' else f'TFEB-{key}-EGFP') if (key,source) in old else key
        display.name = 'Exported'
        display.id = 'Exported'
        display.annotations['keywords'] = [title]
        display.annotations['data_file_division'] = 'SYN'
        ref = Reference()
        ref.authors = '.'
        ref.title = 'Direct Submission'
        ref.journal = 'SnapGene GenBank format; prepared by project script for display only'
        display.annotations['references'] = [ref]
        display.features = []
        selected = []
        for feature in record.features:
            label = feature.qualifiers.get('label',[''])[0]
            if label not in ('CMV promoter','TFEB','EGFP','KanR/NeoR','pUC origin','f1 origin','S142A','S211A','3A','3T'):
                continue
            f = deepcopy(feature)
            color = FEATURE_COLORS.get(label, NPG[0])
            direction = '; direction: RIGHT' if label in ('CMV promoter','TFEB','EGFP','KanR/NeoR') else ''
            if label == 'f1 origin':
                direction = '; direction: LEFT'
            # A nondirectional origin remains a block, never an invented arrow.
            f.qualifiers['note'] = f.qualifiers.get('note',[]) + [f'color: {color}{direction}']
            display.features.append(f)
            selected.append({'name':label,'start_1based':int(f.location.start)+1,
                             'end_1based':int(f.location.end),'strand':f.location.strand,'color':color})
        target = out/f'{title}.gb'
        SeqIO.write(display,target,'genbank')
        reread = SeqIO.read(target,'genbank')
        assert str(reread.seq) == str(record.seq)
        assert len(reread.features) == len(selected)
        entries.append({'key':key,'title':title,'source':str(source.relative_to(ROOT)),
                        'source_sha256':sha(source),'display_file':str(target.relative_to(ROOT)),
                        'display_sha256':sha(target),'bp':len(record),'features':selected})
    payload = {'palette':'ggsci NPG nrc','palette_source':SOURCE,
               'format_spec':'https://support.snapgene.com/hc/en-us/articles/10242682237588-What-is-Genbank-SnapGene-Format',
               'scripts':{str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),Path(__file__).with_name('publication_style.py'))},
               'records':entries,'rendering':'Native SnapGene Viewer GUI; export images separately, preserve vendor credit'}
    (out/'display-inputs.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'records':len(entries),'output':str(out)},ensure_ascii=False))


if __name__ == '__main__':
    main()
