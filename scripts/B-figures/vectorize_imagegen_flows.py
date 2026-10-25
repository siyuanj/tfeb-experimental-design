#!/usr/bin/env python3
"""Trace ImageGen concepts to real SVG paths and add editable exact labels.

The generated raster originals are immutable. Tracing is approximate, and the
finished illustrations require visual/scientific review. No scientific values
or experimental data are inferred from the raster artwork.
"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tools/figure-libs'))
import hashlib
import json
import re
import xml.etree.ElementTree as ET
from datetime import datetime
from zoneinfo import ZoneInfo
from importlib.metadata import version
import cairosvg
import vtracer
from publication_style import NPG, BLUE, TEAL, INK

FIG=ROOT/'报告/images/B'
SOURCE=FIG/'imagegen-v3'
NS='http://www.w3.org/2000/svg'
ET.register_namespace('',NS)
PALETTE=[BLUE,TEAL,NPG[1],NPG[6],'#FFFFFF']
TRACE_OPTIONS=dict(colormode='color',hierarchical='stacked',mode='spline',
                   filter_speckle=12,color_precision=5,layer_difference=48,
                   corner_threshold=60,length_threshold=4.0,max_iterations=10,
                   splice_threshold=45,path_precision=2)

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rgb(c):return tuple(int(c[i:i+2],16) for i in (1,3,5))
def nearest(c):
    # Quantize traced fills to the fixed manuscript palette, never data colors.
    v=rgb(c)
    return min(PALETTE,key=lambda p:sum((a-b)**2 for a,b in zip(v,rgb(p))))

def element(parent,tag,**attrs):
    return ET.SubElement(parent,f'{{{NS}}}{tag}',{k.replace('_','-'):str(v) for k,v in attrs.items()})

def text(parent,x,y,label,size=42,color=INK,bold=False):
    e=element(parent,'text',x=x,y=y,font_family='Arial, Helvetica, sans-serif',
              font_size=size,text_anchor='middle',fill=color,font_weight='700' if bold else '400')
    e.text=label

def arrow(parent,x1,x2,y):
    element(parent,'path',d=f'M {x1} {y} H {x2-14}',fill='none',stroke=INK,stroke_width=3)
    element(parent,'path',d=f'M {x2-16} {y-9} L {x2} {y} L {x2-16} {y+9} Z',fill=INK)

def render(kind):
    png=SOURCE/f'{kind}-source.png';raw=SOURCE/f'{kind}-trace.svg'
    vtracer.convert_image_to_svg_py(str(png),str(raw),**TRACE_OPTIONS)
    traced=ET.fromstring(raw.read_text())
    w=float(traced.attrib['width']);h=float(traced.attrib['height'])
    assert w/h > 2.8 and w/h < 3.2,(w,h)
    height=660 if kind=='cloning' else 940
    root=ET.Element(f'{{{NS}}}svg',{'width':str(int(w)),'height':str(height),
                    'viewBox':f'0 0 {w} {height}','role':'img'})
    title=element(root,'title');title.text='Proposed cloning workflow' if kind=='cloning' else 'Planned evidence structure'
    desc=element(root,'desc');desc.text='ImageGen concept illustrations traced to vector paths; labels are editable SVG text. No experimental data.'
    element(root,'rect',width=w,height=height,fill='#FFFFFF')
    group=element(root,'g',id='imagegen-vector-art',transform=f'translate(0,{-100 if kind=="cloning" else 0})')
    paths=0
    for node in traced:
        if node.tag.endswith('path'):
            fill=node.get('fill')
            assert fill and re.fullmatch(r'#[0-9A-Fa-f]{6}',fill),fill
            node.set('fill',nearest(fill));group.append(node);paths+=1
        else:
            raise ValueError(f'Unexpected traced SVG element: {node.tag}')
    assert 5 < paths < 3000,paths
    overlay=element(root,'g',id='editable-labels-and-connectors')
    if kind=='cloning':
        labels=[('Prepare insert','PCR + purification'),('Assemble plasmid','Digest + ligation'),
                ('Recover clones','Selection + screening'),('Verify identity','Bidirectional sequence')]
        centers=[w*(i+.5)/4 for i in range(4)]
        for i,(x,(head,sub)) in enumerate(zip(centers,labels)):
            text(overlay,x,60,f'0{i+1}',30,TEAL,True)
            text(overlay,x,535,head,42,INK,True);text(overlay,x,595,sub,33)
            if i<3:arrow(overlay,x+w*.11,centers[i+1]-w*.11,h*.48-100)
        name='cloning-workflow-imagegen-v3'
    else:
        labels=[('Identity','Sequence','Intact protein'),
                ('Cellular response','Localization + recovery','Site / pathway controls'),
                ('Functional output','TFEB loss + rescue','Transcription + turnover')]
        centers=[w*(i+.5)/3 for i in range(3)]
        for i,(x,(head,a,b)) in enumerate(zip(centers,labels)):
            text(overlay,x,730,head,44,BLUE,True)
            element(overlay,'path',d=f'M {x-w*.115} 770 H {x+w*.115}',stroke=NPG[6],stroke_width=2)
            text(overlay,x,830,a,34);text(overlay,x,892,b,34)
            if i<2:arrow(overlay,x+w*.13,centers[i+1]-w*.13,h*.5)
        name='evidence-design-imagegen-v3'
    svg=FIG/f'{name}.svg';ET.ElementTree(root).write(svg,encoding='unicode',xml_declaration=True)
    xml=svg.read_text()
    assert '<image' not in xml and 'base64' not in xml and 'foreignObject' not in xml
    assert len(list(root.iter(f'{{{NS}}}text')))==(12 if kind=='cloning' else 9)
    pdf=svg.with_suffix('.pdf');preview=svg.with_suffix('.png')
    cairosvg.svg2pdf(url=str(svg),write_to=str(pdf))
    cairosvg.svg2png(url=str(svg),write_to=str(preview),output_width=1800)
    return [raw,svg,pdf,preview],{'name':name,'traced_paths':paths,'editable_text':12 if kind=='cloning' else 9,
                                'embedded_rasters':0,'source_width':w,'source_height':h}

def main():
    inputs=[SOURCE/'cloning-source.png',SOURCE/'evidence-source.png',
            ROOT/'docs/协作/流程图ImageGen提示词-v3.json',Path(__file__).resolve(),
            Path(__file__).with_name('publication_style.py'),Path(__file__).with_name('requirements-publication.txt')]
    before={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    outputs=[];reviews=[]
    for kind in ('cloning','evidence'):
        files,review=render(kind);outputs+=files;reviews.append(review)
    assert before=={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    payload={'status':'vectorized; visual review required','created_at':datetime.now(ZoneInfo('America/New_York')).isoformat(),
             'sources':before,'outputs':{str(p.relative_to(ROOT)):sha(p) for p in outputs},
             'figures':reviews,'software':{k:version(k) for k in ('vtracer','CairoSVG')},
             'trace_options':TRACE_OPTIONS,'palette':PALETTE,
             'method':'Built-in ImageGen originals -> VTracer spline paths -> NPG fill quantization -> editable SVG labels/arrows -> vector PDF',
             'limitations':['Raster-to-vector tracing approximates contours; visual inspection is required',
                            'Conceptual illustrations only, not molecular geometry or measured data',
                            'Initial cloning source rejected for incomplete ladder motif; corrected with built-in ImageGen']}
    (FIG/'imagegen-flow-v3.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(reviews))

if __name__=='__main__':main()
