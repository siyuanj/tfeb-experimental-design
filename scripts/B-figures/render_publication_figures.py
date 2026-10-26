#!/usr/bin/env python3
"""Compose native SnapGene exports and draw data-free experimental schematics."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'tools/figure-libs'))
import json
import hashlib
import re
import fitz
import cairosvg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, FancyBboxPatch
from publication_style import INK, BLUE, TEAL, ORANGE, NPG, SOURCE
from figure_layouts import save_figure, style
from scientific_schematics import evidence_flow
from academic_panels import schedule_panel, research_design_panel
import subprocess

FIG = ROOT/'报告/images/B'
BOLD_FONT = matplotlib.font_manager.findfont(matplotlib.font_manager.FontProperties(family='DejaVu Sans', weight='bold'))
NATIVE = ROOT/'报告/质粒图-SnapGene'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def compose_native(stems, name, cols=2):
    rows = (len(stems)+cols-1)//cols
    doc = fitz.open()
    page = doc.new_page(width=720, height=rows*290)
    # Pass the font as bytes: PyMuPDF rejects font paths containing non-ASCII characters.
    page.insert_font(fontname='dvbold', fontbuffer=Path(BOLD_FONT).read_bytes())
    for i,stem in enumerate(stems):
        src = fitz.open(NATIVE/f'{stem}-en.pdf')
        col,row = i%cols,i//cols
        if len(stems)%cols and i==len(stems)-1: col=.5
        rect = fitz.Rect(col*360,row*290+15,(col+1)*360,(row+1)*290)
        page.show_pdf_page(rect,src,0)
        page.insert_text((col*360+5,row*290+19),chr(65+i),fontsize=14,fontname='dvbold')
        src.close()
    # Save the PDF directly: a PDF-to-SVG-to-PDF round trip drops the soft mask of the
    # vendor credit badge and turns it into a black box.
    pdf_path = FIG/f'{name}.pdf'
    try:
        doc.subset_fonts()   # keep only the glyphs used by the panel letters
    except Exception:
        pass
    doc.save(str(pdf_path), garbage=3, deflate=True)
    subprocess.run(['pdftocairo', '-svg', str(pdf_path), str(FIG/f'{name}.svg')], check=True)
    with fitz.open(pdf_path) as final:
        final[0].get_pixmap(matrix=fitz.Matrix(2.5,2.5),alpha=False).save(FIG/f'{name}.png')
    doc.close()


def cassette():
    fig,ax=plt.subplots(figsize=(8,4.8))
    ax.set(xlim=(-160,750),ylim=(-.8,4.5));ax.axis('off')
    rows=[('WT',[]),('S211A',[211]),('S142A',[142]),('S142A/S211A',[142,211])]
    for i,(name,sites) in enumerate(rows):
        y=3.55-i
        ax.text(-145,y,name,va='center',fontsize=11,weight='bold')
        for start,end,color in ((0,476,TEAL),(476,483,NPG[4]),(483,722,NPG[1])):
            ax.add_patch(Rectangle((start,y-.14),end-start,.28,color=color,lw=0))
        for pos in sites:
            ax.plot([pos,pos],[y-.18,y+.18],color=ORANGE,lw=2)
            ax.text(pos,y+.26,f'S{pos}A',ha='center',fontsize=9,color=ORANGE)
    ax.text(238,4.2,'TFEB',ha='center',fontsize=11,color=TEAL,weight='bold')
    ax.text(602,4.2,'EGFP',ha='center',fontsize=11,color=NPG[1],weight='bold')
    ax.annotate('GDPPVAT',xy=(479.5,3.75),xytext=(479.5,4.25),ha='center',fontsize=8.5,
                arrowprops={'arrowstyle':'-','color':INK,'lw':.6})
    for pos in (1,142,211,476,483,722):
        ax.plot([pos,pos],[-.15,-.24],color=INK,lw=.6)
        dy=-.43 if pos==476 else (-.65 if pos==483 else -.43)
        ax.text(pos,dy,str(pos),ha='center',fontsize=8)
    ax.text(-145,-.45,'Residue',fontsize=9)
    fig.tight_layout(pad=.4);save_figure(fig,FIG/'fusion-cassettes')


def schedule():
    fig,ax=plt.subplots(figsize=(8.5,4.1));ax.set(xlim=(-44,125),ylim=(-.7,4.6));ax.axis('off')
    rows=[('Fed',3.7,BLUE,[(0,120)],[0,15,30,60,75,90,120]),
          ('AA starvation',2.7,TEAL,[(0,120)],[0,15,30,60,75,90,120]),
          ('Starve / refeed',1.7,TEAL,[(0,60)],[60,75,90,120]),
          ('Torin1',.7,ORANGE,[(0,60)],[0,15,30,60]),
          ('DMSO',-.3,NPG[5],[(0,60)],[0,15,30,60])]
    for name,y,color,segments,times in rows:
        ax.text(-42,y,name,va='center',fontsize=10.5)
        for start,end in segments:ax.add_patch(Rectangle((start,y-.08),end-start,.16,color=color))
        if 'refeed' in name:ax.add_patch(Rectangle((60,y-.08),60,.16,color=BLUE))
        for t in times:
            ax.plot(t,y,'o',color=INK,ms=3)
    for t in (0,15,30,60,75,90,120):ax.text(t,4.3,str(t),ha='center',fontsize=9)
    ax.text(-42,4.3,'Time (min)',fontsize=10)
    ax.annotate('Refeed',xy=(60,1.78),xytext=(75,2.08),fontsize=9,color=BLUE,
                arrowprops={'arrowstyle':'-','color':BLUE,'lw':.7})
    fig.tight_layout(pad=.3);save_figure(fig,FIG/'control-schedule')


def evidence():
    save_figure(evidence_flow(),FIG/'evidence-design')


def main():
    style()
    inputs=json.loads((NATIVE/'display-inputs.json').read_text())
    sources={str((NATIVE/'display-inputs.json').relative_to(ROOT)):sha(NATIVE/'display-inputs.json')}
    outputs=[]
    for row in inputs['records']:
        assert sha(ROOT/row['source'])==row['source_sha256']
        assert sha(ROOT/row['display_file'])==row['display_sha256']
        p=NATIVE/(row['title']+'-publication.svg')
        sources[str(p.relative_to(ROOT))]=sha(p)
        xml=p.read_text()
        assert xml.count('碱基对')==1
        colors={x.upper() for x in re.findall(r'(?:fill|stroke)="(#[0-9a-fA-F]{6})"',xml)}
        chromatic={x for x in colors if not (x[1:3]==x[3:5]==x[5:7])}
        assert chromatic <= set(NPG),chromatic-set(NPG)
        for f in row['features']:assert f['color'] in chromatic,(row['title'],f['name'])
        # Only translate the unit in a derived editable vector document.
        # Native paths, coordinates, colors and the vendor credit are untouched.
        localized=xml.replace('碱基对','bp')
        en=NATIVE/(row['title']+'-en.svg');en.write_text(localized)
        pdf=en.with_suffix('.pdf');cairosvg.svg2pdf(bytestring=localized.encode(),write_to=str(pdf))
        outputs.extend([en,pdf])
    compose_native(['pEGFP-N1','TFEB-WT-EGFP'],'control-plasmids')
    compose_native(['TFEB-S211A-EGFP','TFEB-S142A-EGFP','TFEB-S142A-S211A-EGFP'],'mutant-plasmids')
    compose_native(['TFEB-3A-EGFP','TFEB-3T-EGFP'],'research-tagged-plasmids')
    compose_native(['TFEB-WT-untagged','TFEB-3A-untagged','TFEB-3T-untagged'],'research-untagged-plasmids')
    cassette();save_figure(schedule_panel(),FIG/'control-schedule');evidence()
    save_figure(research_design_panel(),FIG/'research-design')
    names=['control-plasmids','mutant-plasmids','fusion-cassettes','control-schedule','evidence-design',
           'research-tagged-plasmids','research-untagged-plasmids','research-design']
    outputs += [FIG/f'{n}.{ext}' for n in names for ext in ('pdf','png','svg')]
    sources.update({str(p.relative_to(ROOT)):sha(p) for p in (Path(__file__).resolve(),Path(__file__).with_name('publication_style.py'),Path(__file__).with_name('figure_layouts.py'),Path(__file__).with_name('scientific_schematics.py'),Path(__file__).with_name('academic_panels.py'))})
    payload={'status':'rendered; actual page review required','palette_source':SOURCE,
             'sources':sources,'outputs':{str(p.relative_to(ROOT)):sha(p) for p in outputs},
             'native_renderer':'SnapGene Viewer 8.2.3, GUI SVG export, largest font, enzymes hidden',
             'postprocessing':'Translate unit to bp in derivative SVG; preserve paths, colors and SnapGene credit; compose panels',
             'limitations':['No experimental data; schematics describe planned comparisons','Selected native feature labels only; complete annotations remain in design GenBank files']}
    (FIG/'publication-figures.json').write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'figure_families':len(names),'native_maps':len(inputs['records'])}))


if __name__=='__main__':main()
