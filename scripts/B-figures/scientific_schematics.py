"""Original vector scientific objects and compact figure layouts.

Bioicons and D2 informed the design approach; no third-party icon is copied.
These are conceptual objects, never measured microscopy or sequencing traces.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, Arc, Rectangle, FancyArrowPatch
from publication_style import INK, BLUE, TEAL, ORANGE, NPG

def arrow(ax,a,b,color=INK):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=10,lw=.9,color=color))

def dna(ax,x,y,s=1,color=TEAL):
    t=np.linspace(-.55,.55,100)
    for sign in (-1,1):ax.plot(x+t*s,y+sign*.13*np.sin(t*8)*s,color=color,lw=1.6)
    for t0 in np.linspace(-.5,.5,11):
        d=.13*np.sin(t0*8)*s;ax.plot([x+t0*s]*2,[y-d,y+d],color=color,lw=.65)

def plasmid(ax,x,y,s=1):
    ax.add_patch(Circle((x,y),.35*s,ec=BLUE,fc='white',lw=1.25))
    ax.add_patch(Arc((x,y),.7*s,.7*s,theta1=15,theta2=150,color=TEAL,lw=5))
    ax.add_patch(Arc((x,y),.7*s,.7*s,theta1=230,theta2=290,color=NPG[1],lw=4))

def dish(ax,x,y,s=1):
    ax.add_patch(Ellipse((x,y),.95*s,.50*s,fc='white',ec=BLUE,lw=1.3))
    ax.add_patch(Arc((x,y-.09*s),.95*s,.50*s,theta1=180,theta2=360,color=BLUE,lw=1.3))
    for dx,dy in [(-.24,.04),(-.05,.12),(.20,.06),(.03,-.13),(.28,-.07),(-.23,-.11)]:
        ax.add_patch(Circle((x+dx*s,y+dy*s),.025*s,color=TEAL))

def cell(ax,x,y,s=1):
    ax.add_patch(Ellipse((x,y),.94*s,.67*s,ec=BLUE,fc='white',lw=1.3))
    ax.add_patch(Ellipse((x-.08*s,y),.39*s,.39*s,ec=TEAL,fc=NPG[6],lw=1.1))
    for dx,dy in [(-.12,.06),(-.05,-.1),(.05,.02),(.24,-.09),(-.30,.06)]:
        ax.add_patch(Circle((x+dx*s,y+dy*s),.022*s,color=TEAL))

def sequence(ax,x,y,s=1):
    for yy in (-.12,.12):
        ax.plot([x-.48*s,x+.48*s],[y+yy*s]*2,color=BLUE,lw=1)
        for i in range(8):ax.plot([x+(-.42+i*.12)*s]*2,[y+(yy-.04)*s,y+(yy+.04)*s],color=TEAL,lw=2)
    arrow(ax,(x-.48*s,y+.3*s),(x+.45*s,y+.3*s),TEAL)
    arrow(ax,(x+.48*s,y-.3*s),(x-.45*s,y-.3*s),TEAL)

def transcript(ax,x,y,s=1):
    dna(ax,x,y+.18*s,s*.8)
    for j in range(3):
        xx=np.linspace(-.25,.30,40)
        ax.plot(x+xx*s,y+(-.07-j*.12+.025*np.sin(xx*25))*s,color=BLUE,lw=1.3)

def cloning_flow():
    fig,ax=plt.subplots(figsize=(8,2.45));ax.set(xlim=(0,8),ylim=(0,2.45));ax.axis('off')
    centers=[1,3,5,7]
    labels=[('Prepare insert','PCR + purification'),('Assemble plasmid','Digest + ligation'),
            ('Recover clones','Selection + screening'),('Verify identity','Bidirectional sequence')]
    glyphs=[dna,plasmid,dish,sequence]
    for i,(x,(title,sub),glyph) in enumerate(zip(centers,labels,glyphs)):
        ax.text(x-.78,2.13,f'0{i+1}',fontsize=9.5,color=TEAL,weight='bold')
        ax.plot([x-.45,x+.78],[2.17]*2,color=NPG[6],lw=1)
        glyph(ax,x,1.45,1.1)
        ax.text(x,.65,title,fontsize=11,weight='bold',ha='center')
        ax.text(x,.32,sub,fontsize=9,ha='center')
        if i<3:arrow(ax,(x+.67,1.45),(x+1.30,1.45))
    fig.tight_layout(pad=.15);return fig

def evidence_flow():
    fig,ax=plt.subplots(figsize=(8,3.45));ax.set(xlim=(0,8),ylim=(0,3.45));ax.axis('off')
    groups=[(1.3,plasmid,'Identity',['Sequence','Intact protein']),
            (4,cell,'Cellular response',['Localization + recovery','Site / pathway controls']),
            (6.7,transcript,'Functional output',['TFEB loss + rescue','Transcription + turnover'])]
    for x,glyph,title,sub in groups:
        glyph(ax,x,2.35,1.25)
        ax.text(x,1.55,title,ha='center',fontsize=12,weight='bold',color=BLUE)
        ax.plot([x-.98,x+.98],[1.27]*2,color=NPG[6],lw=1)
        ax.text(x,.88,sub[0],ha='center',fontsize=9.5)
        ax.text(x,.49,sub[1],ha='center',fontsize=9.5)
    for x in (1.3,4):arrow(ax,(x+.82,2.35),(x+1.83,2.35))
    fig.tight_layout(pad=.15);return fig

def research_matrix():
    fig,ax=plt.subplots(figsize=(8.4,4.8));ax.set(xlim=(0,8.4),ylim=(0,4.8));ax.axis('off')
    ax.text(.12,4.52,'A',fontsize=12,weight='bold');ax.text(.4,4.52,'C-terminal design',fontsize=11,weight='bold')
    ax.text(4.75,4.52,'B',fontsize=12,weight='bold');ax.text(5.04,4.52,'Matched conditions',fontsize=11,weight='bold')
    for x,t in [(5.35,'Vehicle'),(6.55,'Torin1'),(7.75,'MK-8722')]:ax.text(x,3.93,t,ha='center',fontsize=9)
    rows=[('WT',True),('3A',True),('3T',True),('WT',False),('3A',False),('3T',False)]
    for i,(variant,tag) in enumerate(rows):
        y=3.48-i*.43
        ax.text(.18,y,variant,va='center',fontsize=9.5,weight='bold')
        ax.add_patch(Rectangle((.75,y-.07),2.12,.14,fc=TEAL,lw=0))
        if tag:ax.add_patch(Rectangle((2.88,y-.07),1.06,.14,fc=NPG[1],lw=0))
        if variant!='WT':
            ax.plot([2.80]*2,[y-.14,y+.14],color=ORANGE,lw=1.6)
        for x in (5.35,6.55,7.75):ax.add_patch(Circle((x,y),.067,fc='white',ec=BLUE,lw=1.25))
    ax.text(1.8,3.92,'TFEB',ha='center',fontsize=9,color=TEAL,weight='bold')
    ax.text(3.4,3.92,'EGFP',ha='center',fontsize=9,color=NPG[1],weight='bold')
    ax.plot([.15,8.1],[2.4]*2,color=NPG[6],lw=.6,ls=':')
    ax.plot([.15,8.1],[.95]*2,color=NPG[6],lw=.8)
    for x,t in [(1.45,'Localization'),(4.1,'Transcription'),(6.85,'Phosphorylation')]:
        ax.text(x,.47,t,ha='center',fontsize=10.5,color=BLUE,weight='bold')
    fig.tight_layout(pad=.2);return fig
