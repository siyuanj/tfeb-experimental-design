#!/usr/bin/env python3
"""Draw the three fusion architectures (untagged, linker-only, EGFP fusion) from accepted lengths.

Residue and plasmid lengths are read from the sequence files, not typed in, so the
drawing cannot drift from the designs.
"""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO
from publication_style import INK, TEAL, NPG, SOURCE
from figure_layouts import save_figure, style

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / '报告/images/B'
LINKER_COLOR, EGFP_COLOR = NPG[4], NPG[1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def residues(path):
    return len(SeqIO.read(path, 'fasta').seq)


def main():
    style()
    dna_files = {
        'untagged': ROOT / '设计/research-v3/TFEB-WT-untagged.fasta',
        'linker': ROOT / '设计/linker-v4/TFEB-WT-linker.fasta',
        'EGFP': ROOT / '设计/construct.fasta',
    }
    protein_files = {
        'untagged': ROOT / '设计/research-v3/TFEB-WT-untagged-protein.fasta',
        'linker': ROOT / '设计/linker-v4/TFEB-WT-linker-protein.fasta',
        'EGFP': ROOT / '设计/fusion_protein.fasta',
    }
    aa = {k: residues(v) for k, v in protein_files.items()}
    bp = {k: residues(v) for k, v in dna_files.items()}
    assert (aa['untagged'], aa['linker'], aa['EGFP']) == (476, 483, 722), aa
    assert (bp['untagged'], bp['linker'], bp['EGFP']) == (5396, 5417, 6134), bp
    rows = [('Untagged', 'untagged', 3.0), ('Linker-only', 'linker', 2.0), ('EGFP fusion', 'EGFP', 1.0)]
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.set(xlim=(-190, 900), ylim=(0.0, 3.9))
    ax.axis('off')
    for label, key, y in rows:
        ax.text(-180, y + .08, label, va='center', fontsize=11, weight='bold')
        ax.text(-180, y - .22, f'{bp[key]:,} bp; {aa[key]} aa', va='center', fontsize=8.5)
        ax.add_patch(Rectangle((0, y - .14), 476, .28, color=TEAL, lw=0))
        end = 476
        if key in ('linker', 'EGFP'):
            ax.add_patch(Rectangle((476, y - .14), 7, .28, color=LINKER_COLOR, lw=0))
            end = 483
        if key == 'EGFP':
            ax.add_patch(Rectangle((483, y - .14), 239, .28, color=EGFP_COLOR, lw=0))
            end = 722
        ax.add_patch(Rectangle((end, y - .14), 9, .28, color=INK, lw=0))
    ax.text(238, 3.3, 'TFEB 1–476', ha='center', fontsize=10, color=TEAL, weight='bold')
    ax.text(602, 1.5, 'EGFP 484–722', ha='center', fontsize=10, color=EGFP_COLOR, weight='bold')
    ax.annotate('GDPPVAT 477–483', xy=(480, 2.14), xytext=(560, 2.55), fontsize=9,
                arrowprops={'arrowstyle': '-', 'color': INK, 'lw': .6})
    ax.plot([800, 800], [2.0, 3.0], color=INK, lw=.8)
    ax.text(812, 2.5, '+7 aa linker;\nstop 21 nt later', va='center', fontsize=8.5)
    ax.plot([800, 800], [1.0, 2.0], color=INK, lw=.8)
    ax.text(812, 1.5, '+239 aa EGFP', va='center', fontsize=8.5)
    for pos, dx, dy in ((1, 0, 0), (476, -14, 0), (483, 16, -.2), (722, 0, 0)):
        ax.plot([pos, pos], [0.62, 0.7], color=INK, lw=.6)
        ax.text(pos + dx, 0.45 + dy, str(pos), ha='center', fontsize=8)
    ax.text(-180, 0.45, 'Residue', fontsize=9)
    fig.tight_layout(pad=.4)
    save_figure(fig, FIG / 'architecture-panel')
    outputs = {str(FIG.joinpath('architecture-panel' + s).relative_to(ROOT)): sha(FIG / ('architecture-panel' + s))
               for s in ('.pdf', '.svg', '.png')}
    record = {'status': 'rendered; actual page review required', 'palette_source': SOURCE,
              'sources': {str(p.relative_to(ROOT)): sha(p) for p in [*dna_files.values(), *protein_files.values()]},
              'script': {str(Path(__file__).resolve().relative_to(ROOT)): sha(Path(__file__).resolve())},
              'outputs': outputs, 'values': {'aa': aa, 'bp': bp}}
    (FIG / 'architecture-panel.json').write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(record['values']))


if __name__ == '__main__':
    main()
