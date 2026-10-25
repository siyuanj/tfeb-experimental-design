"""Journal-style panels drawn by code from declared design values.

Both figures describe planned experiments; they contain no measured data.
Contrast weights are derived from the written definitions, not typed in.
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle
from publication_style import INK, BLUE, TEAL, ORANGE, NPG

LINKER, EGFP, MUTED = NPG[4], NPG[1], NPG[5]
RULE = '#BDBDBD'


def schedule_panel():
    """Sampling schedule: bars show exposure, open circles show separate matched wells."""
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.set(xlim=(-50, 124), ylim=(-0.75, 5.7))
    for side in ('top', 'right', 'left'):
        ax.spines[side].set_visible(False)
    ax.spines['bottom'].set_position(('data', -0.25))
    ax.spines['bottom'].set_bounds(0, 120)
    ax.set_yticks([])
    ax.set_xticks([0, 15, 30, 60, 75, 90, 120])
    ax.tick_params(axis='x', labelsize=9)
    ax.set_xlabel('Time (min)', fontsize=10)
    ax.xaxis.set_label_coords((60 + 50) / 174, -0.045)
    rows = [  # label, y, segments (start, end, color), sampling times, group
        ('AA-replete (fed)', 4.55, [(0, 120, BLUE)], [0, 15, 30, 60, 75, 90, 120]),
        ('AA-depleted', 3.70, [(0, 120, TEAL)], [0, 15, 30, 60, 75, 90, 120]),
        ('Starve 60 min, refeed', 2.85, [(0, 60, TEAL), (60, 120, BLUE)], [60, 75, 90, 120]),
        ('Torin1', 1.75, [(0, 60, ORANGE)], [0, 15, 30, 60]),
        ('DMSO (matched)', 0.90, [(0, 60, MUTED)], [0, 15, 30, 60]),
    ]
    for label, y, segments, times in rows:
        ax.text(-3, y, label, ha='right', va='center', fontsize=9.5)
        for start, end, color in segments:
            ax.add_patch(Rectangle((start, y - .17), end - start, .34, fc=color, ec='none'))
        ax.plot(times, [y] * len(times), 'o', mfc='white', mec=INK, ms=5.2, mew=1.0, zorder=3)
    ax.plot([-50, 124], [2.3] * 2, color=RULE, lw=.6)
    ax.text(-49, 5.15, 'Nutrient', fontsize=10, weight='bold', va='center')
    ax.text(-49, 2.05, 'Drug', fontsize=10, weight='bold', va='center')
    ax.plot([60, 60], [-.1, 5.0], color=INK, lw=.7, ls=(0, (3, 3)), zorder=1)
    for x0, x1, text in ((0, 60, 'Acute series'), (60, 120, 'Recovery series')):
        ax.plot([x0 + 1, x1 - 1], [5.3, 5.3], color=INK, lw=.8)
        ax.text((x0 + x1) / 2, 5.5, text, ha='center', va='bottom', fontsize=9.5)
    handles = [Rectangle((0, 0), 1, 1, fc=c, ec='none') for c in (BLUE, TEAL, ORANGE, MUTED)]
    handles.append(Line2D([], [], marker='o', ls='', mfc='white', mec=INK, ms=5.2))
    fig.legend(handles, ['AA-replete', 'AA-depleted', 'Torin1', 'DMSO', 'Sampling time (separate matched wells)'],
               loc='lower center', ncol=5, frameon=False, fontsize=8.2, handlelength=1.4, columnspacing=1.2)
    fig.subplots_adjust(left=.02, right=.99, top=.97, bottom=.27)
    return fig


def _weights():
    """Weights of A_M, A_T and C on the (genotype, architecture) rows of each stimulus."""
    cells = {}
    for stimulus in ('M', 'T'):
        for genotype, g_sign in (('3A', 1), ('WT', -1)):
            for architecture, a_sign in (('E', 1), ('U', -1)):
                cells[(genotype, architecture, stimulus)] = g_sign * a_sign
    return cells  # A_s weight; C weight is +A_M on M cells and -A_T on T cells


def _glyph(ax, x, y, sign):
    ax.add_patch(Circle((x, y), .092, fc=INK if sign > 0 else 'white', ec=INK, lw=1.0, zorder=3))
    ax.add_patch(Rectangle((x - .05, y - .008), .10, .016, fc='white' if sign > 0 else INK, ec='none', zorder=4))
    if sign > 0:
        ax.add_patch(Rectangle((x - .008, y - .05), .016, .10, fc='white', ec='none', zorder=4))


def research_design_panel():
    """Constructs, matched conditions and the weights of the three primary contrasts."""
    fig = plt.figure(figsize=(7.2, 4.9))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 7.2), ylim=(0, 4.9))
    ax.axis('off')
    scale = 1.75 / 722                      # inches per residue
    rows = [  # group, genotype, architecture key, y
        ('EGFP fusion', 'WT', 'E', 3.55), ('EGFP fusion', '3A', 'E', 3.15), ('EGFP fusion', '3T', 'E', 2.75),
        ('Untagged', 'WT', 'U', 2.20), ('Untagged', '3A', 'U', 1.80), ('Untagged', '3T', 'U', 1.40),
        ('Linker-only', 'WT', 'L', 0.85), ('Linker-only', '3A', 'L', 0.45),
    ]
    x0 = 1.0
    xs = {'Vehicle': 3.35, 'Torin1': 3.90, 'MK-8722': 4.45}
    for letter, x, title in (('A', .05, 'Constructs'), ('B', 3.0, 'Conditions'), ('C', 5.1, 'Contrast weights')):
        ax.text(x, 4.78, letter, fontsize=12, weight='bold', va='top')
        ax.text(x + .26, 4.78, title, fontsize=10.5, weight='bold', va='top')
    for name, x in xs.items():
        ax.text(x, 4.30, name, ha='center', va='center', fontsize=8.6, rotation=0)
    weights = _weights()
    columns = [('A$_M$', 5.28), ('A$_T$', 5.78), ('C', 6.28), ('C', 6.78)]
    ax.text(5.53, 4.30, 'by stimulus', ha='center', va='center', fontsize=8.6, style='italic')
    ax.text(6.53, 4.30, '$C=A_M-A_T$', ha='center', va='center', fontsize=8.6)
    for label, x in columns:
        ax.text(x, 4.02, label, ha='center', va='center', fontsize=9, weight='bold')
    for label, x in zip(('', '', 'M', 'T'), [c[1] for c in columns]):
        if label:
            ax.text(x, 3.84, label, ha='center', va='center', fontsize=8)
    for group, genotype, arch, y in rows:
        ax.text(.52, y, genotype, ha='center', va='center', fontsize=9.5, weight='bold')
        length = {'E': 722, 'U': 476, 'L': 483}[arch]
        ax.add_patch(Rectangle((x0, y - .075), 476 * scale, .15, fc=TEAL, ec='none'))
        if arch in ('E', 'L'):
            ax.add_patch(Rectangle((x0 + 476 * scale, y - .075), 7 * scale, .15, fc=LINKER, ec='none'))
        if arch == 'E':
            ax.add_patch(Rectangle((x0 + 483 * scale, y - .075), 239 * scale, .15, fc=EGFP, ec='none'))
        ax.add_patch(Rectangle((x0 + length * scale, y - .075), .018, .15, fc=INK, ec='none'))
        if genotype != 'WT':
            cluster = x0 + 467 * scale
            ax.plot([cluster] * 2, [y - .13, y + .13], color=ORANGE, lw=1.7)
        for name, x in xs.items():
            if arch == 'L':
                ax.add_patch(Circle((x, y), .075, fc='white', ec=MUTED, lw=1.2, ls=(0, (2, 1.4))))
            else:
                ax.add_patch(Circle((x, y), .075, fc=BLUE, ec=BLUE, lw=1.0))
        if arch != 'L' and genotype != '3T':
            for index, (_, x) in enumerate(columns):
                stimulus = 'M' if index in (0, 2) else 'T'
                if index < 2:
                    sign = weights[(genotype, arch, 'M' if index == 0 else 'T')]
                else:
                    base = weights[(genotype, arch, stimulus)]
                    sign = base if stimulus == 'M' else -base
                _glyph(ax, x, y, sign)
    for group, ys in (('EGFP fusion', (2.75, 3.55)), ('Untagged', (1.40, 2.20)), ('Linker-only\n(R06)', (.45, .85))):
        ax.text(.14, sum(ys) / 2, group, rotation=90, ha='center', va='center', fontsize=8.8, weight='bold')
    for y in (2.48, 1.13):
        ax.plot([.05, 7.15], [y, y], color=RULE, lw=.6)
    ax.add_patch(Rectangle((2.98, .22), 1.75, 4.0, fc='none', ec=RULE, lw=.6))
    ax.text(x0 + 476 * scale / 2, 3.98, 'TFEB', ha='center', va='center', fontsize=8.6, color=TEAL, weight='bold')
    ax.text(x0 + (483 + 239 / 2) * scale, 3.98, 'EGFP', ha='center', va='center', fontsize=8.6, color=EGFP, weight='bold')
    ax.text(x0 + 467 * scale, 4.20, 'S466/S467/S469', ha='center', va='center', fontsize=7.8, color=ORANGE)
    return fig
