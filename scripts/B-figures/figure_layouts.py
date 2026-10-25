#!/usr/bin/env python3
"""Reusable vector figure layouts; caller supplies accepted design values."""

import textwrap
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from publication_style import INK, BLUE, TEAL, ORANGE


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "pdf.fonttype": 42,
        "svg.fonttype": "none", "figure.facecolor": "white",
        "axes.facecolor": "white", "text.color": INK,
        "font.size": 10,
    })


def save_figure(fig, prefix):
    """Save only explicit generated figure destinations; sources stay untouched."""
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(prefix.with_suffix("." + suffix), dpi=220,
                    bbox_inches="tight", facecolor="white")
    plt.close(fig)


def workflow(nodes, controls, title):
    """Group the seven accepted operations into four compact scientific stages."""
    if len(nodes) != len(controls) or len(nodes) < 2:
        raise ValueError("Each workflow node requires one control/decision note")
    if len(nodes)!=7:raise ValueError('Review stage grouping when operations change')
    from scientific_schematics import cloning_flow
    style()
    return cloning_flow()


def predicted_gel(panels, marker_sizes, lower_bp=75, upper_bp=10000,
                  marker_tick_labels=None):
    """Draw equal-width schematic bands; no intensity or gel-resolution claim.

    panels: [{title, lanes: [{label, bands_bp, note}]}]
    Small fragments are listed below the panel, never silently discarded.
    """
    if not panels or not marker_sizes:
        raise ValueError("Panels and an explicitly defined idealized marker are required")
    if any(size <= 0 for size in marker_sizes):
        raise ValueError("Marker sizes must be positive")
    style()
    columns = min(2, len(panels))
    rows = math.ceil(len(panels) / columns)
    fig, axes = plt.subplots(rows, columns, figsize=(5 * columns, 3.9 * rows), squeeze=False)
    for ax in axes.flat[len(panels):]:
        ax.axis("off")
    for ax, panel in zip(axes.flat, panels):
        lanes = [{"label": "M", "bands_bp": marker_sizes}, *panel["lanes"]]
        ax.set_yscale("log")
        ax.set_ylim(lower_bp, upper_bp)
        ax.set_xlim(-0.6, len(lanes) - 0.3)
        ax.set_xticks(range(len(lanes)), [x["label"] for x in lanes], fontsize=12)
        ax.xaxis.tick_top()
        tick_labels = marker_tick_labels or {x: f"{x:,}" for x in marker_sizes}
        ticks = [x for x in tick_labels if lower_bp <= x <= upper_bp]
        ax.set_yticks(ticks, [tick_labels[x] for x in ticks], fontsize=11)
        ax.minorticks_off()
        ax.set_ylabel("Fragment length (bp; logarithmic)", fontsize=11)
        ax.grid(axis="y", color="#DDDDDD", lw=0.6)
        ax.set_axisbelow(True)
        ax.set_title(panel["title"], fontsize=13, weight="bold", pad=30)
        omitted, notes = [], []
        for lane_index, lane in enumerate(lanes):
            for size in sorted(set(lane["bands_bp"]), reverse=True):
                if size <= 0:
                    raise ValueError("A fragment cannot have non-positive length")
                if lower_bp <= size <= upper_bp:
                    ax.plot([lane_index - 0.22, lane_index + 0.22], [size, size],
                            color=INK if lane_index == 0 else TEAL,
                            lw=2.7, solid_capstyle="butt")
                else:
                    omitted.append(f"{lane['label']}: {size} bp")
            if lane_index and lane.get("note"):
                notes.append(f"{lane['label']}: {lane['note']}")
        for spine in ax.spines.values():
            spine.set_color("#BBBBBB")
    fig.subplots_adjust(left=.08,right=.98,top=.87,bottom=.05,wspace=.32,hspace=.55)
    return fig
