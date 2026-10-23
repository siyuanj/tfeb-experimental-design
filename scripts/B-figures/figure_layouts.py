#!/usr/bin/env python3
"""Reusable vector figure layouts; caller supplies accepted design values."""

import textwrap
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK = "#243849"
BLUE = "#326EB5"
TEAL = "#148C83"
ORANGE = "#B65B39"


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
    """Draw a top-to-bottom workflow with aligned control/decision annotations."""
    if len(nodes) != len(controls) or len(nodes) < 2:
        raise ValueError("Each workflow node requires one control/decision note")
    style()
    fig, ax = plt.subplots(figsize=(10, 1.08 * len(nodes) + 0.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.4, len(nodes) + 0.6)
    ax.axis("off")
    ax.text(0.1, len(nodes) + 0.33, title, fontsize=13, weight="bold")
    for index, ((heading, detail), note) in enumerate(zip(nodes, controls)):
        y = len(nodes) - index - 0.45
        color = TEAL if index == len(nodes) - 1 else BLUE
        ax.add_patch(FancyBboxPatch((0.1, y - 0.36), 5.75, 0.77,
                     boxstyle="round,pad=0.03,rounding_size=0.06",
                     fc="#F1F6FA", ec=color, lw=1.2))
        ax.text(0.28, y + 0.19, heading, fontsize=11.5, weight="bold", va="center")
        ax.text(0.28, y - 0.12, detail, fontsize=10.5, va="center")
        ax.text(6.14, y + 0.02, textwrap.fill(note, width=41),
                fontsize=10.5, va="center", linespacing=1.3)
        if index < len(nodes) - 1:
            ax.add_patch(FancyArrowPatch((2.95, y - 0.39), (2.95, y - 0.57),
                         arrowstyle="-|>", mutation_scale=12, color=INK, lw=1))
    fig.tight_layout()
    return fig


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
    fig, axes = plt.subplots(rows, columns, figsize=(6 * columns, 5.2 * rows), squeeze=False)
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
        ax.grid(axis="y", color="#E2E8EC", lw=0.6)
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
        footer = "\n".join(textwrap.fill(note, width=63) for note in notes)
        if omitted:
            footer += "\nBelow/above displayed scale: " + "; ".join(omitted) + "."
        ax.text(0, -0.06, footer, transform=ax.transAxes, ha="left", va="top",
                fontsize=11.5, linespacing=1.25)
        for spine in ax.spines.values():
            spine.set_color("#B8C6D0")
    fig.suptitle("PREDICTED BAND POSITIONS - SCHEMATIC, NO EXPERIMENTAL DATA",
                 fontsize=13, weight="bold", color=ORANGE, y=0.98)
    fig.subplots_adjust(left=0.07, right=0.98, top=0.87, bottom=0.17,
                        wspace=0.27, hspace=1.10 if rows > 1 else 0.2)
    return fig
