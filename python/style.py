"""Shared plotting style.

One palette, defined once. Colors are assigned to *species* (the entity), never
to the method that produced the curve -- the SSA and the ODE for the same
species share a hue and are told apart by line style.

Palette slots 1-3 of the reference categorical theme; validated for
colorblind separation on a light surface (worst all-pairs CVD dE 9.2,
normal-vision dE 24.0).
"""

import matplotlib as mpl

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8984"
GRID = "#e8e7e3"

# Categorical slots, in fixed order. Species take them in the order they appear
# in the CSV, so a species keeps its hue across every figure in the project.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]


def color_for(species, index):
    """Stable hue per species. Index is position in the network's species list."""
    return SERIES[index % len(SERIES)]


def apply():
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "text.color": INK,
        "axes.labelcolor": INK_2,
        "axes.titlecolor": INK,
        "axes.titlesize": 11,
        "axes.titleweight": "medium",
        "axes.titlelocation": "left",
        "axes.titlepad": 22,
        "axes.labelsize": 9,
        "axes.edgecolor": GRID,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        # Hairline solid grid, one shade off the surface. Never dashed.
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "grid.linestyle": "-",
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "legend.labelcolor": INK_2,
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.25,
    })


def subtitle(ax, text):
    """One line of context under the title, in secondary ink."""
    ax.text(0.0, 1.02, text, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.5, color=INK_2)
