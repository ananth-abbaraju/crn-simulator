#!/usr/bin/env python3
"""Plot crnsim output.

    python/plot.py data/lv.csv                            # trajectory
    python/plot.py trajectory data/lv.csv --ode data/lv_ode.csv -o figures/lv.png
    python/plot.py phase data/lv.csv --ode data/lv_ode.csv
"""

import argparse
import os
import sys

import matplotlib
import pandas as pd

import style

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


def read_traj(path):
    df = pd.read_csv(path)
    if "t" not in df.columns:
        sys.exit(f"{path}: expected a 't' column, got {list(df.columns)}")
    return df, [c for c in df.columns if c != "t"]


def finish(fig, out, show):
    if show:
        matplotlib.use("MacOSX", force=True)
        plt.show()
        return
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out)
    print(f"wrote {out}")


def headroom(ax, factor=1.30):
    """Leave clear space at the top for the legends, so they never sit on a curve."""
    top = ax.get_ylim()[1]
    ax.set_ylim(0, top * factor)


def method_legend(ax, has_ode, loc="upper right"):
    """A second legend keyed on line style, drawn in neutral ink.

    Identity (which species) is carried by the color legend; this one carries
    which description of the network produced the curve. Keeping them separate
    is what stops a four-entry legend from reading as four species.
    """
    if not has_ode:
        return
    handles = [
        Line2D([], [], color=style.INK_2, lw=1.0, label="stochastic (Gillespie SSA)"),
        Line2D([], [], color=style.INK_2, lw=1.8, ls=(0, (5, 2)), label="deterministic (RK4)"),
    ]
    leg = ax.legend(handles=handles, loc=loc, fontsize=8, labelcolor=style.INK_2)
    ax.add_artist(leg)


def cmd_trajectory(args):
    df, species = read_traj(args.csv)
    fig, ax = plt.subplots(figsize=(8.2, 4.0))

    for i, sp in enumerate(species):
        # Counts jump discretely at reaction events and are constant between
        # them, so the SSA trace is a step function, not a polyline.
        ax.step(df["t"], df[sp], where="post", color=style.color_for(sp, i),
                lw=1.0, alpha=0.9, label=sp, zorder=3)

    if args.ode:
        ode, ode_species = read_traj(args.ode)
        for i, sp in enumerate(ode_species):
            if sp not in species:
                continue
            ax.plot(ode["t"], ode[sp], color=style.color_for(sp, species.index(sp)),
                    lw=1.6, ls=(0, (5, 2)), alpha=0.95, zorder=4)

    ax.set_xlabel("time")
    ax.set_ylabel("molecule count")
    ax.set_xlim(df["t"].min(), df["t"].max())
    ax.set_ylim(bottom=0)
    headroom(ax)
    ax.set_title(args.title or "Molecule counts over time")
    style.subtitle(ax, args.subtitle or os.path.basename(args.csv))

    leg = ax.legend(loc="upper left", ncol=len(species), labelcolor=style.INK_2)
    ax.add_artist(leg)
    method_legend(ax, bool(args.ode))

    finish(fig, args.out, args.show)


def cmd_phase(args):
    df, species = read_traj(args.csv)
    if len(species) < 2:
        sys.exit("phase portrait needs at least two species")
    xs, ys = args.x or species[0], args.y or species[1]

    fig, ax = plt.subplots(figsize=(5.4, 5.0))

    if args.ode:
        ode, _ = read_traj(args.ode)
        ax.plot(ode[xs], ode[ys], color=style.SERIES[1], lw=1.8, ls=(0, (5, 2)),
                zorder=3, label="deterministic (RK4)")

    ax.plot(df[xs], df[ys], color=style.SERIES[0], lw=0.7, alpha=0.75, zorder=2,
            label="stochastic (Gillespie SSA)")

    # Direct-label the start: a phase portrait has no time axis, so without it
    # the reader cannot tell where the orbit begins or which way it turns.
    ax.plot(df[xs].iloc[0], df[ys].iloc[0], "o", ms=7, color=style.INK,
            mec=style.SURFACE, mew=2, zorder=5)
    ax.annotate("start", (df[xs].iloc[0], df[ys].iloc[0]),
                textcoords="offset points", xytext=(9, -13),
                fontsize=8.5, color=style.INK_2, zorder=5)

    ax.set_xlabel(f"{xs} count")
    ax.set_ylabel(f"{ys} count")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    headroom(ax, 1.18)
    ax.set_title(args.title or f"Phase portrait: {ys} vs {xs}")
    style.subtitle(ax, args.subtitle or "the deterministic orbit closes; the stochastic path spirals out")
    ax.legend(loc="upper left", labelcolor=style.INK_2)

    finish(fig, args.out, args.show)


def main():
    style.apply()

    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp, default_out):
        sp.add_argument("-o", "--out", default=default_out)
        sp.add_argument("--show", action="store_true", help="display instead of saving")
        sp.add_argument("--title")
        sp.add_argument("--subtitle")

    t = sub.add_parser("trajectory", help="counts vs time, SSA with optional ODE overlay")
    t.add_argument("csv")
    t.add_argument("--ode")
    common(t, "figures/trajectory.png")
    t.set_defaults(func=cmd_trajectory)

    ph = sub.add_parser("phase", help="one species against another")
    ph.add_argument("csv")
    ph.add_argument("--ode")
    ph.add_argument("--x")
    ph.add_argument("--y")
    common(ph, "figures/phase.png")
    ph.set_defaults(func=cmd_phase)

    # `plot.py data/lv.csv` -- the bare form, with `trajectory` implied.
    argv = sys.argv[1:]
    if argv and argv[0] not in sub.choices and not argv[0].startswith("-"):
        argv.insert(0, "trajectory")

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
