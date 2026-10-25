#!/usr/bin/env python3
"""Plot crnsim output.

    python/plot.py data/lv.csv                            # trajectory
    python/plot.py trajectory data/lv.csv --ode data/lv_ode.csv -o figures/lv.png
    python/plot.py phase data/lv.csv --ode data/lv_ode.csv
    python/plot.py ensemble data/lv_ensemble.csv --species X
    python/plot.py stats data/lv_stats.csv
    python/plot.py scale --pair 1 data/lv.csv data/lv_ode.csv \
                         --pair 25 data/lv_s25.csv data/lv_s25_ode.csv
    python/plot.py majority --point 100 90 data/am_split100.csv
"""

import argparse
import os
import sys

import matplotlib
import numpy as np
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


def cmd_ensemble(args):
    df = pd.read_csv(args.csv)
    col = f"first_zero_{args.species}"
    if col not in df.columns:
        sys.exit(f"{args.csv}: no column '{col}' (have {list(df.columns)})")

    times = df[col].dropna().to_numpy()
    n_total = len(df)
    n_extinct = len(times)
    if n_extinct == 0:
        sys.exit(f"no run in {args.csv} ever reached {args.species} = 0")

    # The runs missing from the histogram did not simply "survive" -- say which.
    rest = df[df[col].isna()]
    reasons = []
    if "stop" in df.columns:
        labels = {"horizon": "still oscillating at the horizon",
                  "population-cap": f"{args.species} blew up after its predators died out",
                  "absorbing": "absorbed without reaching zero",
                  "event-cap": "hit the event budget"}
        for reason, n in rest["stop"].value_counts().items():
            reasons.append(f"{n} {labels.get(reason, reason)}")

    fig, ax = plt.subplots(figsize=(7.6, 4.0))
    # One series: slot 1, no legend box -- the title names it.
    ax.hist(times, bins=args.bins, color=style.SERIES[0], alpha=0.9,
            edgecolor=style.SURFACE, linewidth=1.0, zorder=3)

    median = float(np.median(times))
    ax.axvline(median, color=style.INK, lw=1.2, zorder=4)
    ax.annotate(f"median {median:.1f}", (median, ax.get_ylim()[1]),
                textcoords="offset points", xytext=(6, -12),
                fontsize=8.5, color=style.INK, zorder=5)

    ax.set_xlabel(f"time of first {args.species} = 0")
    ax.set_ylabel("runs")
    ax.set_xlim(left=0)
    ax.set_title(args.title or f"{args.species} extinction time over {n_total} trajectories")
    tail = ("; " + ", ".join(reasons)) if reasons else ""
    style.subtitle(ax, args.subtitle or
                   f"{n_extinct}/{n_total} runs reached {args.species} = 0{tail}")

    finish(fig, args.out, args.show)


def cmd_stats(args):
    df = pd.read_csv(args.csv)
    species = [c[len("mean_"):] for c in df.columns if c.startswith("mean_")]
    if not species:
        sys.exit(f"{args.csv}: no mean_* columns")

    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    for i, sp in enumerate(species):
        c = style.color_for(sp, i)
        m, sd = df[f"mean_{sp}"], df[f"sd_{sp}"]
        ax.fill_between(df["t"], m - sd, m + sd, color=c, alpha=0.16, lw=0, zorder=2)
        ax.plot(df["t"], m, color=c, lw=1.8, label=sp, zorder=3)

    if args.ode:
        ode, ode_species = read_traj(args.ode)
        for sp in ode_species:
            if sp not in species:
                continue
            ax.plot(ode["t"], ode[sp], color=style.color_for(sp, species.index(sp)),
                    lw=1.4, ls=(0, (5, 2)), alpha=0.9, zorder=4)

    ax.set_xlabel("time")
    ax.set_ylabel("molecule count")
    ax.set_xlim(df["t"].min(), df["t"].max())
    ax.set_ylim(bottom=0)
    headroom(ax)
    ax.set_title(args.title or "Ensemble mean, one standard deviation shaded")
    style.subtitle(ax, args.subtitle or
                   "the ensemble mean is not a trajectory: independent orbits dephase and the mean flattens")
    leg = ax.legend(loc="upper left", ncol=len(species), labelcolor=style.INK_2)
    ax.add_artist(leg)
    method_legend(ax, bool(args.ode))

    finish(fig, args.out, args.show)


def cmd_majority(args):
    """Does the network actually compute majority?

    One point per initial split: the fraction of independent trajectories that
    locked in on X. Plotted against what an unbiased coin weighted by the
    initial population would give, which is the thing the network has to beat.
    """
    xs, ys, err = [], [], []
    for x0_s, y0_s, path in args.point:
        x0, y0 = int(x0_s), int(y0_s)
        df = pd.read_csv(path)
        wins = (df["final_X"] > df["final_Y"]).to_numpy()
        n = len(df)
        p = wins.mean()
        xs.append(x0 / (x0 + y0))
        ys.append(p)
        # Binomial standard error on the estimated probability.
        err.append(np.sqrt(max(p * (1.0 - p), 0.0) / n))

    order = np.argsort(xs)
    xs = np.asarray(xs)[order]
    ys = np.asarray(ys)[order]
    err = np.asarray(err)[order]

    fig, ax = plt.subplots(figsize=(6.6, 4.2))

    lo, hi = float(xs.min()), float(xs.max())
    grid = np.linspace(lo, hi, 200)
    ax.plot(grid, grid, color=style.INK_MUTED, lw=1.2, ls=(0, (5, 2)), zorder=2)
    # Exact majority as two horizontal segments rather than a connected step:
    # the riser at 0.5 is not data, and drawn it reads as a second axis.
    ax.hlines(0.0, lo, 0.5, color=style.INK_MUTED, lw=1.2, zorder=2)
    ax.hlines(1.0, 0.5, hi, color=style.INK_MUTED, lw=1.2, zorder=2)
    ax.annotate("exact majority", (hi, 1.0), textcoords="offset points",
                xytext=(-4, -14), ha="right", fontsize=8, color=style.INK_2)
    ax.annotate("weighted coin", (hi, hi), textcoords="offset points",
                xytext=(-4, 6), ha="right", fontsize=8, color=style.INK_2)

    ax.errorbar(xs, ys, yerr=err, color=style.SERIES[0], lw=1.8, marker="o",
                ms=6, mec=style.SURFACE, mew=1.5, capsize=3, elinewidth=1.0,
                zorder=4)

    ax.set_xlabel("initial fraction of X")
    ax.set_ylabel("runs reaching X consensus")
    ax.set_ylim(0, 1.05)
    ax.set_title(args.title or "Approximate majority amplifies a small initial lead")
    style.subtitle(ax, args.subtitle or
                   f"{len(df)} trajectories per point; bars are binomial standard error")

    finish(fig, args.out, args.show)


def cmd_scale(args):
    """Small multiples across volume scalings: the deterministic limit arriving."""
    pairs = args.pair
    n = len(pairs)
    fig, axes = plt.subplots(1, n, figsize=(4.0 * n, 3.8), sharey=True)
    if n == 1:
        axes = [axes]

    for ax, (omega_s, ssa_path, ode_path) in zip(axes, pairs):
        omega = float(omega_s)
        ssa, species = read_traj(ssa_path)
        ode, _ = read_traj(ode_path)
        keep = [sp for sp in species if not args.species or sp in args.species]

        for sp in keep:
            c = style.color_for(sp, species.index(sp))
            # Divide counts by the volume factor so every panel is in the same
            # units. The deterministic curve is then identical across panels and
            # only the noise around it changes -- which is the whole point.
            ax.step(ssa["t"], ssa[sp] / omega, where="post", color=c, lw=0.9,
                    alpha=0.9, label=sp, zorder=3)
            ax.plot(ode["t"], ode[sp] / omega, color=c, lw=1.6, ls=(0, (5, 2)),
                    alpha=0.95, zorder=4)

        ax.set_xlabel("time")
        ax.set_xlim(ssa["t"].min(), ssa["t"].max())
        ax.set_ylim(bottom=0)
        ax.set_title(f"Ω = {omega:g}", fontsize=10)
        style.subtitle(ax, f"{int(round(ssa[keep[0]].iloc[0]))} {keep[0]} at t = 0")

    for ax in axes:
        headroom(ax, 1.22)

    label = keep[0] if len(keep) == 1 else "species"
    axes[0].set_ylabel(f"{label} concentration (count / Ω)")
    # One series per panel needs no color legend -- the axis label names it.
    if len(keep) > 1:
        leg = axes[0].legend(loc="upper left", labelcolor=style.INK_2)
        axes[0].add_artist(leg)
    method_legend(axes[-1], True, loc="upper right")
    # Clear of the panel titles, which sit above the axes, not inside them.
    fig.suptitle(args.title or
                 "The deterministic limit: same network, more molecules",
                 x=0.0, y=1.06, ha="left", va="bottom", fontsize=12, color=style.INK)

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

    e = sub.add_parser("ensemble", help="extinction-time histogram from --ensemble-out")
    e.add_argument("csv")
    e.add_argument("--species", default="X")
    e.add_argument("--bins", type=int, default=40)
    common(e, "figures/ensemble.png")
    e.set_defaults(func=cmd_ensemble)

    st = sub.add_parser("stats", help="ensemble mean +/- sd from --stats-out")
    st.add_argument("csv")
    st.add_argument("--ode")
    common(st, "figures/stats.png")
    st.set_defaults(func=cmd_stats)

    sc = sub.add_parser("scale", help="small multiples across volume scalings")
    sc.add_argument("--pair", nargs=3, action="append", metavar=("OMEGA", "SSA", "ODE"),
                    required=True)
    sc.add_argument("--species", nargs="*", default=None)
    common(sc, "figures/scale.png")
    sc.set_defaults(func=cmd_scale)

    mj = sub.add_parser("majority", help="consensus probability vs initial split")
    mj.add_argument("--point", nargs=3, action="append", metavar=("X0", "Y0", "CSV"),
                    required=True)
    common(mj, "figures/majority.png")
    mj.set_defaults(func=cmd_majority)

    # `plot.py data/lv.csv` -- the bare form, with `trajectory` implied.
    argv = sys.argv[1:]
    if argv and argv[0] not in sub.choices and not argv[0].startswith("-"):
        argv.insert(0, "trajectory")

    args = p.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
