#!/usr/bin/env python3
"""Tests for crnsim.

Run with `make test`, or directly:

    .venv/bin/python tests/test_crnsim.py

Each test is a claim about the simulator that can be checked against something
known independently -- a conservation law, an absorbing state, a repeated seed --
rather than against a previous run of the simulator.
"""

import math
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SIM = os.path.join(ROOT, "crnsim")

FAILURES = []
TMP = None


def run(*args):
    proc = subprocess.run([SIM, *map(str, args)], capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(f"crnsim {' '.join(map(str, args))} failed:\n{proc.stderr}")
    return proc.stdout, proc.stderr


def net(text, name="t"):
    """Write a throwaway .crn file and return its path.

    Most tests below want a network small enough to have a known answer, which
    is exactly what the .crn parser is for.
    """
    path = os.path.join(TMP, f"{name}.crn")
    with open(path, "w") as f:
        f.write(text)
    return path


def out(name):
    return os.path.join(TMP, name)


def check(name, fn):
    try:
        fn()
    except AssertionError as e:
        FAILURES.append(name)
        print(f"FAIL  {name}\n      {e}")
    except Exception as e:  # noqa: BLE001
        FAILURES.append(name)
        print(f"ERROR {name}\n      {type(e).__name__}: {e}")
    else:
        print(f"ok    {name}")


def approx(actual, expected, tol, what):
    assert abs(actual - expected) <= tol, \
        f"{what}: got {actual:.6g}, expected {expected:.6g} +/- {tol:.3g}"


# --- propensity forms ------------------------------------------------------
#
# Each of these pins one branch of the propensity function by measuring the
# mean waiting time to the first (and only) reaction, which is 1/a exactly.

def first_event_mean(reactions, init, runs=4000, seed=11):
    """Mean waiting time for a network in which exactly one reaction can ever fire.

    Each network below is rigged so that firing once empties a reactant and
    leaves an absorbing state, which makes the recorded quiescent time the
    waiting time for a single exponential clock of rate a -- mean exactly 1/a.
    """
    path = net(f"tend 1000\ninit {init}\n{reactions}\n")
    run("--network", path, "--ensemble", runs, "--seed", seed,
        "--ensemble-out", out("e.csv"))
    df = pd.read_csv(out("e.csv"))
    assert df["quiescent_time"].notna().all(), "expected every run to absorb"
    assert (df["events"] == 1).all(), \
        f"the network must admit exactly one event, saw {sorted(df['events'].unique())}"
    return df["quiescent_time"].mean(), df["quiescent_time"].std() / math.sqrt(runs)


def test_unimolecular():
    # X -> 0 with k = 2 and one molecule: a = 2, mean waiting time 1/2.
    m, se = first_event_mean("X -> 0 ; k = 2", "X=1")
    approx(m, 0.5, 4 * se, "mean time to X -> 0")


def test_bimolecular_distinct():
    # X + Y -> 0 with n_X = 3, n_Y = 1: a = k * 3 * 1 = 3, mean 1/3.
    m, se = first_event_mean("X + Y -> 0 ; k = 1", "X=3 Y=1")
    approx(m, 1.0 / 3.0, 4 * se, "mean time to X + Y -> 0")


def test_bimolecular_same_species():
    # 2X + Y -> 0 with n_X = 4, n_Y = 1: a = k * C(4,2) = 6, mean 1/6.
    # Getting this branch wrong gives k * 4 * 3 = 12 or k * 4 * 4 = 16 instead,
    # both of which this test separates from 6 by many standard errors. The
    # spectator Y is there only to make the state absorbing after one event.
    m, se = first_event_mean("2X + Y -> 0 ; k = 1", "X=4 Y=1")
    approx(m, 1.0 / 6.0, 4 * se, "mean time to 2X + Y -> 0")


def test_trimolecular_same_species():
    # 3X + Y -> 0 with n_X = 5: a = k * C(5,3) = 10, mean 1/10.
    m, se = first_event_mean("3X + Y -> 0 ; k = 1", "X=5 Y=1")
    approx(m, 0.1, 4 * se, "mean time to 3X + Y -> 0")


def test_zeroth_order():
    # 0 -> X has a constant propensity k, independent of any count.
    path = net("tend 10\ninit X=0\n0 -> X ; k = 7\n")
    run("--network", path, "--ensemble", 2000, "--seed", 5,
        "--ensemble-out", out("z2.csv"))
    df = pd.read_csv(out("z2.csv"))
    # Poisson(7 * 10): mean 70, and the count equals the event count.
    approx(df["final_X"].mean(), 70.0, 4 * math.sqrt(70.0 / 2000), "mean X at t=10")
    approx(df["final_X"].var(), 70.0, 8.0, "var X at t=10")


# --- the SSA is exact ------------------------------------------------------

def test_ssa_matches_master_equation():
    """Pure decay has a closed-form solution; the SSA must reproduce it.

    X -> 0 with rate k is a binomial thinning: at time t each of the N0 initial
    molecules survives independently with probability exp(-kt), so the count is
    Binomial(N0, exp(-kt)) -- mean and variance both known exactly.
    """
    path = net("tend 1\ninit X=1000\nX -> 0 ; k = 1\n")
    run("--network", path, "--ensemble", 2000, "--seed", 99,
        "--stats-out", out("s.csv"), "--grid-dt", 0.1)
    df = pd.read_csv(out("s.csv"))
    p = np.exp(-df["t"].to_numpy())
    mean = 1000 * p
    sd = np.sqrt(1000 * p * (1 - p))
    # 2000 runs: the standard error on the mean is sd/sqrt(2000) <= 0.36.
    assert np.abs(df["mean_X"] - mean).max() < 1.5, \
        f"mean deviates by {np.abs(df['mean_X'] - mean).max():.3f}"
    assert np.abs(df["sd_X"] - sd).max() < 0.8, \
        f"sd deviates by {np.abs(df['sd_X'] - sd).max():.3f}"


def test_thread_count_does_not_change_results():
    """Each run owns its RNG stream, so the answer cannot depend on scheduling."""
    for t in (1, 4):
        run("--network", "lotka-volterra", "--t-end", 30, "--ensemble", 200,
            "--seed", 7, "--threads", t, "--ensemble-out", out(f"t{t}.csv"))
    a = pd.read_csv(out("t1.csv"))
    b = pd.read_csv(out("t4.csv"))
    pd.testing.assert_frame_equal(a, b)


# --- the deterministic integrator ------------------------------------------

def test_rk4_matches_closed_form():
    """Pure decay integrates to 1000 * exp(-t); RK4 must land on it."""
    path = net("tend 5\ninit X=1000\nX -> 0 ; k = 1\n")
    run("--network", path, "--ode-out", out("o.csv"),
        "--ode-step", 0.001, "--ode-sample-dt", 0.1)
    df = pd.read_csv(out("o.csv"))
    exact = 1000 * np.exp(-df["t"].to_numpy())
    rel = np.abs(df["X"] - exact) / np.maximum(exact, 1e-12)
    assert rel.max() < 1e-9, f"worst relative error {rel.max():.3g}"


def test_rk4_is_fourth_order():
    """Halving the step must cut the error by about 2^4."""
    errs = []
    for h in (0.2, 0.1):
        path = net("tend 2\ninit X=1000\nX -> 0 ; k = 1\n")
        run("--network", path, "--ode-out", out(f"o{h}.csv"),
            "--ode-step", h, "--ode-sample-dt", 2.0)
        df = pd.read_csv(out(f"o{h}.csv"))
        errs.append(abs(df["X"].iloc[-1] - 1000 * math.exp(-2.0)))
    ratio = errs[0] / max(errs[1], 1e-300)
    assert 10 < ratio < 24, f"error ratio {ratio:.1f}, expected ~16 for 4th order"


def test_ode_respects_the_conservation_law():
    """The integrator reads the same stoichiometry, so it conserves the same sum."""
    run("--network", "approximate-majority", "--ode-out", out("amo.csv"))
    df = pd.read_csv(out("amo.csv"))
    total = df["X"] + df["Y"] + df["B"]
    assert (total - 190).abs().max() < 1e-6, \
        f"total drifts by {(total - 190).abs().max():.3g}"


# --- invariants ------------------------------------------------------------

def test_conservation_law():
    """Approximate majority conserves X + Y + B on every single event.

    Every reaction in that network has two molecules on each side, so the total
    is an invariant of the stoichiometry -- and any indexing slip in the update
    loop would break it.
    """
    run("--network", "approximate-majority", "--seed", 4, "--out", out("am.csv"))
    df = pd.read_csv(out("am.csv"))
    total = df["X"] + df["Y"] + df["B"]
    assert (total == 190).all(), f"total ranges over {total.min()}..{total.max()}"


def test_counts_never_go_negative():
    """A reaction must never fire without enough reactants present."""
    for network in ("lotka-volterra", "approximate-majority"):
        run("--network", network, "--seed", 2, "--out", out("n.csv"))
        df = pd.read_csv(out("n.csv"))
        species = [c for c in df.columns if c != "t"]
        assert (df[species] >= 0).all().all(), f"{network} produced a negative count"


def test_absorbing_state_stops_the_run():
    """Approximate majority always reaches consensus, and consensus is absorbing.

    Both reactants of every reaction are gone in that state, so the total
    propensity is zero. Without the a0 == 0 guard the next waiting time is
    -log(u)/0 and the loop spins forever.
    """
    _, err = run("--network", "approximate-majority", "--seed", 4, "--out", out("a.csv"))
    assert "absorbing state reached at t =" in err, err

    df = pd.read_csv(out("a.csv"))
    last = df.iloc[-1]
    # Consensus: one species holds the whole population.
    assert last["B"] == 0 and (last["X"] == 0 or last["Y"] == 0), \
        f"expected a consensus state, got X={last['X']} Y={last['Y']} B={last['B']}"


def test_absorbing_at_t_zero():
    """A network with nothing to fire must stop at t = 0, not spin."""
    path = net("tend 100\ninit X=0\nX -> 0 ; k = 1\n")
    _, err = run("--network", path, "--out", out("z.csv"))
    assert "absorbing state reached at t = 0" in err, err


# --- reproducibility -------------------------------------------------------

def test_same_seed_same_trajectory():
    for name in ("r1.csv", "r2.csv"):
        run("--network", "lotka-volterra", "--seed", 42, "--out", out(name))
    with open(out("r1.csv")) as f1, open(out("r2.csv")) as f2:
        assert f1.read() == f2.read(), "identical seeds produced different output"


def test_different_seeds_differ():
    for seed, name in ((1, "s1.csv"), (2, "s2.csv")):
        run("--network", "lotka-volterra", "--seed", seed, "--out", out(name))
    with open(out("s1.csv")) as f1, open(out("s2.csv")) as f2:
        assert f1.read() != f2.read(), "different seeds produced identical output"


# --- the CLI ---------------------------------------------------------------

def test_unknown_network_is_an_error():
    proc = subprocess.run([SIM, "--network", "nope"], capture_output=True, text=True)
    assert proc.returncode != 0
    assert "lotka-volterra" in proc.stderr, "error should list the built-ins"


# --- the .crn parser -------------------------------------------------------

def test_parser_handles_comments_and_forms():
    path = net(
        "# a comment\n"
        "name custom     # trailing comment\n"
        "tend 3\n"
        "\n"
        "init A=10 B=5\n"
        "0      -> A      ; k = 1.5\n"
        "rxn 2A + B -> 3A ; k = 2e-3\n"
        "A -> 0 ; k=0.25\n",
        "parsed")
    _, err = run("--network", path, "--t-end", 0, "--out", out("c.csv"))
    assert "network: custom" in err, err
    assert "A=10" in err and "B=5" in err, err
    assert "2A + B -> 3A   k = 0.002" in err, err
    assert list(pd.read_csv(out("c.csv")).columns) == ["t", "A", "B"]


def test_parser_reports_the_failing_line():
    path = net("init X=1\nX -> ; k = 1\nY -> 0\n", "bad")
    proc = subprocess.run([SIM, "--network", path], capture_output=True, text=True)
    assert proc.returncode != 0, "expected a nonzero exit"
    assert ":3:" in proc.stderr, f"expected the failing line number: {proc.stderr}"


def test_parser_reads_the_brusselator():
    """The shipped example file parses and its trimolecular step fires."""
    path = os.path.join(ROOT, "networks", "brusselator.crn")
    _, err = run("--network", path, "--t-end", 40, "--sample-dt", 0.05,
                 "--seed", 1, "--out", out("br.csv"))
    assert "2X + Y -> 3X" in err, err

    df = pd.read_csv(out("br.csv"))
    # A limit cycle, not a fixed point: the oscillation has to be wide.
    assert df["X"].max() > 4 * df["X"].min(), \
        f"X ranged only over {df['X'].min()}..{df['X'].max()}"


def test_grid_sampling_matches_the_horizon():
    """--sample-dt lays down an exact uniform grid from 0 to t_end."""
    run("--network", "lotka-volterra", "--t-end", 10, "--sample-dt", 0.5,
        "--seed", 1, "--out", out("g.csv"))
    df = pd.read_csv(out("g.csv"))
    assert len(df) == 21, f"expected 21 grid points, got {len(df)}"
    assert abs(df["t"].iloc[-1] - 10.0) < 1e-9, df["t"].iloc[-1]
    steps = df["t"].diff().dropna()
    assert (steps - 0.5).abs().max() < 1e-9, "grid is not uniform"


def main():
    global TMP
    if not os.path.exists(SIM):
        sys.exit("crnsim not built -- run `make` first")

    with tempfile.TemporaryDirectory() as tmp:
        TMP = tmp
        for name, fn in sorted(globals().items()):
            if name.startswith("test_") and callable(fn):
                check(name[len("test_"):].replace("_", " "), fn)

    print()
    if FAILURES:
        print(f"{len(FAILURES)} failed: {', '.join(FAILURES)}")
        return 1
    print("all tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
