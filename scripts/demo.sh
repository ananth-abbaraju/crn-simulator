#!/usr/bin/env bash
# Regenerate every committed sample dataset and figure from scratch.
#
#   ./scripts/demo.sh
#
# Everything here is seeded, so the outputs are reproducible byte for byte.
set -euo pipefail

cd "$(dirname "$0")/.."

SIM=./crnsim
PY=${PY:-.venv/bin/python}
DATA=data/samples
FIG=figures

[ -x "$SIM" ] || { echo "build first: make" >&2; exit 1; }
[ -x "$PY" ]  || { echo "no venv: python3 -m venv .venv && .venv/bin/pip install numpy matplotlib pandas" >&2; exit 1; }

mkdir -p "$DATA" "$FIG"

echo "== Lotka-Volterra: single trajectory, stochastic vs deterministic =="
$SIM --network lotka-volterra --seed 1 \
     --out "$DATA/lv.csv" --ode-out "$DATA/lv_ode.csv"

echo "== Lotka-Volterra: the deterministic limit, three system volumes =="
for omega in 1 5 25; do
    $SIM --network lotka-volterra --seed 1 --scale "$omega" --sample-dt 0.01 \
         --out "$DATA/lv_omega$omega.csv" --ode-out "$DATA/lv_omega${omega}_ode.csv"
done

echo "== Lotka-Volterra: 1000-trajectory ensemble =="
# A long horizon, because prey extinction is a rare-ish event and the histogram
# needs its tail. Runs whose predators die out first blow up instead and are
# stopped by --max-count; they show up as stop=population-cap in the CSV.
$SIM --network lotka-volterra --seed 1000 --t-end 200 --ensemble 1000 \
     --ensemble-out "$DATA/lv_ensemble.csv"

# A shorter horizon for the mean/sd figure, before the amplitude random walk
# makes the ensemble mean meaningless.
$SIM --network lotka-volterra --seed 1000 --t-end 40 --ensemble 1000 \
     --stats-out "$DATA/lv_stats.csv" --ode-out "$DATA/lv_stats_ode.csv"

echo "== Approximate majority: consensus =="
$SIM --network approximate-majority --seed 3 --t-end 0.12 \
     --out "$DATA/am.csv" --ode-out "$DATA/am_ode.csv"

echo "== Approximate majority: does it actually compute majority? =="
# One ensemble per initial split, total population held at 190. The question is
# whether a lead of a handful of molecules is reliably amplified to consensus.
AM_POINTS=()
for x0 in 95 96 98 100 103 107 112 120; do
    y0=$((190 - x0))
    $SIM --network approximate-majority --seed $((7000 + x0)) --t-end 1.0 \
         --init "X=$x0" --init "Y=$y0" --ensemble 400 \
         --ensemble-out "$DATA/am_split$x0.csv" >/dev/null 2>&1
    AM_POINTS+=(--point "$x0" "$y0" "$DATA/am_split$x0.csv")
done

echo "== Brusselator, loaded from networks/brusselator.crn =="
$SIM --network networks/brusselator.crn --seed 1 --sample-dt 0.005 \
     --out "$DATA/brusselator.csv" --ode-out "$DATA/brusselator_ode.csv"

echo
echo "== Figures =="
$PY python/plot.py trajectory "$DATA/lv.csv" --ode "$DATA/lv_ode.csv" \
    -o "$FIG/lv_trajectory.png" \
    --title "Lotka-Volterra: 120 prey, 140 predators" \
    --subtitle "the same network under both descriptions; the ODE is periodic, the SSA is not"

$PY python/plot.py phase "$DATA/lv.csv" --ode "$DATA/lv_ode.csv" \
    -o "$FIG/lv_phase.png"

$PY python/plot.py scale \
    --pair 1  "$DATA/lv_omega1.csv"  "$DATA/lv_omega1_ode.csv" \
    --pair 5  "$DATA/lv_omega5.csv"  "$DATA/lv_omega5_ode.csv" \
    --pair 25 "$DATA/lv_omega25.csv" "$DATA/lv_omega25_ode.csv" \
    --species X -o "$FIG/lv_scale.png"

$PY python/plot.py ensemble "$DATA/lv_ensemble.csv" --species X \
    -o "$FIG/lv_extinction.png"

$PY python/plot.py stats "$DATA/lv_stats.csv" --ode "$DATA/lv_stats_ode.csv" \
    -o "$FIG/lv_ensemble_mean.png"

$PY python/plot.py trajectory "$DATA/am.csv" --ode "$DATA/am_ode.csv" \
    -o "$FIG/am_trajectory.png" \
    --title "Approximate majority: 100 X against 90 Y" \
    --subtitle "the SSA reaches the absorbing consensus at a definite time; the ODE only approaches it"

$PY python/plot.py majority "${AM_POINTS[@]}" -o "$FIG/am_majority.png"

$PY python/plot.py trajectory "$DATA/brusselator.csv" --ode "$DATA/brusselator_ode.csv" \
    -o "$FIG/brusselator_trajectory.png" \
    --title "Brusselator, loaded from a .crn file" \
    --subtitle "the ODE relaxes onto a limit cycle; the SSA follows it but drifts in phase"

echo
echo "done."
