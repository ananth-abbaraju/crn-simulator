# crn-simulator

A chemical reaction network simulator in C++, with a Python analysis layer.

Chemical reaction networks are usually taught macroscopically, as continuous
deterministic ODEs. That description breaks down at the scale of a single cell,
where you are dealing with tens or hundreds of discrete molecules and the
dynamics are genuinely stochastic. This project implements both descriptions of
the same network and compares them.

## Contents

- **Gillespie SSA** (`src/`) — an exact stochastic simulation of the chemical
  master equation via the direct method. Samples the time to the next reaction
  and which reaction fires, rather than stepping on a fixed grid.
- **Deterministic integrator** — RK4 on the mass-action ODEs for the same
  network, for side-by-side comparison.
- **Analysis** (`python/`) — trajectory overlays, phase portraits, and ensemble
  statistics.

## Networks

- **Lotka–Volterra** (`A + X -> 2X`, `X + Y -> 2Y`, `Y -> B`) — predator–prey
  oscillation. The ODE orbits forever; the stochastic version random-walks in
  amplitude and eventually goes extinct.
- **Approximate majority** (`X + Y -> 2B`, `B + X -> 2X`, `B + Y -> 2Y`) — a CRN
  that computes the majority of two populations.

## Build

```
make
./crnsim --network lotka-volterra --out data/lv.csv
```

## Plot

```
python3 -m venv .venv
.venv/bin/pip install numpy matplotlib pandas
.venv/bin/python python/plot.py data/lv.csv
```

## Background

Built off material from ECE 381V (Unconventional Computation), UT Austin.
