#!/usr/bin/env bash
# Thread scaling of the ensemble runner.
#
# Every run owns its RNG stream, so the total event count is identical at every
# thread count -- which is both the correctness check and what makes the
# events/second numbers comparable.
set -euo pipefail
cd "$(dirname "$0")/.."

RUNS=${RUNS:-1000}
TEND=${TEND:-40}

echo "$RUNS Lotka-Volterra trajectories to t = $TEND"
echo
for t in 1 2 4 8; do
    ./crnsim --network lotka-volterra --seed 1000 --t-end "$TEND" \
             --ensemble "$RUNS" --threads "$t" --ensemble-out /dev/null 2>&1 \
        | grep -E 'threads in|events/s' | sed "s/^/  /"
done
