#!/bin/bash

# Distributed under the MIT License.
# See LICENSE.txt for details.

# Idle + contended sweep of SphericalTransformBenchmark for the SPHEREPACK
# replacement study. Run on an EXCLUSIVELY allocated compute node.
#
# Usage: tools/SpherepackClusterSweep.sh <build_dir> [ncopies]
#   build_dir  spectre build directory containing bin/SphericalTransformBenchmark
#   ncopies    concurrent copies for the contended runs (default: all cores
#              reported by nproc; on SMT machines pass the number of PHYSICAL
#              cores and check the core numbering with lscpu -e first).
#
# Writes: spherepack-sweep-results-<host>.txt  (idle sweeps + contended medians)
#         spherepack-sweep-contended-<host>/   (raw per-copy outputs)
set -euo pipefail

BUILD_DIR=${1:?usage: $0 <build_dir> [ncopies]}
NCOPIES=${2:-$(nproc)}
BIN="$BUILD_DIR/bin/SphericalTransformBenchmark"
[[ -x "$BIN" ]] || {
  echo "ERROR: $BIN not found or not executable." \
       "Build it with: ninja -C $BUILD_DIR SphericalTransformBenchmark" >&2
  exit 1
}
HOST=$(hostname -s)
OUT="spherepack-sweep-results-$HOST.txt"
RAW="spherepack-sweep-contended-$HOST"
mkdir -p "$RAW"

# The benchmark must be single-threaded per copy.
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

{
  echo "host: $HOST"
  date
  if [[ -r /proc/cpuinfo ]]; then
    grep -m1 "model name" /proc/cpuinfo
  fi
  command -v lscpu > /dev/null &&
    lscpu | grep -E "^CPU\(s\)|Thread|Core|Socket|L1d|L2|L3"
  echo "binary: $BIN"
  echo "contended copies: $NCOPIES"
  echo
} > "$OUT"

run_idle() {
  local title=$1
  shift
  echo "=== IDLE: $title ===" >> "$OUT"
  "$@" >> "$OUT" 2>&1
}

echo "--- idle sweeps ---"
# Evolution regime: 50 GH-like components, 8 radial points.
for l in 8 10 12 16 20 24 32; do
  if [[ $l -le 24 ]]; then iters=150; else iters=60; fi
  for mode in 0 1 2; do
    run_idle "evolution l=$l mode=$mode" "$BIN" "$l" 8 50 "$iters" 10 "$mode"
  done
done
# Horizon/Shape regime: single field, single radius, up to l=100.
for l in 8 12 16 20 24 32 48 64 80 100; do
  if [[ $l -le 32 ]]; then
    iters=300
  elif [[ $l -le 64 ]]; then
    iters=80
  else
    iters=30
  fi
  for mode in 0 1 2 3; do
    run_idle "horizon l=$l mode=$mode" "$BIN" "$l" 1 1 "$iters" 5 "$mode"
  done
done

# Contended: NCOPIES identical pinned copies; the per-copy min-of-iters is the
# contended time. Representative sizes: smooth + prime n_phi in each regime.
median_summary() {
  # Median per op label across the per-copy outputs given as arguments.
  awk -F: '
    / us\/callset| us$| bytes$/ {
      label = $1
      gsub(/^[ \t]+/, "", label)
      split($2, a, " ")
      vals[label] = vals[label] " " a[1]
      count[label]++
    }
    END {
      for (label in vals) {
        n = split(vals[label], v, " ")
        asort_n = n
        # insertion sort (tiny n)
        for (i = 2; i <= n; i++) {
          x = v[i] + 0
          j = i - 1
          while (j >= 1 && v[j] + 0 > x) { v[j + 1] = v[j]; j-- }
          v[j + 1] = x
        }
        mid = int((n + 1) / 2)
        printf "  %s: median %s over %d copies\n", label, v[mid], n
      }
    }' "$@"
}

echo "--- contended sweeps ($NCOPIES copies) ---"
echo "=== CONTENDED ($NCOPIES pinned copies, per-copy min, medians) ===" \
  >> "$OUT"
# l n_r n_comp mode iters
for spec in "12 8 50 0 60" "12 8 50 2 60" "20 8 50 0 40" "20 8 50 2 40" \
            "48 1 1 0 120" "48 1 1 2 120" "100 1 1 0 25" "100 1 1 1 25" \
            "100 1 1 2 25"; do
  read -r l nr nc mode iters <<< "$spec"
  tag="l${l}_nr${nr}_nc${nc}_m${mode}"
  echo "  running $tag ..."
  pids=()
  for ((i = 0; i < NCOPIES; i++)); do
    if command -v taskset > /dev/null; then
      taskset -c "$i" "$BIN" "$l" "$nr" "$nc" "$iters" 5 "$mode" \
        > "$RAW/${tag}_copy$i.txt" 2>&1 &
    else
      "$BIN" "$l" "$nr" "$nc" "$iters" 5 "$mode" \
        > "$RAW/${tag}_copy$i.txt" 2>&1 &
    fi
    pids+=($!)
  done
  # Wait on explicit PIDs (bare `wait` hangs on bash >= 5.1 with process
  # substitutions; explicit PIDs are safe everywhere).
  wait "${pids[@]}"
  {
    echo "--- $tag ---"
    median_summary "$RAW/${tag}"_copy*.txt
  } >> "$OUT"
done

echo
echo "Done. Results in $OUT (raw contended outputs in $RAW/)."
