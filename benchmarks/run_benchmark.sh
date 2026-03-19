#!/bin/bash
# Performance benchmark: m3_light (Python) vs rnamotifs_search (C++)
#
# Usage: bash examples/run_benchmark.sh [/path/to/theRNAmars]
#
# Run from the RNAmotifs2 root directory. Requires the mm9 genome
# .string files in genomes/mm9/.
set -e

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MARS="${1:-$(dirname "$ROOT")/theRNAmars}"
INPUT="$ROOT/examples/NOVA.txt"
BINARY="$ROOT/build/rnamotifs_search"
GENOME="mm9"
HW=15
HMIN=4
PTH=0.5

echo "============================================"
echo "RNAmotifs2 Performance Benchmark"
echo "============================================"
echo "Input: $INPUT ($(wc -l < "$INPUT") lines)"
echo "Genome: $GENOME"
echo "Parameters: hw=$HW, min_height=$HMIN, pth=$PTH"
echo "C++ binary: $BINARY"
echo "m3_light: $MARS/m3_light/"
echo ""

# Build C++ if needed
if [ ! -f "$BINARY" ]; then
    echo "Building C++ rnamotifs_search..."
    mkdir -p "$ROOT/build"
    ( cd "$ROOT/build" && cmake .. && make -j"$(nproc)" rnamotifs_search )
    echo ""
fi

run_cpp() {
    local cores=$1
    local tag="BENCH_CPP${cores}"
    echo "============================================"
    echo "C++ rnamotifs_search ($cores core$([ "$cores" -gt 1 ] && echo s))"
    echo "============================================"
    rm -rf "$ROOT/results/$tag" "$ROOT/regions/$tag"
    { time "$BINARY" "$INPUT" "$ROOT" "$GENOME" "$tag" "$ROOT" \
        $HW $HMIN $PTH "$cores" 2>/dev/null ; } 2>&1
    rm -rf "$ROOT/results/$tag" "$ROOT/regions/$tag"
    echo ""
}

# C++ benchmarks
for c in 1 4 8 12; do
    run_cpp "$c"
done

# Python m3_light benchmark
echo "============================================"
echo "Python m3_light (1 core)"
echo "============================================"
if [ -d "$MARS/m3_light" ]; then
    export PYTHONPATH="$MARS"
    rm -rf "$MARS/m3_light/results/BENCH_PY" "$MARS/m3_light/regions/BENCH_PY"
    mkdir -p "$MARS/m3_light/regions/BENCH_PY"
    cp "$INPUT" "$MARS/m3_light/regions/BENCH_PY/BENCH_PY.tab"
    ( cd "$MARS" && { time python3 -c "
import m3_light
m3_light.find_tetramers('BENCH_PY', '$GENOME', '$HW', '$HMIN', '$PTH')
" ; } 2>&1 )
    rm -rf "$MARS/m3_light/results/BENCH_PY" "$MARS/m3_light/regions/BENCH_PY"
else
    echo "  Skipped: $MARS/m3_light not found."
    echo "  Pass the path to theRNAmars as argument: bash examples/run_benchmark.sh /path/to/theRNAmars"
fi

echo ""
echo "============================================"
echo "Benchmark complete."
echo "============================================"
