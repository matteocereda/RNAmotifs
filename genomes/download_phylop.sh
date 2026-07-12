#!/bin/bash
# Download PhyloP conservation scores and convert to per-chromosome binary arrays
# Usage: ./download_phylop.sh <genome>
#
# Supports: hg19, hg38 (BigWig), mm9, mm10 (wigFix per chromosome)
set -euo pipefail

GENOME="${1:-}"
if [ -z "$GENOME" ]; then
    echo "Usage: $0 <hg19|hg38|mm9|mm10>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DIR="${SCRIPT_DIR}/${GENOME}"
BWTOBG="${SCRIPT_DIR}/bigWigToBedGraph"

if [ ! -d "$DIR" ]; then
    echo "Error: genome directory $DIR not found. Download the genome first."
    exit 1
fi

# Get chromosome list from .string files
CHROMS=$(ls "$DIR"/chr*.string 2>/dev/null | sed 's/.*\///' | sed 's/\.string//' | sort -V)
if [ -z "$CHROMS" ]; then
    echo "Error: no .string files in $DIR"
    exit 1
fi

# --- Converter: wigFix.gz -> binary float array ---
convert_wigfix() {
    local INFILE="$1"
    local CHR_LEN="$2"
    local OUTFILE="$3"
    python3 -c "
import sys, gzip, array

chrom_len = int(sys.argv[1])
scores = array.array('f', [0.0] * chrom_len)

pos = 0
step = 1
with gzip.open(sys.argv[2], 'rt') as f:
    for line in f:
        line = line.strip()
        if line.startswith('fixedStep'):
            parts = dict(p.split('=') for p in line.split() if '=' in p)
            pos = int(parts['start']) - 1  # 0-based
            step = int(parts.get('step', '1'))
        elif line and not line.startswith('#'):
            if pos < chrom_len:
                scores[pos] = float(line)
            pos += step

with open(sys.argv[3], 'wb') as f:
    scores.tofile(f)
" "$CHR_LEN" "$INFILE" "$OUTFILE"
}

# --- Converter: bedGraph (from bigWig) -> binary float array ---
convert_bedgraph() {
    local BW_FILE="$1"
    local CHR="$2"
    local CHR_LEN="$3"
    local OUTFILE="$4"
    "$BWTOBG" -chrom="$CHR" "$BW_FILE" /dev/stdout 2>/dev/null | \
    python3 -c "
import sys, array

chrom_len = int(sys.argv[1])
scores = array.array('f', [0.0] * chrom_len)
for line in sys.stdin:
    parts = line.strip().split('\t')
    start, end, val = int(parts[1]), int(parts[2]), float(parts[3])
    for i in range(start, min(end, chrom_len)):
        scores[i] = val
with open(sys.argv[2], 'wb') as f:
    scores.tofile(f)
" "$CHR_LEN" "$OUTFILE"
}

# --- Genome-specific download logic ---

case "$GENOME" in
    hg19)
        URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw"
        FORMAT="bigwig"
        ;;
    hg38)
        URL="https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/hg38.phyloP100way.bw"
        FORMAT="bigwig"
        ;;
    mm9)
        URL_BASE="https://hgdownload.cse.ucsc.edu/goldenpath/mm9/phyloP30way/vertebrate"
        FORMAT="wigfix"
        ;;
    mm10)
        URL_BASE="https://hgdownload.cse.ucsc.edu/goldenpath/mm10/phyloP60way/mm10.60way.phyloP60way.bw"
        FORMAT="bigwig"
        ;;
    *)
        echo "Error: unknown genome '$GENOME'"
        exit 1
        ;;
esac

echo "PhyloP for $GENOME (format: $FORMAT)"
echo "Output: $DIR/chr*.phylop.bin"
echo ""

if [ "$FORMAT" = "bigwig" ]; then
    # BigWig: download once, extract per chromosome
    if [ ! -x "$BWTOBG" ]; then
        echo "Error: bigWigToBedGraph not found at $BWTOBG"
        echo "Download from: https://hgdownload.cse.ucsc.edu/admin/exe/"
        exit 1
    fi

    BW_FILE="$DIR/phyloP.bw"
    if [ ! -f "$BW_FILE" ]; then
        echo "Downloading BigWig..."
        wget -q --show-progress -O "$BW_FILE" "$URL"
    fi

    for CHR in $CHROMS; do
        OUTFILE="$DIR/${CHR}.phylop.bin"
        [ -f "$OUTFILE" ] && echo "  $CHR: exists, skip" && continue
        CHR_LEN=$(wc -c < "$DIR/${CHR}.string")
        echo -n "  $CHR ($CHR_LEN bp)..."
        convert_bedgraph "$BW_FILE" "$CHR" "$CHR_LEN" "$OUTFILE"
        echo " done"
    done

elif [ "$FORMAT" = "wigfix" ]; then
    # wigFix.gz: download per chromosome
    for CHR in $CHROMS; do
        OUTFILE="$DIR/${CHR}.phylop.bin"
        [ -f "$OUTFILE" ] && echo "  $CHR: exists, skip" && continue

        WIG_FILE="$DIR/${CHR}.phyloP.wigFix.gz"
        WIG_URL="${URL_BASE}/${CHR}.phyloP30way.wigFix.gz"

        echo -n "  $CHR: downloading..."
        if ! wget -q -O "$WIG_FILE" "$WIG_URL" 2>/dev/null; then
            echo " not available, skipping"
            rm -f "$WIG_FILE"
            continue
        fi

        CHR_LEN=$(wc -c < "$DIR/${CHR}.string")
        echo -n " converting ($CHR_LEN bp)..."
        convert_wigfix "$WIG_FILE" "$CHR_LEN" "$OUTFILE"
        rm -f "$WIG_FILE"
        echo " done"
    done
fi

echo ""
echo "PhyloP scores ready in $DIR/"
ls -lh "$DIR"/*.phylop.bin 2>/dev/null | head -5
echo "..."
ls "$DIR"/*.phylop.bin 2>/dev/null | wc -l
echo "total .phylop.bin files"
