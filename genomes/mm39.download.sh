#!/bin/bash
# Download and prepare mm39 (GRCm39) genome for RNAmotifs
set -euo pipefail

GENOME="mm39"
DIR="$(cd "$(dirname "$0")" && pwd)/${GENOME}"
mkdir -p "$DIR"
cd "$DIR"

echo "Downloading ${GENOME} 2bit file..."
wget -q "https://hgdownload.cse.ucsc.edu/goldenPath/${GENOME}/bigZips/${GENOME}.2bit"

echo "Converting to FASTA and splitting by chromosome..."
TWOBITTOFA="$(dirname "$0")/twoBitToFa"
if [ ! -x "$TWOBITTOFA" ]; then
    echo "Error: twoBitToFa not found at $TWOBITTOFA"
    exit 1
fi

for CHR in $(seq 1 19) X Y; do
    echo "  chr${CHR}..."
    "$TWOBITTOFA" -seq="chr${CHR}" "${GENOME}.2bit" "chr${CHR}.fa"
    grep -v "^>" "chr${CHR}.fa" | tr -d '\n' > "chr${CHR}.string"
    rm "chr${CHR}.fa"
done

# Clean up 2bit
rm -f "${GENOME}.2bit"

echo "Done. Genome ready in ${DIR}"
