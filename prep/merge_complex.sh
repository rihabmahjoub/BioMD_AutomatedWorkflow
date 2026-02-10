#!/bin/bash
set -e

echo "Merging protein and ligand (robust)..."

PROT=outputs/md/protein.gro
LIG=outputs/md/ligand.gro
OUT=outputs/md/complex_raw.gro

# Extract counts
PROT_ATOMS=$(sed -n '2p' $PROT)
LIG_ATOMS=$(sed -n '2p' $LIG)

TOTAL_ATOMS=$((PROT_ATOMS + LIG_ATOMS))

echo "Protein atoms: $PROT_ATOMS"
echo "Ligand atoms: $LIG_ATOMS"
echo "Total atoms: $TOTAL_ATOMS"

# Write header
head -n 1 $PROT > $OUT
echo $TOTAL_ATOMS >> $OUT

# Write protein atoms (skip header + atom count, skip box)
sed '1,2d;$d' $PROT >> $OUT

# Write ligand atoms (skip header + atom count, skip box)
sed '1,2d;$d' $LIG >> $OUT

# Write box from protein
tail -n 1 $PROT >> $OUT

echo "Complex created: $OUT"

