#!/bin/bash
set -e

echo "Preparing protein..."

cd outputs/md

gmx pdb2gmx \
-f ../../inputs/protein_clean.pdb \
-o protein.gro \
-p topol.top \
-water tip3p \
-ff amber99sb-ildn \
-ignh

cd ../..
echo "Protein preparation completed."


