#!/usr/bin/env python3

import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

if len(sys.argv) != 3:
    print("Usage: python clean_protein.py input.pdb output.pdb")
    sys.exit(1)

input_pdb = sys.argv[1]
output_pdb = sys.argv[2]

print("Cleaning protein structure with PDBFixer...")

fixer = PDBFixer(filename=input_pdb)

# Remove heterogens (keeps protein, removes ligands/ions/waters)
fixer.removeHeterogens(keepWater=False)

# Fix missing atoms and residues
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Add hydrogens at physiological pH
fixer.addMissingHydrogens(pH=7.0)

# Write cleaned protein
with open(output_pdb, "w") as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print(f"Cleaned protein saved to: {output_pdb}")

