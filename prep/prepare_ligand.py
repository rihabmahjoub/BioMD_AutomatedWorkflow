import subprocess
import shutil

print("Preparing ligand...")

# Convert SDF → MOL2 with OpenBabel
subprocess.run([
    "obabel",
    "inputs/ligand.sdf",
    "-O", "outputs/md/ligand.mol2",
    "-h"
], check=True)

print("Running ACPYPE...")
subprocess.run([
    "acpype",
    "-i", "outputs/md/ligand.mol2",
    "-b", "LIG",
    "-o", "gmx"
], check=True)

print("Fixing residue name to LIG...")

# Copy gro and rename residue UNL → LIG
with open("LIG.acpype/LIG_GMX.gro") as f:
    lines = f.readlines()

for i in range(2, len(lines)-1):
    lines[i] = lines[i].replace("UNL", "LIG")

with open("outputs/md/ligand.gro", "w") as f:
    f.writelines(lines)

shutil.copy("LIG.acpype/LIG_GMX.itp", "outputs/md/ligand.itp")

print("Ligand ready.")

topol = "outputs/md/topol.top"
lig_itp = "outputs/md/ligand.itp"

# 1. Ensure ligand.itp is included
with open(topol, "r") as f:
    content = f.read()

if '#include "ligand.itp"' not in content:
    content = content.replace(
        '#include "amber99sb-ildn.ff/forcefield.itp"',
        '#include "amber99sb-ildn.ff/forcefield.itp"\n#include "ligand.itp"'
    )

# 2. Ensure ligand appears in molecules section
if "LIG 1" not in content:
    content = content.replace(
        "[ molecules ]",
        "[ molecules ]\nLIG 1"
    )

with open(topol, "w") as f:
    f.write(content)

print("Topology updated with ligand.")

