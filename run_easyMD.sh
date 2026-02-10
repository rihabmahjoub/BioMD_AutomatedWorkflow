#!/bin/bash
set -e

echo "========== EasyMD v2 (Stable) =========="

# Activate conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate easyMD_v2

mkdir -p outputs/md outputs/analysis

################################
# STEP 1: Prepare protein
################################
echo "1️⃣ Preparing protein..."
bash prep/prepare_protein.sh

START_STRUCTURE="outputs/md/protein.gro"

################################
# STEP 2: Ligand detection
################################
if [ -f inputs/ligand.sdf ]; then
    echo "Ligand detected → preparing ligand"
    python prep/prepare_ligand.py
    bash prep/merge_complex.sh
    START_STRUCTURE="outputs/md/complex_raw.gro"
else
    echo "No ligand detected → protein-only MD"
fi

################################
# STEP 3: Box (larger box = more stable)
################################
echo "2️⃣ Creating box..."
gmx editconf \
-f $START_STRUCTURE \
-o outputs/md/boxed.gro \
-c -d 1.8 -bt cubic

################################
# STEP 4: Solvate
################################
echo "3️⃣ Solvating..."
gmx solvate \
-cp outputs/md/boxed.gro \
-cs spc216.gro \
-o outputs/md/solv.gro \
-p outputs/md/topol.top

################################
# STEP 5: Add ions
################################
echo "4️⃣ Adding ions..."

gmx grompp \
-f mdp/minim.mdp \
-c outputs/md/solv.gro \
-p outputs/md/topol.top \
-o outputs/md/ions.tpr \
-maxwarn 1

echo "SOL" | gmx genion \
-s outputs/md/ions.tpr \
-o outputs/md/solv_ions.gro \
-p outputs/md/topol.top \
-neutral

################################
# STEP 6: Energy Minimization 1
################################
echo "5️⃣ Energy Minimization (1)..."

gmx grompp \
-f mdp/minim.mdp \
-c outputs/md/solv_ions.gro \
-p outputs/md/topol.top \
-o outputs/md/em.tpr \
-maxwarn 1

gmx mdrun -deffnm outputs/md/em

################################
# STEP 6b: Energy Minimization 2 (stabilization)
################################
echo "Energy Minimization (2)..."

gmx grompp \
-f mdp/minim.mdp \
-c outputs/md/em.gro \
-p outputs/md/topol.top \
-o outputs/md/em2.tpr \
-maxwarn 1

gmx mdrun -deffnm outputs/md/em2

################################
# STEP 7: Create index
################################
echo "Creating index..."

gmx make_ndx -f outputs/md/em2.gro -o outputs/md/index.ndx <<EOF
q
EOF

################################
# STEP 8: NVT
################################
echo "6️⃣ NVT equilibration..."

gmx grompp \
-f mdp/nvt.mdp \
-c outputs/md/em2.gro \
-r outputs/md/em2.gro \
-p outputs/md/topol.top \
-n outputs/md/index.ndx \
-o outputs/md/nvt.tpr \
-maxwarn 1

gmx mdrun -deffnm outputs/md/nvt

################################
# STEP 9: NPT
################################
echo "7️⃣ NPT equilibration..."

gmx grompp \
-f mdp/npt.mdp \
-c outputs/md/nvt.gro \
-r outputs/md/nvt.gro \
-p outputs/md/topol.top \
-n outputs/md/index.ndx \
-o outputs/md/npt.tpr \
-maxwarn 1

gmx mdrun -deffnm outputs/md/npt

################################
# STEP 10: Production MD
################################
echo "8️⃣ Production MD..."

gmx grompp \
-f mdp/md.mdp \
-c outputs/md/npt.gro \
-p outputs/md/topol.top \
-n outputs/md/index.ndx \
-o outputs/md/md.tpr \
-maxwarn 1

gmx mdrun -deffnm outputs/md/md

################################
# STEP 11: Analysis
################################
echo "9️⃣ Analysis..."
python analysis/run_analysis.py
python analysis/generate_report.py

echo ""
echo "========== EasyMD FINISHED =========="
echo "Results in outputs/analysis/"

