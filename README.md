# BioMD_AutomatedWorkflow
### Automated Molecular Dynamics Pipeline for Protein–Ligand Systems

![Linux](https://img.shields.io/badge/platform-Linux-blue)
![Python](https://img.shields.io/badge/python-3.10+-green)
![GROMACS](https://img.shields.io/badge/GROMACS-2021+-orange)
![License](https://img.shields.io/badge/license-MIT-lightgrey)



BioMD_AutomatedWorkflow is a fully automated molecular dynamics pipeline designed to prepare, simulate, and analyze protein–ligand systems using GROMACS and MDAnalysis.

#The workflow performs:

- Protein preparation
- Ligand preparation
- Complex assembly
- Solvation and ion addition
- Energy minimization
- NVT and NPT equilibration
- Production molecular dynamics
- Automated analysis and reporting

This pipeline is designed for:
- Reproducible research
- Thesis work
- Rapid MD screening
- Educational use

---


# Workflow Overview

Protein PDB + Ligand SDF
│
▼
Preparation (Protein + Ligand)
│
▼
Complex assembly
│
▼
Box + Solvation + Ions
│
▼
Energy Minimization
│
▼
NVT Equilibration
│
▼
NPT Equilibration
│
▼
Production MD
│
▼
Automated Analysis
│
▼
Figures + PDF Report


---

# Project Structure

BioMD_AutomatedWorkflow/
│
├── analysis/ # Analysis and report generation
├── prep/ # Preparation scripts
├── mdp/ # Simulation parameters
├── example_run/ # Reproducible example
├── inputs/ # User input folder
├── run_easyMD.sh # Main pipeline
└── environment.yml # Conda environment


---

# Requirements

- Linux
- Conda / Miniconda
- GROMACS 2021+
- Python 3.10+

---

# Installation

Clone repository:

git clone https://github.com/rihabmahjoub/BioMD_AutomatedWorkflow.git

cd BioMD_AutomatedWorkflow


Create environment:

conda env create -f environment.yml
conda activate easyMD_v2


---

# Quick Start

1. Place your protein and ligand in:

inputs/


Example:

inputs/protein.pdb
inputs/ligand.sdf


2. Run pipeline:

bash run_easyMD.sh


3. Results will appear in:

outputs/


---

# Output Files

The pipeline automatically generates:

- RMSD plot
- RMSF plot
- Radius of gyration
- Free energy landscape
- Protein–ligand distance
- PDF simulation report

All results are saved in:

outputs/analysis/


---

# Example Run

---

```markdown
## Reproducibility

All steps of the workflow are fully automated and reproducible.
An example dataset is provided in:

example_run/

This ensures reproducibility of results and enables rapid testing.

A reproducible example is available:

example_run/


Contains:
- Input structures
- Output figures
- Generated report

---

# Analysis Pipeline

The analysis module performs:

- Trajectory alignment
- RMSD calculation
- RMSF calculation
- Radius of gyration
- Free energy landscape (PCA)
- Ligand distance monitoring
- Automatic PDF report generation

---

# Citation

If you use this workflow in academic work, please cite:

Mahjoub, R. (2026)
BioMD_AutomatedWorkflow: Automated Molecular Dynamics Pipeline
GitHub Repository


---

# Author

Rihab Mahjoub  
PhD Researcher – Computational Biology  

GitHub:
https://github.com/rihabmahjoub

---

# License

MIT License
