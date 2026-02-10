# Example Run

This folder contains a minimal example demonstrating the BioMD Automated Workflow.

A reproducible example is provided in `example_run/`.

To test:

cp example_run/inputs/* inputs/
bash run_easyMD.sh

## Inputs
- protein.pdb : ubiquitin structure
- ligand.sdf : test ligand (glycine)

## Outputs
Example analysis results:
- RMSD
- RMSF
- Free Energy Landscape
- Radius of gyration

## How to reproduce

Copy inputs to main folder:

cp example_run/inputs/* inputs/

Run:

bash run_easyMD.sh
