#!/usr/bin/env python3
"""
md_quick_plot.py
================
Quick MD analysis for protein or protein–ligand simulations.

Inputs:
  - topology (.pdb or .gro)
  - trajectory (.xtc, .trr, or .dcd)

Outputs:
  - RMSD
  - RMSF
  - Radius of gyration
  - Free Energy Landscape (PCA-based)
  - Protein–ligand distance (if ligand present)

Author: EasyMD v2 (MD_quick_plot-inspired)
"""

import os
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde


class MDAnalyzer:
    def __init__(
        self,
        topology,
        trajectory,
        protein_selection="protein",
        ligand_selection=None,
        output_dir="outputs/analysis",
    ):
        self.topology = topology
        self.trajectory = trajectory
        self.protein_selection = protein_selection
        self.ligand_selection = ligand_selection
        self.output_dir = output_dir

        os.makedirs(self.output_dir, exist_ok=True)
        self.traj = md.load(self.trajectory, top=self.topology)

    # =========================
    # RMSD
    # =========================
    def plot_rmsd(self):
        protein_atoms = self.traj.topology.select(self.protein_selection)
        rmsd = md.rmsd(self.traj, self.traj, 0, atom_indices=protein_atoms)

        plt.figure()
        plt.plot(rmsd)
        plt.xlabel("Frame")
        plt.ylabel("RMSD (nm)")
        plt.title("Protein RMSD")
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/rmsd.png", dpi=300)
        plt.close()

    # =========================
    # RMSF
    # =========================
    def plot_rmsf(self):
        protein_atoms = self.traj.topology.select(self.protein_selection)
        rmsf = md.rmsf(self.traj, self.traj, 0, atom_indices=protein_atoms)

        plt.figure()
        plt.plot(rmsf)
        plt.xlabel("Atom index")
        plt.ylabel("RMSF (nm)")
        plt.title("Protein RMSF")
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/rmsf.png", dpi=300)
        plt.close()

    # =========================
    # Radius of Gyration
    # =========================
    def plot_rg(self):
        protein_atoms = self.traj.topology.select(self.protein_selection)
        rg = md.compute_rg(self.traj.atom_slice(protein_atoms))

        plt.figure()
        plt.plot(rg)
        plt.xlabel("Frame")
        plt.ylabel("Rg (nm)")
        plt.title("Radius of Gyration")
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/rg.png", dpi=300)
        plt.close()

    # =========================
    # Free Energy Landscape (PCA)
    # =========================
    def plot_fel(self):
        protein_atoms = self.traj.topology.select(self.protein_selection)
        xyz = self.traj.atom_slice(protein_atoms).xyz.reshape(self.traj.n_frames, -1)

        pca = PCA(n_components=2)
        proj = pca.fit_transform(xyz)

        kde = gaussian_kde(proj.T)
        xi, yi = np.mgrid[
            proj[:, 0].min():proj[:, 0].max():100j,
            proj[:, 1].min():proj[:, 1].max():100j,
        ]
        zi = kde(np.vstack([xi.flatten(), yi.flatten()]))

        plt.figure()
        plt.contourf(xi, yi, -np.log(zi.reshape(xi.shape)), levels=20)
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("Free Energy Landscape")
        plt.colorbar(label="Free energy (a.u.)")
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/fel.png", dpi=300)
        plt.close()

    # =========================
    # Protein–Ligand Distance
    # =========================
    def plot_protein_ligand_distance(self):
        if self.ligand_selection is None:
            return

        prot = self.traj.topology.select(self.protein_selection)
        lig = self.traj.topology.select(self.ligand_selection)

        distances = md.compute_contacts(
            self.traj, contacts=[(i, j) for i in prot for j in lig]
        )[0]

        min_dist = distances.min(axis=1)

        plt.figure()
        plt.plot(min_dist)
        plt.xlabel("Frame")
        plt.ylabel("Min distance (nm)")
        plt.title("Protein–Ligand Distance")
        plt.tight_layout()
        plt.savefig(f"{self.output_dir}/protein_ligand_distance.png", dpi=300)
        plt.close()

    # =========================
    # Run all
    # =========================
    def run_complete_analysis(self):
        self.plot_rmsd()
        self.plot_rmsf()
        self.plot_rg()
        self.plot_fel()
        self.plot_protein_ligand_distance()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--top", required=True)
    parser.add_argument("--traj", required=True)
    parser.add_argument("--ligand", default=None)
    args = parser.parse_args()

    analyzer = MDAnalyzer(
        topology=args.top,
        trajectory=args.traj,
        ligand_selection=args.ligand,
    )
    analyzer.run_complete_analysis()

