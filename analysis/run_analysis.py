#!/usr/bin/env python3

import MDAnalysis as mda
from MDAnalysis.analysis import rms, distances
from MDAnalysis.analysis.align import AlignTraj
from MDAnalysis.analysis.rms import RMSF
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from sklearn.decomposition import PCA
from scipy.stats import gaussian_kde
import os

matplotlib.use("Agg")

os.makedirs("outputs/analysis", exist_ok=True)

print("Loading trajectory...")

u = mda.Universe("outputs/md/md.tpr", "outputs/md/md.xtc")

protein = u.select_atoms("protein")
ligand = u.select_atoms("resname LIG or resname UNL")

print("Protein atoms:", len(protein))
print("Ligand atoms detected:", len(ligand))

# ================= ALIGNMENT =================
print("Aligning trajectory...")
AlignTraj(u, u, select="protein and name CA", in_memory=True).run()

# ================= RMSD =================
print("Calculating RMSD...")
R = rms.RMSD(protein, protein, ref_frame=0).run()

time = R.results.rmsd[:, 1]  # ps
rmsd_values = R.results.rmsd[:, 2]

np.savetxt(
    "outputs/analysis/rmsd.dat",
    np.column_stack((time, rmsd_values)),
    header="time(ps) rmsd(A)"
)

plt.figure()
plt.plot(time, rmsd_values)
plt.title("Protein RMSD")
plt.xlabel("Time (ps)")
plt.ylabel("RMSD (Å)")
plt.savefig("outputs/analysis/rmsd.png")
plt.close()

# ================= RMSF =================
print("Calculating RMSF...")
rmsf = RMSF(protein).run()

np.savetxt(
    "outputs/analysis/rmsf.dat",
    rmsf.results.rmsf,
    header="RMSF(A)"
)

plt.figure()
plt.plot(rmsf.results.rmsf)
plt.title("Protein RMSF")
plt.xlabel("Residue index")
plt.ylabel("RMSF (Å)")
plt.savefig("outputs/analysis/rmsf.png")
plt.close()

# ================= Radius of gyration =================
print("Calculating Radius of Gyration...")
rg = []
for ts in u.trajectory:
    rg.append(protein.radius_of_gyration())

np.savetxt("outputs/analysis/rg.dat", rg, header="Rg(nm)")

plt.figure()
plt.plot(rg)
plt.title("Radius of Gyration")
plt.xlabel("Frame")
plt.ylabel("Rg (nm)")
plt.savefig("outputs/analysis/rg.png")
plt.close()

# ================= Free Energy Landscape =================
print("Calculating Free Energy Landscape...")

coords = []
ca = protein.select_atoms("name CA")

for ts in u.trajectory:
    coords.append(ca.positions.flatten())

coords = np.array(coords)

pca = PCA(n_components=2)
proj = pca.fit_transform(coords)

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
plt.savefig("outputs/analysis/fel.png")
plt.close()

# ================= Ligand distance =================
if len(ligand) > 0:
    print("Calculating protein–ligand distance...")

    distances_list = []

    for ts in u.trajectory:
        d = distances.distance_array(
            protein.center_of_mass().reshape(1, 3),
            ligand.center_of_mass().reshape(1, 3),
        )[0][0]

        distances_list.append(d)

    np.savetxt(
        "outputs/analysis/ligand_distance.dat",
        distances_list,
        header="distance(A)"
    )

    plt.figure()
    plt.plot(distances_list)
    plt.title("Protein–Ligand Distance")
    plt.xlabel("Frame")
    plt.ylabel("Distance (Å)")
    plt.savefig("outputs/analysis/ligand_distance.png")
    plt.close()
else:
    print("No ligand detected → skipping distance.")

print("Analysis complete.")

