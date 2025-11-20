# Principal Component Analysis (PCA), Free Energy Landscapes (FEL), and Extreme Conformations

## 1. PCA Workflow (MDAnalysis)

We performed PCA on the protein–ligand trajectory using **MDAnalysis**, replicating the standard GROMACS workflow.

- **Selection:** `backbone or resname PYC`
- **Alignment:** RMSD-align all frames to backbone
- **PCA:** Covariance matrix → eigenvalues/eigenvectors → projections

```python
import MDAnalysis as mda
from MDAnalysis.analysis import align, pca

u = mda.Universe("start.pdb", "prod.xtc")
sel = u.select_atoms("backbone or resname PYC")

align.AlignTraj(u, u, select="backbone", in_memory=True).run()
pc = pca.PCA(u, select="backbone or resname PYC", align=False)
pc.run()
proj = pc.transform(sel)
```

---

## 2. Free Energy Landscapes (FEL)

We computed FELs by projecting onto PC pairs and applying the Boltzmann relation:

```python
import numpy as np, matplotlib.pyplot as plt
kBT = 2.494  # kJ/mol at 300K

def free_energy_surface(pcX, pcY, bins=100):
    H, xedges, yedges = np.histogram2d(pcX, pcY, bins=bins, density=True)
    F = -kBT * np.log(H + 1e-10)
    X, Y = np.meshgrid(0.5*(xedges[:-1]+xedges[1:]),
                       0.5*(yedges[:-1]+yedges[1:]))
    return X, Y, F.T

X, Y, F12 = free_energy_surface(proj[:,0], proj[:,1])
plt.contourf(X, Y, F12, levels=50, cmap="viridis")
plt.colorbar(label="Free Energy (kJ/mol)")
plt.xlabel("PC1"); plt.ylabel("PC2")
plt.title("FEL: PC1 vs PC2")
plt.show()
```

- Repeat for PC2–PC3 and PC3–PC1.
- FEL plots reveal conformational basins and transitions.

---

## 3. Extreme Conformations Along PC1–PC3

We generated extreme conformations displaced along ±PC directions, analogous to `gmx anaeig -extr`.

```python
mean_coords = pc.mean.flatten()
eigvecs = pc.eigenvectors

def extreme_conformation(pc_index, displacement=5.0):
    vec = eigvecs[pc_index]
    plus_coords = mean_coords + displacement * vec
    minus_coords = mean_coords - displacement * vec
    return plus_coords.reshape(-1, 3), minus_coords.reshape(-1, 3)

for i in range(3):
    plus, minus = extreme_conformation(i, displacement=5.0)
    sel.positions = plus; sel.write(f"extreme_PC{i+1}_plus.pdb")
    sel.positions = minus; sel.write(f"extreme_PC{i+1}_minus.pdb")
```

- Six `.pdb` files generated: `extreme_PC1_plus/minus.pdb`, etc.
- Load into PyMOL to visualize conformational extremes.

---

## 4. Morph Trajectories Along PCs

To visualize **continuous motion** along PCs, we interpolate between mean and extreme conformations:

```python
def morph_trajectory(pc_index, displacement=5.0, nframes=20):
    vec = eigvecs[pc_index]
    frames = []
    for alpha in np.linspace(-displacement, displacement, nframes):
        coords = mean_coords + alpha * vec
        sel.positions = coords.reshape(-1, 3)
        fname = f"morph_PC{pc_index+1}_{alpha:.2f}.pdb"
        sel.write(fname)
        frames.append(fname)
    return frames

# Example: morph along PC1
morph_files_PC1 = morph_trajectory(0, displacement=5.0, nframes=20)
```

- Produces a series of `.pdb` snapshots interpolating from −extreme → mean → +extreme.
- Load sequentially in PyMOL/VMD to animate the motion along PC1–PC3.

---

## 5. Observations

- **PC1:** Large-scale backbone motions dominate.
- **PC2:** Ligand conformational rearrangements visible.
- **PC3:** Tail “whipping” motion consistent with FEL valley transitions.
- **Morph trajectories:** Provide smooth visualization of conformational pathways, complementing static FEL plots.

---

## 6. Summary

- PCA via MDAnalysis reproduces GROMACS results when selections, alignment, and weighting are matched.
- FEL surfaces reveal conformational basins and transitions.
- Extreme conformations and morph trajectories along PC1–PC3 provide structural insight into motions driving FEL transitions.
- Workflow is fully automated in Python/Jupyter, eliminating manual `.xvg → .xpm → .dat` conversions.

## TODO:

- Add code to interpolate a polynomial for the FEL and sample the interpolant at a chosen cross-section.

