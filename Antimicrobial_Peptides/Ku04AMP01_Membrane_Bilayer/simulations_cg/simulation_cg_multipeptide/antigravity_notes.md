Based on your system setup (MARTINI 3 coarse-grained / CHARMM36 multi-peptide membrane bilayer system with N copies of the KU04 peptide, md.xtc / traj_nojump.xtc, and ku04_gim_multi.top / .tpr / .gro), here is a
  complete guide and clean Python code to compute the individual z-COM (Center of Mass) coordinates for each peptide molecule.
  ──────
  ### Key Technical Considerations for Multi-Peptide Bilayers

  1. Splitting into Individual Peptides: Calling .center_of_mass() on all protein atoms at once returns a single overall center of mass. To track each peptide independently, you must split the selection into individual
  AtomGroup objects (via .fragments or .split("molecule")).
  2. Periodic Boundary Conditions (PBC Unwrapping): If a peptide spans across the periodic box boundaries (z = 0 or
    z = L
         z
    
  ), computing standard COM without unwrapping will distort the COM calculation. MDAnalysis handles this natively with unwrap=True (or by using your pre-processed traj_nojump.xtc).
  3. Relative z-COM vs Absolute z: Bilayers undergo undulations and box fluctuations in NPT ensembles. It is usually best practice to compute relative height:
    Δzᵢ(t) = z            (t) - z                   (t)
              peptideᵢ,COM       membrane,\text{COM}
    
  or relative to the phosphate headgroups.
  ──────
  ### Option 1: MDAnalysis (Recommended for Flexibility & PBC Handling)
  Below is a complete, production-ready script using MDAnalysis:

    import numpy as np
    import pandas as pd
    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    
    # 1. Load Universe (Use topology TPR/GRO/PDB + XTC trajectory)
    # Standard practice: use traj_nojump.xtc if available, or pass unwrap=True
    top_file = "ku04_gim_multi.tpr"   # or "ku04_gim_multi.gro" / "ku04_gim_multi.top"
    traj_file = "traj_nojump.xtc"     # or "md.xtc"
    
    u = mda.Universe(top_file, traj_file)
    
    # 2. Select protein and membrane groups
    protein = u.select_atoms("protein")  # or "resname KU04"
    membrane = u.select_atoms("resname DVPE POPE POPG DOPE TOCL") # Adjust lipid resnames if AA
    
    # Split protein into individual peptide molecules
    # .fragments automatically groups covalently bonded/molecule-connected beads
    peptides = protein.fragments  # List of AtomGroups, one per peptide molecule
    num_peptides = len(peptides)
    print(f"Found {num_peptides} individual peptide molecules.")
    
    # 3. Iterate over trajectory and record z-COM coordinates
    n_frames = len(u.trajectory)
    times = np.zeros(n_frames)
    z_com_abs = np.zeros((n_frames, num_peptides)) # Absolute z-COM
    z_com_rel = np.zeros((n_frames, num_peptides)) # Relative to membrane z-COM
    
    for ts_idx, ts in enumerate(u.trajectory):
        times[ts_idx] = u.trajectory.time  # ps (or ns)
        
        # Calculate membrane center of mass z-coordinate
        memb_z = membrane.center_of_mass(unwrap=True)[2]
        
        # Calculate each peptide's z-COM
        for pep_idx, peptide in enumerate(peptides):
            # unwrap=True prevents PBC split artifacts
            pep_com = peptide.center_of_mass(unwrap=True)
            z_abs = pep_com[2]
            
            z_com_abs[ts_idx, pep_idx] = z_abs
            z_com_rel[ts_idx, pep_idx] = z_abs - memb_z
    
    # 4. Convert to Pandas DataFrame for easy saving and plotting
    cols = [f"Peptide_{i+1}" for i in range(num_peptides)]
    df_rel = pd.DataFrame(z_com_rel, columns=cols)
    df_rel.insert(0, "Time_ps", times)
    
    # Save results
    df_rel.to_csv("peptide_z_com_relative.csv", index=False)
    print("Saved z-COM coordinates to 'peptide_z_com_relative.csv'")
    
    # 5. Quick Visualization
    plt.figure(figsize=(10, 6))
    for i in range(num_peptides):
        plt.plot(df_rel["Time_ps"] / 1000.0, df_rel[f"Peptide_{i+1}"], label=f"Pep {i+1}", alpha=0.8)
    
    plt.axhline(0, color="black", linestyle="--", linewidth=1.5, label="Membrane Midplane")
    plt.xlabel("Time (ns)")
    plt.ylabel(r"Relative $z$-COM ($\text{\AA}$ or nm)")
    plt.title("Individual Peptide Center of Mass ($z$) vs. Time")
    plt.legend(ncol=3, loc="upper right")
    plt.tight_layout()
    plt.savefig("peptide_z_com_trajectories.png", dpi=300)
    plt.show()
    ──────
  ### Exploring Alternatives to MDAnalysis

  Depending on your performance needs, workflow integration, or post-processing requirements, here are three strong alternatives:

   Tool                   | Best For                              | Speed / Performance                | Pros                                                    | Cons
  ------------------------|---------------------------------------|------------------------------------|---------------------------------------------------------|----------------------------------------------------------
   MDAnalysis             | General analysis, flexible selections | Moderate (Python loop)             | Native unwrap=True, robust selection syntax, handles    | Slower on very long trajectories (≥10⁵ frames)
                          |                                       |                                    | MARTINI & AA seamlessly                                 |
   MDTraj                 | High-speed batch processing           | High (C-accelerated vectorization) | Fast frame loading (md.iterload), low memory footprint  | Requires manual mass-weighting if topology lacks atomic
                          |                                       |                                    |                                                         | weights
   GROMACS CLI (gmx traj) | Native C++ execution                  | Very High (Compiled C++)           | No Python overhead, natively handles TPR topology       | Generates .xvg text files; less flexible for inline
                          |                                       |                                    | masses and PBC (-pbc mol)                               | Python pipelines
   LiPyphilic             | Bilayer-specific analysis             | High (Built on MDAnalysis)         | Built specifically for membrane systems (lipid          | Additional dependency
                          |                                       |                                    | headgroup planes, leaflet assignment, insertion depth)  |
  ──────
  ### Alternative 1: MDTraj (High Performance / Chunked Processing)

  If trajectory reading in MDAnalysis feels slow, MDTraj loads XTC files significantly faster using chunked iteration:

    import mdtraj as md
    import numpy as np
    
    # Load topology (PDB or GRO) and stream trajectory in chunks
    top_file = "ku04_gim_multi.gro"
    chunk_size = 1000
    
    # Identify atom indices for each peptide
    top = md.load_topology(top_file)
    # Assuming 12 peptide chains/molecules or selecting by residue range
    peptide_indices = [top.select(f"resname KU04 and moleculetype {i}") for i in range(12)] 
    # Or slice by atom indices per molecule
    
    z_com_list = []
    for chunk in md.iterload("md.xtc", top=top_file, chunk=chunk_size):
        # Image molecules to fix PBC wrapping
        chunk.image_molecules(inplace=True)
        
        # Compute COM for each frame in chunk
        chunk_z = np.zeros((chunk.n_frames, len(peptide_indices)))
        for idx, atom_indices in enumerate(peptide_indices):
            # Compute center of mass (nm)
            pep_traj = chunk.atom_slice(atom_indices)
            # md.compute_center_of_mass computes COM across all frames in chunk
            chunk_z[:, idx] = md.compute_center_of_mass(pep_traj)[:, 2]
        
        z_com_list.append(chunk_z)
    
    z_com_all = np.vstack(z_com_list) # Array of shape (total_frames, num_peptides)
    ──────
  ### Alternative 2: GROMACS Native (gmx traj / gmx select)

  If you want maximum calculation speed without Python memory overhead, use GROMACS' built-in tools to extract individual COMs directly:

    # 1. Create an index file containing separate groups for each peptide molecule
    # (If index.ndx has [ KU04_1 ], [ KU04_2 ], ..., [ KU04_12 ])
    
    # 2. Run gmx traj with PBC molecule treatment
    gmx traj -f traj_nojump.xtc -s md.tpr -n index.ndx -com -ox z_com.xvg -pbc mol
    ──────
  ### Alternative 3: LiPyphilic (Membrane-Specific Python Package)

  LiPyphilic https://lipyphilic.oxacml.com/ is a specialized Python library built on top of MDAnalysis specifically designed for membrane-protein systems. It includes dedicated modules to compute insertion depth relative
  to upper/lower lipid leaflets rather than just a single z plane, which accounts for local membrane curvature.
  ──────
  ### Recommendation Summary

  • Use MDAnalysis (Option 1 above) if you want an easy, self-contained Python script with built-in PBC unwrapping.
  • Ensure you use relative coordinates (

    z        - z
     peptide    bilayer

  ) to eliminate rigid-body translation/drift of the bilayer along the Z-axis.