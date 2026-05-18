# Simulation of Polymixin-B Cyclic protein with outer/inner Cell Membrane of Gram-Negative bacteria
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
1. All topology information of the Gram Negative Bacterial Membranes can be downloaded from the [Charmm-Gui website](https://charmm-gui.org/?doc=archive&lib=biomembrane):
   1. **G-OM (outer membranes of Gram-negative bacteria)**
        Download PDB file of last snapshot: [last.pdb](https://charmm-gui.org/archive/biomembrane/G-OM/last.pdb)  
        Download whole simulation input files: [G-OM.tar.gz](https://charmm-gui.org/archive/biomembrane/G-OM/G-OM.tar.gz)

        Size (Å³): 81.44 x 81.44 x 129.69  
        \# water / Na+ / Cl- / Ca2+: 14934 / 70 / 40 / 175  
        Temperature (K): 310.15  
        Lipid composition:

        | Lipid Name | Lipid Head/Tail | #Lipids in Leaflets | | APL (Å) | | Charge(e) |
        |---|---|---|---|---|---|---|
        | | | Outer | Inner | Outer | Inner | |
        | ECLIPA | E.coli R1 core-Lipid A* | 35 | 0 | 183 | N/A | -10 |
        | PPPE | PE (16:0/16:1 (9Z)) | 0 | 75 | N/A | 60.54 | 0 |
        | PVPG | PG (16:0/18:1 (11Z)) | 0 | 20 | N/A | 62.1 | -1 |
        | PVCL2 | CL (1'-[16:0/18:1 (11Z)],3'-[16:0/18:1 (11Z)]) | 0 | 5 | N/A | 124.47 | -2 |
        | **TOTAL** | | **35** | **100** | | | |
    2. - **G-IM (inner membranes of Gram-negative bacteria)**
        Download PDB file of last snapshot: [last.pdb](https://charmm-gui.org/archive/biomembrane/G-IM/last.pdb)  
        Download whole simulation input files: [G-IM.tar.gz](https://charmm-gui.org/archive/biomembrane/G-IM/G-IM.tar.gz)

            Size (Å³): 79.79 x 79.79 x 100  
            \# water / Na+ / Cl-: 11550 / 76 / 30  
            Temperature (K): 310.15  
            Lipid composition: ([show chemical structures](#))

            | Lipid Name | Lipid Head/Tail | #Lipids in Leaflets | APL (Å) | Charge(e) |
            |---|---|---|---|---|
            | | | Outer or Inner | Inner or Inner | |
            | PMPE | PE (16:0/cy17:0) | 46 | 62.25 | 0 |
            | POPE | PE (16:0/18:1 (9Z)) | 13 | 62.18 | 0 |
            | QMPE | PE (15:0/cy17:0) | 12 | 62.19 | 0 |
            | OYPE | PE (18:1 (9Z)/16:1 (9Z)) | 8 | 62.35 | 0 |
            | PMPG | PG (16:0/cy17:0) | 10 | 64.22 | -1 |
            | PYPG | PE (16:0/16:1 (9Z)) | 9 | 64.04 | -1 |
            | PVCL2 | CL (1'-[16:0/18:1 (11Z)],3'-[16:0/18:1 (11Z)]) | 2 | 127.43 | -2 |
            | **TOTAL** | | **100** | | |
2. Structure of Pomymyxin-B from [PubChem id: 49800004](https://pubchem.ncbi.nlm.nih.gov/compound/49800004)
   1. The [pubchem sdf file](pmbn_pubchem.sdf) does not contain residue ids, but the residue sequence is available on their website:![Polymyxin-B residue sequence](pmbn_seq.svg)
   2. The residue-labelled PDB generated from the PubChem SDF is [polymyxin_b_pubchem_residues.pdb](polymyxin_b_pubchem_residues.pdb). It uses the residue sequence `MOA-DAB-THR-DAB-DAB-DDB-DPH-LEU-DAB-DAB-THR`, where `MOA` is the 6-methyloctanoyl acyl group, `DDB` is D-Dab, and `DPH` is D-Phe. The lactam ring closure is preserved in the PDB `CONECT` records between `DAB 5 ND` and `THR 11 C`.

      Keep the following fragment map before sending the structure to CHARMM-GUI/CGenFF. If topology generation treats polymyxin B as one ligand residue, this table preserves the atom groups needed later for per-fragment membrane interaction-energy analysis. The `PDB atom serials` column applies to `polymyxin_b_pubchem_residues.pdb` and `polymyxin_b_3d.pdb`; the `SDF atom IDs` column applies to `pmbn_pubchem.sdf` and `polymyxin_b_3d.sdf`.

      | Fragment group | Chemical fragment | PDB atom serials | SDF atom IDs | Atom names |
      |---|---|---:|---|---|
      | `PMB_MOA1` | 6-methyloctanoyl tail | 1-27 | 73,12,58,43,39,34,33,44,71,57,90-92,97,98,103-106,127-131,146-148 | C,O,C2,C3,C4,C5,C6,C7,C8,CM,H01-H17 |
      | `PMB_DAB2` | Dab 2 | 28-42 | 24,61,75,13,74,79,29,134,151,152,154,162,163,182,183 | N,CA,C,O,CB,CG,ND,H01-H08 |
      | `PMB_THR3` | Thr 3 | 43-56 | 23,59,69,11,70,8,80,132,145,153,164-166,174 | N,CA,C,O,CB,OG1,CG2,H01-H07 |
      | `PMB_DAB4` | Dab 4 | 57-71 | 20,49,60,7,62,76,27,115,135,136,139,155,156,177,178 | N,CA,C,O,CB,CG,ND,H01-H08 |
      | `PMB_DAB5` | Dab 5, lactam side-chain N | 72-85 | 18,37,52,3,41,54,22,95,100,101,122,123,125,142 | N,CA,C,O,CB,CG,ND,H01-H07 |
      | `PMB_DAB6` | Dab 6, verify D/L assignment | 86-100 | 17,40,53,4,51,72,26,99,120,121,124,149-150,171,172 | N,CA,C,O,CB,CG,ND,H01-H08 |
      | `PMB_DPH7` | D-Phe 7 | 101-120 | 16,35,42,2,45,67,81-85,93,107,108,119,167,168,175,176,181 | N,CA,C,O,CB,CG,CD1,CD2,CE1,CE2,CZ,H01-H09 |
      | `PMB_LEU8` | Leu 8 | 121-139 | 14,30,38,1,31,32,46,47,86-89,96,109-114 | N,CA,C,O,CB,CG,CD1,CD2,H01-H11 |
      | `PMB_DAB9` | Dab 9 | 140-154 | 15,36,55,5,48,66,25,94,102,116,117,143,144,169,170 | N,CA,C,O,CB,CG,ND,H01-H08 |
      | `PMB_DAB10` | Dab 10 | 155-169 | 19,50,64,9,63,77,28,118,133,137,138,157,158,179,180 | N,CA,C,O,CB,CG,ND,H01-H08 |
      | `PMB_THR11` | Thr 11, lactam carbonyl | 170-183 | 21,56,68,10,65,6,78,126,140,141,159-161,173 | N,CA,C,O,CB,OG1,CG2,H01-H07 |

      For GROMACS interaction-energy reruns, make one index group per fragment using the final topology atom numbers. The successful charged CHARMM-GUI Ligand Reader run in `PMB_topol/charmm-gui-7883126705` uses the SDF atom order for atoms 1-183 and appends the five added ammonium hydrogens as atoms 184-188. Also verify residue 6 before production: canonical polymyxin B usually has `D-Phe` as the D-amino-acid residue, while the current template labels residue 6 as `DDB`/D-Dab.

      Final PMB fragment groups for postprocessing with `PMB_topol/charmm-gui-7883126705/gromacs/PMB.itp`:

      | Fragment group | Final `PMB.itp` atom numbers | Notes |
      |---|---|---|
      | `PMB_MOA1` | 12,33,34,39,43,44,57,58,71,73,90-92,97,98,103-106,127-131,146-148 | 6-methyloctanoyl tail |
      | `PMB_DAB2` | 13,24,29,61,74,75,79,134,151,152,154,162,163,182,183,188 | Free Dab side-chain ammonium; added H is 188 |
      | `PMB_THR3` | 8,11,23,59,69,70,80,132,145,153,164-166,174 | Thr 3 |
      | `PMB_DAB4` | 7,20,27,49,60,62,76,115,135,136,139,155,156,177,178,186 | Free Dab side-chain ammonium; added H is 186 |
      | `PMB_DAB5` | 3,18,22,37,41,52,54,95,100,101,122,123,125,142 | Lactam side-chain N; no added ammonium H |
      | `PMB_DAB6` | 4,17,26,40,51,53,72,99,120,121,124,149,150,171,172,185 | Free Dab side-chain ammonium; added H is 185 |
      | `PMB_DPH7` | 2,16,35,42,45,67,81-85,93,107,108,119,167,168,175,176,181 | D-Phe 7 |
      | `PMB_LEU8` | 1,14,30-32,38,46,47,86-89,96,109-114 | Leu 8 |
      | `PMB_DAB9` | 5,15,25,36,48,55,66,94,102,116,117,143,144,169,170,184 | Free Dab side-chain ammonium; added H is 184 |
      | `PMB_DAB10` | 9,19,28,50,63,64,77,118,133,137,138,157,158,179,180,187 | Free Dab side-chain ammonium; added H is 187 |
      | `PMB_THR11` | 6,10,21,56,65,68,78,126,140,141,159-161,173 | Thr 11; lactam carbonyl |

      When PMB is inserted into a solvated membrane system, add the PMB molecule's atom-number offset to these final `PMB.itp` atom numbers. For example, if PMB begins at atom `50001` in the system, then `PMB_DAB9` is system atoms `50005,50015,50025,...,50184`. Save the exact system-level groups in the index file used for analysis.

      Interaction-energy analysis after MD:

      1. Build fragment and membrane index groups from the final run input or final coordinates. Use the fragment map above, but convert it to the final solvated-system atom numbers after topology generation. For the membrane, make groups that match the question being asked, for example `LPS`, `LipidA`, `Phosphate`, `InnerLeaflet`, `OuterLeaflet`, `PPPE`, `PVPG`, or `PVCL2`.

          ```bash
          gmx make_ndx -f md.tpr -o interaction_groups.ndx
          ```

          Name the polymyxin groups consistently with the table, for example `PMB_DAB2`, `PMB_THR3`, `PMB_DPH7`, and so on. Save a copy of the exact final atom-number ranges used, because these are the analysis definitions.

      2. Prepare a rerun `.mdp` that writes pair energies between the selected groups. Keep the physical settings compatible with the production system, but set `energygrps` to the fragment and membrane groups:

          ```ini
          energygrps = PMB_MOA1 PMB_DAB2 PMB_THR3 PMB_DAB4 PMB_DAB5 PMB_DAB6 PMB_DPH7 PMB_LEU8 PMB_DAB9 PMB_DAB10 PMB_THR11 LPS PPPE PVPG PVCL2
          ```

          For a focused calculation, use fewer groups, such as one PMB fragment and one membrane group. Too many `energygrps` can make the rerun slow and creates a large `.edr` file.

      3. Recompute energies on the saved trajectory. This does not change the trajectory; it only recalculates energies for the frames in `md.xtc`.

          ```bash
          gmx grompp -f rerun_energy.mdp -c md.gro -p topol.top -n interaction_groups.ndx -t md.cpt -o rerun_energy.tpr
          gmx mdrun -s rerun_energy.tpr -rerun md.xtc -deffnm rerun_energy
          gmx energy -f rerun_energy.edr -o pmb_membrane_interactions.xvg
          ```

          Extract terms such as `Coul-SR:PMB_DAB2-LPS` and `LJ-SR:PMB_DAB2-LPS`. The usual descriptive interaction energy is:

          ```text
          E_interaction = Coul-SR(fragment, membrane) + LJ-SR(fragment, membrane)
          ```

      4. Interpret the result as a nonbonded energy decomposition, not a binding free energy. The rerun reports force-field short-range Coulomb and Lennard-Jones terms between the chosen groups under the simulation settings. It does not include entropy, conformational reorganization, solvent free energy, long-range electrostatic decomposition in a simple residue-pair sense, or membrane deformation cost. It is useful for ranking which PMB fragments contact the membrane strongly, especially when paired with contacts, distances, and H-bond occupancies.

      5. Analyze hydrogen bonds separately. H-bonds are not an energy term in the CHARMM/GROMACS nonbonded decomposition; they emerge from electrostatics and geometry. Use MDAnalysis or `gmx hbond` to calculate PMB-fragment to membrane H-bond counts, occupancies, donor/acceptor identities, and lifetimes. The fragment map above can be reused to build the MDAnalysis atom selections.

      6. Treat MM/PBSA or MM/GBSA as a separate follow-up analysis. It can estimate a more binding-like free-energy decomposition, but it has stronger assumptions for membranes, charged ligands, ions, and lipid/water boundaries. Read the method details before using it for final claims, and use the GROMACS rerun energies plus H-bond/contact analysis as the first-pass diagnostic.
   3. A 3D conformer can be generated with the RDKit script [../scripts/generate_3d_conformer.py](../scripts/generate_3d_conformer.py):

        ```bash
        conda run -n htvs python ../scripts/generate_3d_conformer.py \
          --input pmbn_pubchem.sdf \
          --output-prefix polymyxin_b_3d \
          --template-pdb polymyxin_b_pubchem_residues.pdb \
          --keep-existing-hydrogens \
          --num-conformers 200 \
          --num-threads 0
        ```

        The script writes `polymyxin_b_3d.pdb` and `polymyxin_b_3d.sdf`. The `--template-pdb` option keeps the residue labels from `polymyxin_b_pubchem_residues.pdb`; `--keep-existing-hydrogens` keeps the PubChem SDF atom order aligned with that template. On an HPC node, increase `--num-conformers` for a broader conformer search and set `--num-threads` to the scheduler-allocated CPU count, or leave it as `0` to let RDKit use all visible cores.
   4. Protonate PMB and generate the GROMACS topology with CHARMM-GUI Ligand Reader & Modeler using the charged PMB state.

        The neutral PubChem SDF gives total charge `0`, which is not appropriate for the usual PMB membrane-binding simulation because the Dab side chains should be cationic and interact strongly with LPS/lipid A phosphate groups. Before uploading to CHARMM-GUI, generate a charged `+5` SDF with [protonate_pmb_dab.py](protonate_pmb_dab.py):

        ```bash
        conda run -n hpc python ./protonate_pmb_dab.py \
          --input polymyxin_b_3d.sdf \
          --output polymyxin_b_3d_plus5.sdf
        ```

        The script preserves the original 183-atom order and appends one new hydrogen to each of the five free Dab side-chain amines. Upload `polymyxin_b_3d_plus5.sdf` to Ligand Reader & Modeler, not PDB Reader.

        Protonated these free side-chain nitrogens:

        | Fragment | Atom in template | SDF atom ID | Expected state |
        |---|---|---:|---|
        | `PMB_DAB2` | `DAB 2 ND` | 29 | `-NH3+` |
        | `PMB_DAB4` | `DAB 4 ND` | 27 | `-NH3+` |
        | `PMB_DAB6` | `DAB/DDB 6 ND` | 26 | `-NH3+` |
        | `PMB_DAB9` | `DAB 9 ND` | 25 | `-NH3+` |
        | `PMB_DAB10` | `DAB 10 ND` | 28 | `-NH3+` |

        Did not protonate `DAB 5 ND` (`SDF atom ID 22`) as a free side-chain amine: it forms the lactam ring closure with `THR 11 C`. The backbone/amide nitrogens should remain neutral amides.

        > **Note:** The charged SDF has 188 atoms rather than 183 because five new side-chain amine hydrogens are appended at the end. The original SDF atom IDs in the fragment map remain valid for the original heavy atoms and hydrogens. After CHARMM-GUI topology generation, rebuild the final analysis atom map from the generated `PMB.itp`/coordinates because CHARMM-GUI may rename atoms.

        After downloading the CHARMM-GUI output, check the generated files before using them:

        ```bash
        awk '/^\[ atoms \]/{a=1;next} /^\[ bonds \]/{a=0} a && $1 ~ /^[0-9]+$/ {q+=$7} END {printf "PMB total charge = %.6f\n", q}' gromacs/PMB.itp
        grep -n '^RESI' pmb/pmb.rtf
        ```

        The accepted topology should report approximately `+5.000000` in `PMB.itp`, and the CHARMM residue line should be `RESI pmb 5.000` or equivalent. If it reports `0.000`, discard that topology and rerun Ligand Reader with the side-chain amines protonated.

        The successful charged topology is in `PMB_topol/charmm-gui-7883126705`. It contains `188` PMB atoms and total charge `+5.000000` in `gromacs/PMB.itp`. CHARMM-GUI completed normally with no unknown coordinates, `ligandrm.str` reports `SET NCHARGE = 5`, and `pmb/pmb.rtf` reports `RESI pmb 5.000 ! param penalty=4.600 ; charge penalty=4.252`. CGenFF penalties below 10 indicate fair analogies, so these values are acceptable for this first model. CHARMM-GUI renamed atoms generically (`O1`, `N12`, `C56`, etc.); use the final fragment table above when building future GROMACS index groups.
        
        **Update:** The vacuum simulation ran okay. 
        * [rungmx.sh](PMB_topol/charmm-gui-7883126705/gromacs/rungmx.sh)
        * [video](PMB_topol/charmm-gui-7883126705/gromacs/gromacs_sim.mp4)
3. Insert PMB into the equilibrated G-OM outer membrane and run GROMACS MD.

   ### G-OM archive stages

   The downloaded [G-OM.tar.gz](G-OM.tar.gz) contains CHARMM-format files only (`.inp`, `.psf`, `.crd`, `.dcd`). There are no GROMACS input files in the archive. The stages are:

   | Stage | Key files | What it represents |
   |---|---|---|
   | step3 | `step3_packing.pdb`, `step3_glipid.inp` | LPS/Lipid-A packing and glycolipid ring building |
   | step4 | `step4_lipid.psf`, `step4.2_waterbox.*`, `step4.3_ion.*` | Lipid bilayer assembly, solvation, and ion placement (including Ca²⁺ for LPS) |
   | step5 | `step5_assembly.pdb`, `step5_assembly.psf`, `step5_assembly.crd` | Fully assembled system (lipids + water + ions), before any MD; energy-minimized with `step5_input.inp` |
   | step6.1–6.6 | `step6.x_equilibration.inp` | Six-stage NVT/NPT equilibration in CHARMM with progressively relaxed restraints (see table below) |
   | step7 | `step7_production.inp`, `step7_1000.dcd` | Production MD, 250,000 steps at 310.15 K |
   | last snapshot | [G-OM-last.pdb](G-OM-last.pdb) | Last frame of step7 — the fully equilibrated membrane ready for use |

   The six CHARMM equilibration stages use the following restraint schedule. `wforce`/`tforce`/`mforce` are harmonic force constants (kcal/mol/Å²) for water exclusion from the hydrophobic core, lipid tail order, and lipid head-group position, respectively. `fcis`/`fc2`/`cring` are dihedral restraint constants (kcal/mol/rad²) for cis double bonds, C2 chirality, and sugar ring conformation (critical for the LPS glycan rings):

   | Stage | Timestep | Steps | Duration | wforce / tforce / mforce | fcis / fc2 / cring |
   |---|---|---|---|---|---|
   | EM (step6.1) | — | 8 × SD+ABNR rounds (250–500 steps each) | — | 2.5 / 2.5 / 2.5 | 250 / 250 / 250 |
   | step6.1 | 1 fs | 125,000 | 125 ps | 2.5 / 2.5 / 2.5 | 250 / 250 / 250 |
   | step6.2 | 1 fs | 125,000 | 125 ps | 2.5 / 2.5 / 2.5 | 100 / 100 / 100 |
   | step6.3 | 1 fs | 125,000 | 125 ps | 1.0 / 1.0 / 1.0 | 50 / 50 / 50 |
   | step6.4 | 1 fs | 250,000 | 250 ps | 0.5 / 0.5 / 0.5 | 50 / 50 / 50 |
   | step6.5 | 1 fs | 250,000 | 250 ps | 0.1 / 0.1 / 0.1 | 25 / 25 / 25 |
   | step6.6 | 1 fs | 250,000 | 250 ps | 0 (none) | 0 (none) |
   | step7 | 1 fs | 250,000 | 250 ps | 0 (none) | 0 (none) |

   **PMB insertion point:** use `G-OM-last.pdb` (the step7 final snapshot). The membrane is already equilibrated; there is no need to re-run CHARMM steps 3–7. The relevant CHARMM force field parameter files are in `tarball/toppar/`, in particular `toppar_all36_lipid_lps.str` and `toppar_all36_lipid_bacterial.str` for ECLIPA.

   ### Option A: CHARMM-GUI "Add Molecules to Membrane" (recommended)

   This is the preferred route because CHARMM-GUI correctly handles the LPS/Lipid-A glycan topology (ECLIPA) during CHARMM→GROMACS conversion, which is very error-prone to do manually.

   > **Note on "Add Molecules to Bilayer System":** This Bilayer Builder sub-option was attempted first but does not work for the G-OM system. Uploading `GOM_last.pdb` causes CHARMM-GUI to run CGenFF on all residues in the PDB. The LPS glycan residues (AGAL, AGLC, AHEP, AKDO, ECLI) are not in CGenFF's database, producing the error *"CGenFF failed to parameterize force field for your molecule"*. There is no PSF upload field in the interface to suppress this. The solution is to rebuild the G-OM membrane fresh using **Build New Bilayer** and request GROMACS output directly, then add PMB manually to the resulting GROMACS files.

   #### Step A1 — Rebuild G-OM in CHARMM-GUI Bilayer Builder with GROMACS output

   Go to [CHARMM-GUI Membrane Builder → Bilayer Builder](https://charmm-gui.org/?doc=input/membrane.bilayer). On the landing page, select **"Build New Bilayer"**. Configure as follows:

   **Membrane composition** — match the original G-OM system exactly:

   | Leaflet | Lipid | Count |
   |---|---|---|
   | Outer | ECLIPA (E. coli R1 core-Lipid A) | 35 |
   | Inner | PPPE (PE 16:0/16:1 9Z) | 75 |
   | Inner | PVPG (PG 16:0/18:1 11Z) | 20 |
   | Inner | PVCL2 (CL 16:0/18:1 11Z) | 5 |

   **System settings:**
   - Box size: 81.44 × 81.44 Å lateral (X/Y); Z will be determined automatically (~130 Å)
   - Water model: TIP3P
   - Temperature: 310.15 K
   - Ion concentration: 0.15 M NaCl + neutralising Ca²⁺ for LPS (CHARMM-GUI will suggest amounts based on LPS charge)

   **Force field / output page:**
   - Force field: **CHARMM36**
   - Output format: **GROMACS**
   
   Download the resulting package. The key GROMACS files will be in `gromacs/`:

   | File | Contents |
   |---|---|
   | `gromacs/step5_input.pdb` | Assembled system coordinates in GROMACS-compatible PDB |
   | `gromacs/topol.top` | Master topology (ECLIPA, PPPE, PVPG, PVCL2, water, ions) |
   | `gromacs/step6.x_equilibration.mdp` | Six staged equilibration `.mdp` files |
   | `gromacs/step7_production.mdp` | Production `.mdp` |

   #### Step A2 — Manually insert PMB into the GROMACS system
   First, copy the contents of the `gromacs/` directory mentioned above to a new 
   `simulation` directory. It has a [README](simulation/README).

   PMB is added by merging its coordinates into the membrane GRO file and updating the topology. Use the pre-built CHARMM-GUI Ligand Reader output from `PMB_topol/charmm-gui-7883126705` — do not re-run CGenFF or re-upload the SDF, as that would produce a different atom ordering and potentially different partial charges.

   **2a. Use `step5_input.gro` as the membrane GRO:**

   CHARMM-GUI writes both `step5_input.pdb` and `step5_input.gro`. The PDB has no `CRYST1` record, so converting it with `gmx editconf` produces a GRO with zero box vectors, which causes a fatal `Empty diagonal for a 3-dimensional periodic box` error in `grompp`. Use the GRO file directly instead:

   ```bash
   cp step5_input.gro membrane.gro
   ```

   Verify the box line (last line) is non-zero:
   ```bash
   tail -1 membrane.gro
   # Expected e.g.:   8.14401   8.14401  19.47255   0.00000 ...
   ```

   **2b. Place PMB above the outer leaflet:**

   Determine the z-coordinate of the outer leaflet phosphate plane. The outer leaflet of G-OM contains ECLIPA; its phosphate groups are the topmost layer. Estimate from the membrane GRO or from the original box (z ≈ 65 Å above the box midpoint for a ~130 Å tall box):

   ```bash
   # Find approximate z of outer leaflet phosphate atoms (e.g. PB atom in ECLI)
   $ grep " PB " membrane.gro | awk '{print $NF}' | sort -n | tail -5
    11.666
    11.667
    11.668
    11.668
    11.669
   ```

   Translate the PMB PDB so its acyl tail (MOA1, atom 1 in `ligandrm.pdb`) sits ~15–20 Å above that z value ($11.669+2=13.669$):

   ```bash
   gmx editconf -f ../PMB_topol/charmm-gui-7883126705/ligandrm.pdb \
       -o pmb_placed.gro -translate 0 0 13.669
   ```

   (GROMACS uses nm; divide the Å value by 10.)

   **2c. Merge coordinates** by concatenating the atom blocks. The simplest approach is a Python or shell merge:

   ```bash
   # Strip END/ENDMDL lines, concatenate, then fix atom count on line 2
   python3 - <<'EOF'
   import re

   with open("membrane.gro") as f:
       lines = f.readlines()
   title = lines[0]
   natoms_mem = int(lines[1].strip())
   mem_atoms = lines[2:2+natoms_mem]
   box = lines[2+natoms_mem]

   with open("pmb_placed.gro") as f:
       plines = f.readlines()
   natoms_pmb = int(plines[1].strip())
   pmb_atoms = plines[2:2+natoms_pmb]

   with open("system.gro", "w") as out:
       out.write(title)
       out.write(f"{natoms_mem + natoms_pmb}\n")
       out.writelines(mem_atoms)
       out.writelines(pmb_atoms)
       out.write(box)
   EOF
   ```
   
   The merge script (`merge_coordinates.py`) copies the box vector line from `membrane.gro` into `system.gro`, so using `step5_input.gro` as the source guarantees the correct box is propagated.

   **2d. Update the topology** — add `PMB.itp` and one PMB molecule entry to `topol.top`:

   ```bash
   # Copy the validated PMB itp and its charmm36 parameter file into the simulation directory
   cp ../PMB_topol/charmm-gui-7883126705/gromacs/PMB.itp .
   cp ../PMB_topol/charmm-gui-7883126705/gromacs/charmm36.itp toppar/
   ```

   The membrane `toppar/forcefield.itp` is generated for lipid atom types only and does not contain the 16 CGenFF atom types used by PMB (`OG2D1`, `NG2S1`, `CG2O1`, etc.). Copying `charmm36.itp` from the PMB Ligand Reader output supplies those missing types — all PMB atom types are absent from the membrane forcefield so there are no duplicate definitions.

   However, `charmm36.itp` as generated by CHARMM-GUI contains a `[ defaults ]` section at the top, which will cause a fatal `grompp` error (`Invalid order for directive defaults`) because `forcefield.itp` already defines it. Strip it before use:

   ```bash
   sed -i '/^\[ defaults \]/,/^$/d' toppar/charmm36.itp
   ```

   Verify the `[ defaults ]` block is gone and `[ atomtypes ]` is the first directive:

   ```bash
   grep -n "^\[" toppar/charmm36.itp | head -3
   # Expected: first match is [ atomtypes ]
   ```

   Edit `topol.top`:
   - Immediately after `#include "toppar/forcefield.itp"` (before any lipid molecule `.itp` includes), add:
     ```
     #include "toppar/charmm36.itp"
     ```
     **Order is critical:** `[ atomtypes ]` must precede all `[ moleculetype ]` sections. The lipid `.itp` files (LPSA, PVPG, etc.) contain `[ moleculetype ]`, so placing `charmm36.itp` after them causes a fatal `Invalid order for directive atomtypes` error.
   - After the last lipid/water `#include` line, add:
     ```
     #include "PMB.itp"
     ```
   - At the end of the `[ molecules ]` block, add:
     ```
     PMB               1
     ```

   #### Step A3 — Verify the combined topology

   ```bash
   # Check total PMB charge
   awk '/^\[ atoms \]/{a=1;next} /^\[ bonds \]/{a=0} a && $1 ~ /^[0-9]+$/ {q+=$7} END {printf "PMB total charge = %.6f\n", q}' PMB.itp
   ```

   The PMB charge must be `+5.000000`.

   **Generate the index file.** The CHARMM-GUI equilibration MDP files reference two custom thermostat/COM groups — `MEMB` (all lipid, LPS glycan, and PMB atoms) and `SOLV` (water + ions) — that are not in the default `make_ndx` output. Use the provided script to build them:

   ```bash
   python3 make_index.py
   # Prints: MEMB: N atoms, SOLV: M atoms, Total: N+M / natoms
   # WARNING lines mean unknown residue names — add them to MEMB_RES or SOLV_RES in the script
   ```

   > **Note:** ECLIPA (LPSA in `topol.top`) is split into multiple glycan residue names in the GRO (`ECLIP`, `ARHM`, `AGAL`, `BMANN`, `BGLCN`, `AHEP`, `AKDO`, `AGLC`, `BGLC`); all are assigned to `MEMB`.

   **Neutralise the system.** The pre-existing ions in the membrane system were balanced for lipid charges only. PMB adds +5, leaving the system with a net charge of +5. Ewald electrostatics (PME) requires a neutral system; running with net charge causes a uniform background charge correction and serious artefacts. Add 5 Cl⁻ counter-ions by replacing water molecules.

   First do a trial `grompp` to produce a TPR for `genion`:

   ```bash
   gmx grompp -f step6.1_equilibration.mdp -c system.gro -r system.gro \
       -p topol.top -n index.ndx -o genion.tpr -maxwarn 5
   ```

   Then add the neutralising ions. Do **not** pass `-n index.ndx` here — `genion` requires a single molecule type and will error if given a combined group like `SOLV`. Use the default groups instead and select `TIP3` (water) when prompted:

   ```bash
   gmx genion -s genion.tpr -o system.gro -p topol.top \
       -pname SOD -nname CLA -neutral
   # At the prompt: type the number corresponding to TIP3 and press Enter
   ```

   `topol.top` is updated automatically (5 fewer TIP3, 5 more CLA). Regenerate the index file since atom numbers have changed:

   ```bash
   python3 make_index.py
   ```

   **Quick sanity check** — `grompp` should now complete without the net-charge warning. All three extra flags are required:
   - `-r system.gro` — reference coordinates for position restraints (all equilibration MDP files use `-DPOSRES`)
   - `-n index.ndx` — the custom `MEMB`/`SOLV` groups referenced by `tc_grps` and `comm_grps`

   ```bash
   gmx grompp -f step6.1_equilibration.mdp -c system.gro -r system.gro \
       -p topol.top -n index.ndx -o test.tpr -maxwarn 5
   ```

   #### Step A4 — Run the GROMACS equilibration pipeline

   CHARMM-GUI generates six staged GROMACS equilibration `.mdp` files that progressively relax position restraints on lipid head groups, tails, and the LPS glycan rings. Because PMB is newly placed, these stages allow it to relax into position without disrupting lipid packing.

   **Equilibration protocol summary (from `simulation/step6.x_equilibration.mdp`):** The six stages total ~1.875 ns of MD. Steps 6.1 and 6.2 (125 ps each, dt = 1 fs) are **NVT** runs using V-rescale temperature coupling at 310.15 K applied separately to the `MEMB` and `SOLV` groups; velocities are generated at the start of step 6.1. Steps 6.3–6.6 switch to **NPT** with C-rescale semi-isotropic pressure coupling (τ_p = 5 ps, P_ref = 1 bar), and the timestep is doubled to 2 fs for 500 ps per stage. Lipid head-group position restraints (`POSRES_FC_LIPID`, kJ mol⁻¹ nm⁻²) and dihedral restraints (`DIHRES_FC`, kJ mol⁻¹ rad⁻²) are relaxed progressively across the six stages: 1000/1000 → 400/400 → 400/200 → 200/200 → 40/100 → 0/0 (unrestrained). PMB carries no `POSRES` block in its ITP and is therefore unconstrained throughout all equilibration stages. All stages use the Verlet cutoff scheme, PME electrostatics with r_c = 1.2 nm, a force-switched Lennard-Jones cutoff (r_switch = 1.0 nm, r_c = 1.2 nm), and LINCS constraints on all hydrogen bonds.

   For better PP, I chose to run the EM locally on workstn:
   ```bash
   # Energy minimization
   gmx grompp \
       -f step6.0_minimization.mdp \
       -c system.gro -r system.gro \
       -p topol.top \
       -n index.ndx \
       -o em.tpr -maxwarn 5
   gmx mdrun -v -deffnm em -nt 12
   ``` 
   Check for convergence of potential energy
   ```bash
   echo "Potential 0" | gmx energy -f em.edr -o em_potential.xvg
   ```

   ```python
    import matplotlib.pyplot as plt
    x, y = [], []
    with open('em_potential.xvg') as f:
        for line in f:
            if line.startswith(('#','@')): continue
            cols = line.split()
            x.append(float(cols[0])); y.append(float(cols[1]))
    plt.plot(x, y)
    plt.xlabel('Step'); plt.ylabel('Potential Energy (kJ/mol)')
    plt.title('EM convergence')
    plt.tight_layout(); plt.show()
   ```


   EM plateaued cleanly before step 2000, well within the 5000-step limit — the EM fully converged to around $-1.1\times 10^{6}\;kJ/mol$, or around $-9.6\;kJ/mol/atom$.

   Prepare a publishable image of this initial condition by removing PBC artifacts 
  
   ```bash
   # Make molecules whole (repair PBC breaks), center on MEMB
    gmx trjconv -s em.tpr -f em.gro -o em_whole.gro \
    -pbc mol -center -n index.ndx 
   ```

   Then transfer the output (`em.gro, em.edr, em.log`) to HPC and start the equilibrations/production from there. 

   > **Note on `-r` (reference coordinates):** Following the CHARMM-GUI [README](simulation/README) convention, all equilibration stages use `system.gro` as the fixed `-r` reference, not the output of the previous stage. Position restraints should always pull atoms back to the original well-placed starting structure; using `${prev}.gro` as reference would let the restraint target drift each stage, defeating the purpose of the progressive force-constant schedule.

   > **Note on PMB restraints:** PMB is **unrestrained throughout all equilibration stages**. `PMB.itp` contains no `#ifdef POSRES` block, so the `-DPOSRES` flag has no effect on PMB. This is intentional: PMB is placed ~2.76 nm above the outer leaflet and is free to approach the membrane naturally. The equilibration restrains only the lipid head groups and ring dihedrals (via `POSRES_FC_LIPID` and `DIHRES_FC` in the lipid ITPs).

   ```bash
   nprocs=48

   # Staged equilibration (steps 6.1 through 6.6)
   for step in 6.1 6.2 6.3 6.4 6.5 6.6; do
       prev=$(echo $step | awk -F. '{if ($2==1) print "em"; else printf "step%s.%s", $1, $2-1}')
       echo ">>> Starting equilibration step${step} (prev: ${prev})"
       gmx grompp \
           -f step${step}_equilibration.mdp \
           -c ${prev}.gro -r system.gro \
           -p topol.top \
           -n index.ndx \
           -o step${step}.tpr -maxwarn 5
       gmx mdrun -v -deffnm step${step} -nb gpu -pme gpu -bonded gpu -nt ${nprocs} -gpu_id 0
       echo ">>> Completed equilibration step${step}. Checkpoint: step${step}.cpt"
   done
   echo ">>> All equilibration stages complete."

   # Production MD — run in chunks to fit HPC queue time limits.
   # Each chunk restarts from the checkpoint of the previous one (-t flag).
   # Adjust set_cntmax to match your target simulation time given the per-chunk nsteps in step7_production.mdp.
   # To restart after interruption, set cnt_start to the last completed chunk number + 1.
   cnt_start=1
   set_cntmax=10
   for cnt in $(seq ${cnt_start} ${set_cntmax}); do
       pcnt=$((cnt - 1))
       istep=step7_${cnt}
       echo ">>> Starting production chunk ${cnt}/${set_cntmax}: ${istep}"
       if [ ${cnt} -eq 1 ]; then
           pstep=step6.6
           gmx grompp \
               -f step7_production.mdp \
               -c ${pstep}.gro \
               -p topol.top \
               -n index.ndx \
               -o ${istep}.tpr -maxwarn 2
       else
           pstep=step7_${pcnt}
           gmx grompp \
               -f step7_production.mdp \
               -c ${pstep}.gro -t ${pstep}.cpt \
               -p topol.top \
               -n index.ndx \
               -o ${istep}.tpr -maxwarn 2
       fi
       gmx mdrun -v -deffnm ${istep} -nb gpu -pme gpu -bonded gpu -nt ${nprocs} -gpu_id 0
       echo ">>> Completed production chunk ${cnt}/${set_cntmax}. Checkpoint: ${istep}.cpt"
   done
   echo ">>> All production chunks complete."
   ```

   On an HPC cluster, replace the `gmx mdrun` calls with the appropriate scheduler script (PBS/Slurm). Use NPT ensemble with semi-isotropic pressure coupling (`Pcoupltype = semiisotropic`) for a membrane system.

##### Simulation Times
For observing PMB insertion/translocation through a membrane, **1 ns is far too short**. Typical timescales:

| Goal | Time needed |
|---|---|
| PMB docking to outer leaflet | 10–100 ns |
| Partial insertion into LPS headgroup region | 100–500 ns |
| Full translocation through outer membrane | 1–10 µs (rarely seen in unbiased MD) |

**Practical recommendation: 100–500 ns** for a first test.

- At 100 ns you should see whether PMB adsorbs and begins inserting into the LPS/Lipid A headgroup region — the most relevant event for its antimicrobial mechanism
- At 500 ns you have a reasonable chance of observing deeper penetration

If you want to observe full translocation within a feasible compute budget, consider **enhanced sampling methods** instead of plain MD:
- **Steered MD (SMD)** — pull PMB through the membrane along Z with a constant force; fast but biased
- **Umbrella sampling** — compute the PMB translocation free energy profile (PMF) along Z
- **REST2 or metadynamics** — accelerate conformational sampling without a fixed reaction coordinate

**Chunked production runs:** Production is split into `set_cntmax` chunks (default 10), each driven by the same `step7_production.mdp`. The total simulation time is therefore:

$$t_{\text{total}} = \text{set\_cntmax} \times \text{nsteps} \times \text{dt}$$

For example, with `nsteps = 5000000` (10 ns/chunk, since 5,000,000 × 0.002 ps = 10,000 ps) and `set_cntmax = 10`, the total is 100 ns. Adjust `nsteps` in `step7_production.mdp` and `set_cntmax` in the run script to reach your target. To extend an existing run, simply increase `set_cntmax` and resubmit — the loop skips already-completed chunks if you set `cnt_start` to the next chunk number.

For a first exploratory run on a cluster, `step7_production.mdp` is already set to `nsteps = 5000000` (10 ns/chunk). With the default `set_cntmax = 10` this gives 100 ns total — enough to see whether PMB contacts the outer leaflet phosphates. That result will tell you whether full translocation is feasible with unbiased MD or whether enhanced sampling is needed.

##### Runtime Estimate
Rough estimate requires knowing the system atom count first. From the README:

| Component | Count | ~Atoms |
|---|---|---|
| ECLIPA (35 × ~280 atoms) | 35 lipids | ~9,800 |
| PPPE/PVPG/PVCL2 (100 lipids) | 100 lipids | ~5,500 |
| Water (14,934 molecules) | — | ~44,800 |
| Ions (Na⁺/Cl⁻/Ca²⁺) | ~285 | ~285 |
| PMB | 188 atoms | 188 |
| **Total** | | **~60,000–70,000 atoms** |

**Expected performance on V100 + 48 CPUs:** ~100–200 ns/day for a system this size with CHARMM36 (LPS topology is heavier than a simple bilayer due to complex sugar bonded terms).

| Target simulation | Estimated wall time |
|---|---|
| 100 ns | **12–24 hours** |
| 500 ns | **3–5 days** |

**Caveats:**
- Use `-nb gpu -pme gpu -bonded gpu` in `mdrun` to offload as much as possible to the V100
- If PME stays on CPU (common with 48 CPUs available), performance may land closer to the lower end
- Run a short benchmark first (`gmx mdrun -nsteps 5000 -resetstep 4000 -v`) — it will print `Performance: X ns/day` so you get an exact number before committing to a long queue submission

   #### Step A5 — Record PMB atom-number offset for analysis

   After topology generation, PMB will be the last molecule block. Find the first atom number of PMB in the final system:

   ```bash
   gmx make_ndx -f md.tpr -o interaction_groups.ndx
   # At the prompt, identify the PMB group and note the atom range
   ```

   Add the offset to the `PMB.itp` atom numbers in the fragment map (section 2.ii above) to get the system-level atom numbers for each PMB fragment. Save the final index groups to `interaction_groups.ndx` before post-processing.
