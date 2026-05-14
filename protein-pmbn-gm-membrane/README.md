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
