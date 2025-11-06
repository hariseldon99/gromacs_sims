#!/usr/bin/env python
"""
# Protein Ligand Complex MD Setup based on tutorial using BioExcel Building Blocks (biobb)
### Based on the official Gromacs tutorial: http://www.mdtutorials.com/gmx/complex/index.html
***
This tutorial aims to illustrate the process of **setting up a simulation system** containing a **protein in complex with a ligand**, step by step, using the **BioExcel Building Blocks library (biobb)**. 
***
**Biobb modules** used:

 - [biobb_io](https://github.com/bioexcel/biobb_io): Tools to fetch biomolecular data from public databases.
 - [biobb_model](https://github.com/bioexcel/biobb_model): Tools to model macromolecular structures.
 - [biobb_chemistry](https://github.com/bioexcel/biobb_chemistry): Tools to manipulate chemical data.
 - [biobb_gromacs](https://github.com/bioexcel/biobb_gromacs): Tools to setup and run Molecular Dynamics simulations.
 - [biobb_analysis](https://github.com/bioexcel/biobb_analysis): Tools to analyse Molecular Dynamics trajectories.
 - [biobb_structure_utils](https://github.com/bioexcel/biobb_structure_utils):  Tools to modify or extract information from a PDB structure file.
 
***
### Pipeline steps:
 1. [Input Parameters](#input)
 2. [Fetching PDB Structure](#fetch)
 3. [Fix Protein Structure](#fix)
 4. [Create Protein System Topology](#top)
 5. [Create ligand system topology](#ligtop)
 6. [Preparing Ligand Restraints](#restraints)
 7. [Create new protein-ligand complex structure file](#complex)
 8. [Create new protein-ligand complex topology file](#complextop)
 9. [Create Solvent Box](#box)
 10. [Fill the Box with Water Molecules](#water)
 11. [Adding Ions](#ions)
 12. [Energetically Minimize the System](#min)
 13. [Equilibrate the System (NVT)](#nvt)
 14. [Equilibrate the System (NPT)](#npt)
 15. [Free Molecular Dynamics Simulation](#free)
 16. [Post-processing Resulting 3D Trajectory](#post)
"""
import argparse
import os
import sys
from mpi4py import MPI

from biobb_structure_utils.utils.extract_heteroatoms import extract_heteroatoms
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.cat_pdb import cat_pdb

from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_gromacs.gromacs.pdb2gmx import pdb2gmx

from biobb_chemistry.ambertools.reduce_add_hydrogens import reduce_add_hydrogens

from biobb_chemistry.babelm.babel_minimize import babel_minimize

from biobb_chemistry.acpype.acpype_params_gmx import acpype_params_gmx

from biobb_gromacs.gromacs.make_ndx import make_ndx

from biobb_gromacs.gromacs.genrestr import genrestr

from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
from biobb_structure_utils.utils.cat_pdb import cat_pdb
from biobb_gromacs.gromacs_extra.append_ligand import append_ligand

from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.solvate import solvate
from biobb_gromacs.gromacs.grompp import grompp
from biobb_gromacs.gromacs.genion import genion
from biobb_gromacs.gromacs.mdrun import mdrun
# Temporary workaround for Slurm resume + GROMACS 2024 GPU thread requirements:
# Use subprocess-backed mdrun_env to scope OMP_NUM_THREADS per-run and pass -ntmpi/-ntomp explicitly.
# TODO: Remove this import and revert to stock biobb mdrun calls once biobb_gromacs fixes this upstream.
from gromacs_mdrun_env import mdrun_env

from biobb_analysis.gromacs.gmx_energy import gmx_energy

from biobb_gromacs.gromacs.make_ndx import make_ndx

from biobb_analysis.gromacs.gmx_image import gmx_image

from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str
import shutil
from tabulate import tabulate
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel as ob



def deprotonate_pdb(input_file, output_file):
    """
    Remove hydrogen atoms from a PDB file.

    Parameters:
        input_file (str): Path to the input PDB file.
        output_file (str): Path where the deprotonated PDB file will be saved.

    The function reads the input PDB file line by line. For lines that start with "ATOM" 
    or "HETATM", it checks the element column (columns 77-78) to see if it is a hydrogen (H).
    In cases where the element field is missing or ambiguous, the atom name (columns 13-16)
    is examined. If a hydrogen is detected, the line is skipped; otherwise, the line is written 
    to the output file unchanged. Note that the PDB format is fixed-width, so this script assumes 
    that the file follows standard formatting.
    """
    with open(input_file, "r") as f:
        lines = f.readlines()

    with open(output_file, "w") as outf:
        for line in lines:
            # Only check lines that define atoms.
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract the element symbol (columns 77-78) if possible.
                element = line[76:78].strip() if len(line) >= 78 else ""
                # If the element symbol is missing, fallback to the atom name (columns 13-16).
                atom_name = line[12:16].strip() if len(line) >= 16 else ""
                # Determine if this atom is hydrogen:
                if element.upper() == "H" or (not element and atom_name and atom_name[0].upper() == "H"):
                    continue  # Skip writing hydrogen atoms.
            # Write all non-hydrogen lines (and non-ATOM/HETATM lines) to the output.
            outf.write(line)



def molecular_dynamics(complex, protonated=True):
    """
    Main function to set up a molecular dynamics simulation for a protein-ligand complex.
    This function performs the following steps:
    1. Extracts the ligand and protein from the input PDB file.
    2. Fixes the protein structure by modeling missing side-chain atoms and checking for issues.
    3. Creates the protein system topology using GROMACS tools.
    4. Creates the ligand system topology using AmberTools and Babel.
    5. Prepares ligand restraints for the simulation.
    6. Creates a new protein-ligand complex structure file.
    7. Creates a new protein-ligand complex topology file.
    8. Creates a solvent box around the complex.
    9. Fills the box with water molecules.
    10. Adds ions to neutralize the system.
    11. Performs energy minimization of the system.
    12. Equilibrates the system in NVT ensemble.
    13. Equilibrates the system in NPT ensemble.
    14. Runs free molecular dynamics simulation.
    15. Post-processes the resulting trajectory.
    """
    # Redirect standard output and error to a log file with a timestamped name
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    outdir = complex['outdir']
    os.makedirs(outdir, exist_ok=True)
    log_filename = os.path.join(complex['outdir'], f"simulation_log_{timestamp}.log")
    
    # Tabulate the complex dictionary
    complex_table = [
        ["Input Structure", complex['input_structure']],
        ["Ligand Code", complex['ligand_code']],
        ["Ligand Charge", complex['ligand_charge']],
        ["Output Directory", complex['outdir']],
        ["Number of Processors", complex['nprocs']],
        ["Use GPU", complex['usegpu']],
        ["GPU ID", complex['gpuid']],
        ["Energy Minimization Steps", complex['em_steps']],
        ["NVT Time (ps)", int(complex['nvt_steps'])*0.002],
        ["NPT Time (ps)", int(complex['npt_steps'])*0.002],
        ["MD Time (ns)", int(complex['md_steps'])*0.000002],
        ["Log File", log_filename]
    ]
    print("=" * 50)
    print("Start of Protein-Ligand Dynamics Simulation")
    print("=" * 50)
    print(tabulate(complex_table, headers=["Parameter", "Value"], tablefmt="grid"))
    
    
    # Pause for 5 seconds
    time.sleep(5)
    sys.stdout = open(log_filename, "w")
    sys.stderr = sys.stdout
    complex_pdbfile = complex['input_structure']
    
    complex_pdb = os.path.basename(complex['input_structure'])
    ligandCode = complex['ligand_code']
    mol_charge = complex['ligand_charge']
    
    # Copy the complex_pdbfile into the output directory
    shutil.copy(complex_pdbfile, os.path.join(outdir, complex_pdb))
    # Store the current working directory
    current_working_directory = os.getcwd()
    os.chdir(outdir)

    usegpu = complex['usegpu']
    nprocs = complex['nprocs']
    gpuid = complex['gpuid']

    # Ensure all outputs are directed to the specified 'outdir'
    proteinFile = "prot_" + complex_pdb
    ligandFile = ligandCode + '.pdb'
    complexFile = "prot_" + complex_pdb + '_' + ligandCode + '.pdb'

    prop = {
     'heteroatoms' : [{"name": ligandCode}]
    }
    
    extract_heteroatoms(input_structure_path=complex_pdb,
        output_heteroatom_path=ligandFile,
        properties=prop)

    extract_molecule(input_structure_path=complex_pdb,
        output_molecule_path=proteinFile)

    print(proteinFile, ligandFile, complexFile)
    cat_pdb(input_structure1=proteinFile,
        input_structure2=ligandFile,
        output_structure_path=complexFile)

    # ***
    # ## Fix protein structure
    # **Checking** and **fixing** (if needed) the protein structure:<br>
    # - **Modeling** **missing side-chain atoms**, modifying incorrect **amide assignments**, choosing **alternative locations**.<br>
    # - **Checking** for missing **backbone atoms**, **heteroatoms**, **modified residues** and possible **atomic clashes**.
    # 
    # ***
    # **Building Blocks** used:
    #  - [FixSideChain](https://biobb-model.readthedocs.io/en/latest/model.html#module-model.fix_side_chain) from **biobb_model.model.fix_side_chain**
    # ***

    fixed_pdb = proteinFile.removesuffix('.pdb')+'_fixed.pdb'

    fix_side_chain(input_pdb_path=proteinFile,
                output_pdb_path=fixed_pdb)

    # ***
    # ## Create protein system topology
    # 
    # **Extra Step** The protein needs to be deprotonated for this to work.
    # 

    # Deprotonate the fixed PDB file, but only if it is already protonated.
    deprotonated_pdb = fixed_pdb.removesuffix('.pdb')+'_deprotonated.pdb'
    if protonated:
        deprotonate_pdb(fixed_pdb, deprotonated_pdb)
    else:
        deprotonated_pdb = fixed_pdb    

    # **Building GROMACS topology** corresponding to the protein structure.<br>
    # Force field used in this tutorial is [**amber99sb-ildn**](https://dx.doi.org/10.1002%2Fprot.22711): AMBER **parm99** force field with **corrections on backbone** (sb) and **side-chain torsion potentials** (ildn). Water molecules type used in this tutorial is [**spc/e**](https://pubs.acs.org/doi/abs/10.1021/j100308a038).<br>
    # Adding **hydrogen atoms** if missing. Automatically identifying **disulfide bridges**. <br>
    # 
    # Generating two output files: 
    # - **GROMACS structure** (gro file)
    # - **GROMACS topology** ZIP compressed file containing:
    #     - *GROMACS topology top file* (top file)
    #     - *GROMACS position restraint file/s* (itp file/s)
    # ***
    # **Building Blocks** used:
    #  - [Pdb2gmx](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.pdb2gmx) from **biobb_gromacs.gromacs.pdb2gmx**

    output_pdb2gmx_gro = proteinFile.removesuffix('.pdb')+'_pdb2gmx.gro'
    output_pdb2gmx_top_zip = proteinFile.removesuffix('.pdb')+'_pdb2gmx_top.zip'
    prop = {
        'force_field' : 'amber99sb-ildn',
        'water_type': 'spce'
    }

    # Create and launch bb
    pdb2gmx(input_pdb_path=deprotonated_pdb,
            output_gro_path=output_pdb2gmx_gro,
            output_top_zip_path=output_pdb2gmx_top_zip,
            properties=prop, restart=True)

    # ## Create ligand system topology
    # **Building GROMACS topology** corresponding to the ligand structure.<br>
    # Force field used in this tutorial step is **amberGAFF**: [General AMBER Force Field](http://ambermd.org/antechamber/gaff.html), designed for rational drug design.<br>
    # - [Step 1](#ligandTopologyStep1): Add **hydrogen atoms** if missing.
    # - [Step 2](#ligandTopologyStep2): **Energetically minimize the system** with the new hydrogen atoms. 
    # - [Step 3](#ligandTopologyStep3): Generate **ligand topology** (parameters). 
    # ***
    # **Building Blocks** used:
    #  - [ReduceAddHydrogens](https://biobb-chemistry.readthedocs.io/en/latest/ambertools.html#module-ambertools.reduce_add_hydrogens) from **biobb_chemistry.ambertools.reduce_add_hydrogens**
    #  - [BabelMinimize](https://biobb-chemistry.readthedocs.io/en/latest/babelm.html#module-babelm.babel_minimize) from **biobb_chemistry.babelm.babel_minimize** 
    #  - [AcpypeParamsGMX](https://biobb-chemistry.readthedocs.io/en/latest/acpype.html#module-acpype.acpype_params_gmx) from **biobb_chemistry.acpype.acpype_params_gmx** 
    # ***
    # ### Step 1: Add **hydrogen atoms**
    # 
    # Note that this is not needed if ligand was already protonated during docking and not deprotonated after.

    # Create Ligand system topology, STEP 1
    # Reduce_add_hydrogens: add Hydrogen atoms to a small molecule (using Reduce tool from Ambertools package)
    # Import module
    
    # Create prop dict and inputs/outputs
    output_reduce_h =  ligandCode +'.reduce.H.pdb'
    prop = {
    'nuclear' : 'true'
}   
    # Note: reduce_add_hydrogens() (AmberTools) has been observed to produce
    # molecules with odd electron counts that break downstream tools (e.g. acpype).
    # Handle its failures or prefer alternative H-addition methods when possible.
    #If ligand is protonated, then do nothing
    if protonated:
        output_reduce_h = ligandFile
    else:
        try:
            # Use OpenBabel OBConversion API to add hydrogens and produce a PDB string.
            obconv = ob.OBConversion()
            obconv.SetInAndOutFormats("pdb", "pdb")
            mol_ob = ob.OBMol()
            if not obconv.ReadFile(mol_ob, ligandFile):
                raise RuntimeError("OpenBabel could not read ligand PDB")
            # Add hydrogens (keeps coordinates when present)
            mol_ob.AddHydrogens()
            # Attempt to perceive bond orders/assign coordinates if needed (best-effort)
            try:
                mol_ob.PerceiveBondOrders()
            except Exception:
                pass
            pdb_block = obconv.WriteString(mol_ob)
            if not pdb_block:
                raise RuntimeError("OpenBabel produced empty PDB output")

            # Force residue name to ligand code (first 3 chars, upper) for all ATOM/HETATM lines
            resname = ligandCode[:3].upper()
            out_lines = []
            for ln in pdb_block.splitlines():
                if ln.startswith(("ATOM  ", "HETATM")):
                    prefix = ln[:17] if len(ln) >= 17 else ln.ljust(17)
                    suffix = ln[20:] if len(ln) > 20 else ""
                    new_ln = prefix + ("{:<3}".format(resname)) + suffix
                    out_lines.append(new_ln)
                else:
                    out_lines.append(ln)
            with open(output_reduce_h, "w", encoding="utf-8") as fh:
                fh.write("\n".join(out_lines) + "\n")
            print("Ligand protonated with OpenBabel AddH (resname forced):", output_reduce_h)
        except Exception as ob_exc:
            print("OpenBabel AddH failed:", ob_exc)
            print("Falling back to reduce_add_hydrogens (AmberTools).")
            try:
                reduce_add_hydrogens(input_path=ligandFile,
                                    output_path=output_reduce_h,
                                    properties=prop)
                print("Ligand protonated with reduce_add_hydrogens:", output_reduce_h)
            except Exception as reduce_exc:
                raise RuntimeError("Failed to add hydrogens with OpenBabel and reduce_add_hydrogens: "
                                f"OpenBabel error: {ob_exc}; Reduce error: {reduce_exc}")
    
    # ### Step 2: **Energetically minimize the system** with the new hydrogen atoms. 

    # Create Ligand system topology, STEP 2
    # Babel_minimize: Structure energy minimization of a small molecule after being modified adding hydrogen atoms

    output_babel_min = ligandCode +'.H.min.mol2'                              
    prop = {
        'method' : 'sd',
        'criteria' : '1e-10',
        'force_field' : 'GAFF'
    }

    # Create and launch bb
    babel_minimize(input_path=output_reduce_h,
                output_path=output_babel_min,
                properties=prop)

    # ### Step 3: Generate **ligand topology** (parameters).

    # Create Ligand system topology, STEP 3
    # Acpype_params_gmx: Generation of topologies for GROMACS with ACPype
    # Import module

    # Create prop dict and inputs/outputs
    output_acpype_gro = ligandCode+'params.gro'
    output_acpype_itp = ligandCode+'params.itp'
    output_acpype_top = ligandCode+'params.top'
    output_acpype = ligandCode+'params'
    prop = {
        'basename' : output_acpype,
        'charge' : mol_charge
    }

    # Create and launch bb
    # Ensure the charge property is an integer to avoid issues with isinstance checks
    prop['charge'] = int(prop['charge'])

    acpype_params_gmx(input_path=output_babel_min, 
                    output_path_gro=output_acpype_gro,
                    output_path_itp=output_acpype_itp,
                    output_path_top=output_acpype_top,
                    properties=prop)

    # <a id="restraints"></a>
    # ***
    # ## Preparing Ligand Restraints
    # In subsequent steps of the pipeline, such as the equilibration stages of the **protein-ligand complex** system, it is recommended to apply some **restraints** to the small molecule, to avoid a possible change in position due to protein repulsion. **Position restraints** will be applied to the ligand, using a **force constant of 1000 KJ/mol\*nm^2** on the three coordinates: x, y and z. In this steps the **restriction files** will be created and integrated in the **ligand topology**.
    # - [Step 1](#restraintsStep1): Creating an index file with a new group including just the **small molecule heavy atoms**.
    # - [Step 2](#restraintsStep2): Generating the **position restraints** file.
    # ***
    # **Building Blocks** used:
    #  - [MakeNdx](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.make_ndx) from **biobb_gromacs.gromacs.make_ndx** 
    #  - [Genrestr](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.genrestr) from **biobb_gromacs.gromacs.genrestr** 
    # ***

    # <a id="restraintsStep1"></a>
    # ### Step 1: Creating an index file for the small molecule heavy atoms

    # MakeNdx: Creating import sys
    # Create prop dict and inputs/outputs
    output_ligand_ndx = ligandCode+'_index.ndx'
    prop = {
        'selection': "0 & ! a H*"
    }

    # Create and launch bb
    make_ndx(input_structure_path=output_acpype_gro,
            output_ndx_path=output_ligand_ndx,
            properties=prop)

    # <a id="restraintsStep2"></a>
    # ### Step 2: Generating the position restraints file

    # Genrestr: Generating the position restraints file

    # Create prop dict and inputs/outputs
    output_restraints_top = ligandCode+'_posres.itp'
    prop = {
        'force_constants': "1000 1000 1000",
        'restrained_group': "System_&_!H*"
    }

    # Create and launch bb
    genrestr(input_structure_path=output_acpype_gro,
            input_ndx_path=output_ligand_ndx,
            output_itp_path=output_restraints_top,
            properties=prop)


    # <a id="complex"></a>
    # ***
    # ## Create new protein-ligand complex structure file
    # Building new **protein-ligand complex** PDB file with:
    # - The new **protein system** with fixed problems from *Fix Protein Structure* step and hydrogens atoms added from *Create Protein System Topology* step.
    # - The new **ligand system** with hydrogens atoms added from *Create Ligand System Topology* step. 
    # 
    # This new structure is needed for **GROMACS** as it is **force field-compliant**, it **has all the new hydrogen atoms**, and the **atom names are matching the newly generated protein and ligand topologies**.
    # ***
    # **Building Blocks** used:
    #  - [GMXTrjConvStr](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_trjconv_str) from **biobb_analysis.gromacs.gmx_trjconv_str**
    #  - [CatPDB](https://biobb-structure-utils.readthedocs.io/en/latest/utils.html#module-utils.cat_pdb) from **biobb_structure_utils.utils.cat_pdb**
    # ***

    # Convert gro (with hydrogens) to pdb (PROTEIN)
    proteinFile_H = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_complex_H.pdb'
    prop = {
        'selection' : 'System'
    }

    # Create and launch bb
    gmx_trjconv_str(input_structure_path=output_pdb2gmx_gro,
                input_top_path=output_pdb2gmx_gro,
                output_str_path=proteinFile_H, 
                properties=prop)

    # Convert gro (with hydrogens) to pdb (LIGAND)
    ligandFile_H = ligandCode+'_complex_H.pdb'
    prop = {
        'selection' : 'System'
    }

    # Create and launch bb
    gmx_trjconv_str(input_structure_path=output_acpype_gro,
                input_top_path=output_acpype_gro,
                output_str_path=ligandFile_H, 
                properties=prop)

    # Concatenating both PDB files: Protein + Ligand
    complexFile_H = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_H.pdb'

    # Create and launch bb
    cat_pdb(input_structure1=proteinFile_H,
        input_structure2=ligandFile_H,
        output_structure_path=complexFile_H)

    # <a id="complextop"></a>
    # ***
    # ## Create new protein-ligand complex topology file
    # Building new **protein-ligand complex** GROMACS topology file with:
    # - The new **protein system** topology generated from *Create Protein System Topology* step.
    # - The new **ligand system** topology generated from *Create Ligand System Topology* step. 
    # 
    # NOTE: From this point on, the **protein-ligand complex structure and topology** generated can be used in a regular MD setup.
    # ***
    # **Building Blocks** used:
    #  - [AppendLigand](https://biobb-gromacs.readthedocs.io/en/latest/gromacs_extra.html#gromacs-extra-append-ligand-module) from **biobb_gromacs.gromacs_extra.append_ligand**
    # ***

    # AppendLigand: Append a ligand to a GROMACS topology
    # Import module

    # Create prop dict and inputs/outputs
    output_complex_top = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_complex.top.zip'

    posresifdef = "POSRES_"+ligandCode.upper()
    prop = {
        'posres_name': posresifdef
    }

    # Create and launch bb
    append_ligand(input_top_zip_path=output_pdb2gmx_top_zip,
                input_posres_itp_path=output_restraints_top,
                input_itp_path=output_acpype_itp, 
                output_top_zip_path=output_complex_top,
                properties=prop)


    # <a id="box"></a>
    # ***
    # ## Create solvent box
    # Define the unit cell for the **protein-ligand complex** to fill it with water molecules.<br> **Truncated octahedron** box is used for the unit cell. This box type is the one which best reflects the geometry of the solute/protein, in this case a **globular protein**, as it approximates a sphere. It is also convenient for the computational point of view, as it accumulates **less water molecules at the corners**, reducing the final number of water molecules in the system and making the simulation run faster.<br> A **protein to box** distance of **0.8 nm** is used, and the protein is **centered in the box**.  
    # 
    # ***
    # **Building Blocks** used:
    #  - [Editconf](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.editconf) from **biobb_gromacs.gromacs.editconf** 
    # ***

    # Editconf: Create solvent box
    # Import module

    # Create prop dict and inputs/outputs
    output_editconf_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_complex_editconf.gro'

    prop = {
        'box_type': 'octahedron',
        'distance_to_molecule': 0.8
    }

    # Create and launch bb
    editconf(input_gro_path=complexFile_H, 
            output_gro_path=output_editconf_gro,
            properties=prop)

    # <a id="water"></a>
    # ***
    # ## Fill the box with water molecules
    # Fill the unit cell for the **protein-ligand complex** with water molecules.<br>
    # The solvent type used is the default **Simple Point Charge water (SPC)**, a generic equilibrated 3-point solvent model. 
    # 
    # ***
    # **Building Blocks** used:
    #  - [Solvate](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.solvate) from **biobb_gromacs.gromacs.solvate** 
    # ***

    # Solvate: Fill the box with water molecules

    # Create prop dict and inputs/outputs
    output_solvate_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_solvate.gro'
    output_solvate_top_zip = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_solvate_top.zip'

    # Create and launch bb
    solvate(input_solute_gro_path=output_editconf_gro,
            output_gro_path=output_solvate_gro,
            input_top_zip_path=output_complex_top,
            output_top_zip_path=output_solvate_top_zip)

    # <a id="ions"></a>
    # ***
    # ## Adding ions
    # Add ions to neutralize the **protein-ligand complex** and reach a desired ionic concentration.
    # - [Step 1](#ionsStep1): Creating portable binary run file for ion generation
    # - [Step 2](#ionsStep2): Adding ions to **neutralize** the system and reach a **0.05 molar ionic concentration**
    # ***
    # **Building Blocks** used:
    #  - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp** 
    #  - [Genion](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.genion) from **biobb_gromacs.gromacs.genion** 
    # ***

    # <a id="ionsStep1"></a>
    # ### Step 1: Creating portable binary run file for ion generation

    # Grompp: Creating portable binary run file for ion generation

    # Create prop dict and inputs/outputs
    prop = {
        'mdp':{
            'nsteps':'5000',
        },
        'simulation_type':'ions',
        'maxwarn': 1
    }
    output_gppion_tpr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_complex_gppion.tpr'

    # Create and launch bb
    grompp(input_gro_path=output_solvate_gro,
        input_top_zip_path=output_solvate_top_zip, 
        output_tpr_path=output_gppion_tpr,
        properties=prop)

    # <a id="ionsStep2"></a>
    # ### Step 2: Adding ions to neutralize the system
    # Replace **solvent molecules** with **ions** to **neutralize** the system.

    # Genion: Adding ions to reach a 0.05 molar concentration

    # Create prop dict and inputs/outputs
    prop={
        'neutral':True,
        'concentration':0
    }
    output_genion_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_genion.gro'
    output_genion_top_zip = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_genion_top.zip'

    # Create and launch bb
    genion(input_tpr_path=output_gppion_tpr,
        output_gro_path=output_genion_gro, 
        input_top_zip_path=output_solvate_top_zip,
        output_top_zip_path=output_genion_top_zip, 
        properties=prop)

    # <a id="min"></a>
    # ***
    # ## Energetically minimize the system
    # Energetically minimize the **protein-ligand complex** till reaching a desired potential energy.
    # - [Step 1](#emStep1): Creating portable binary run file for energy minimization
    # - [Step 2](#emStep2): Energetically minimize the **protein-ligand complex** till reaching a force of 500 kJ/mol*nm.
    # - [Step 3](#emStep3): Checking **energy minimization** results. Plotting energy by time during the **minimization** process.
    # ***
    # **Building Blocks** used:
    #  - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp** 
    #  - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun** 
    #  - [GMXEnergy](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy) from **biobb_analysis.gromacs.gmx_energy** 
    # ***

    # <a id="emStep1"></a>
    # ### Step 1: Creating portable binary run file for energy minimization
    # Method used to run the **energy minimization** is a **steepest descent**, with a **maximum force of 500 kJ/mol\*nm^2**, and a minimization **step size of 1fs**. The **maximum number of steps** to perform if the maximum force is not reached is **5,000 steps**. 

    # Grompp: Creating portable binary run file for mdrun

    # Create prop dict and inputs/outputs
    prop = {
        'mdp':{
            'nsteps':complex['em_steps'],
            'emstep': 0.01,
            'emtol':'500'
        },
        'simulation_type':'minimization'
    }
    output_gppmin_tpr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_gppmin.tpr'

    # Create and launch bb
    grompp(input_gro_path=output_genion_gro,
        input_top_zip_path=output_genion_top_zip,
        output_tpr_path=output_gppmin_tpr,
        properties=prop)

    # <a id="emStep2"></a>
    # ### Step 2: Running Energy Minimization
    # Running **energy minimization** using the **tpr file** generated in the previous step.

    # Mdrun: Running minimization

    # Create prop dict and inputs/outputs
    output_min_trr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_min.trr'
    output_min_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_min.gro'
    output_min_edr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_min.edr'
    output_min_log = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_min.log'

    # Create and launch bb
    mdrun_env(input_tpr_path=output_gppmin_tpr,
        output_trr_path=output_min_trr,
        output_gro_path=output_min_gro,
        output_edr_path=output_min_edr,
        output_log_path=output_min_log,
        use_gpu=usegpu,
        num_threads_omp=nprocs,
        num_threads_mpi=1,
        gpu_id=gpuid)


    # <a id="emStep3"></a>
    # ### Step 3: Checking Energy Minimization results
    # Checking **energy minimization** results. Plotting **potential energy** by time during the minimization process. 

    # GMXEnergy: Getting system energy by time  

    # Create prop dict and inputs/outputs
    output_min_ene_xvg = proteinFile.removesuffix('.pdb') +'_'+ligandCode+'_min_ene.xvg'
    prop = {
        'terms':  ["Potential"]
    }

    # Create and launch bb
    gmx_energy(input_energy_path=output_min_edr, 
            output_xvg_path=output_min_ene_xvg, 
            properties=prop)

    # <a id="nvt"></a>
    # ***
    # ## Equilibrate the system (NVT)
    # Equilibrate the **protein-ligand complex** system in NVT ensemble (constant Number of particles, Volume and Temperature). To avoid temperature coupling problems, a new *"system"* group will be created including the **protein** + the **ligand** to be assigned to a single thermostatting group.
    # 
    # - [Step 1](#eqNVTStep1): Creating an index file with a new group including the **protein-ligand complex**.
    # - [Step 2](#eqNVTStep3): Creating portable binary run file for system equilibration
    # - [Step 3](#eqNVTStep3): Equilibrate the **protein-ligand complex** with NVT ensemble.
    # - [Step 4](#eqNVTStep4): Checking **NVT Equilibration** results. Plotting **system temperature** by time during the **NVT equilibration** process. 
    # ***
    # **Building Blocks** used:
    # - [MakeNdx](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.make_ndx) from **biobb_gromacs.gromacs.make_ndx** 
    # - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp** 
    # - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun** 
    # - [GMXEnergy](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy) from **biobb_analysis.gromacs.gmx_energy** 
    # ***

    # <a id="eqNVTStep1"></a>
    # ### Step 1: Creating an index file with a new group including the **protein-ligand complex**

    # MakeNdx: Creating index file with a new group (protein-ligand complex)

    # Create prop dict and inputs/outputs
    output_complex_ndx = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_index.ndx'
    prop = {
        'selection': "\"Protein\"|\"Other\"" 
    }

    # Create and launch bb
    make_ndx(input_structure_path=output_min_gro,
            output_ndx_path=output_complex_ndx,
            properties=prop)


    # <a id="eqNVTStep2"></a>
    # ### Step 2: Creating portable binary run file for system equilibration (NVT)
    # Note that for the purposes of temperature coupling, the **protein-ligand complex** (*Protein_Other*) is considered as a single entity.

    # Grompp: Creating portable binary run file for NVT System Equilibration

    # Create prop dict and inputs/outputs
    output_gppnvt_tpr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'gppnvt.tpr'
    prop = {
        'mdp':{
            'nsteps':complex['nvt_steps'],
            'tc-grps': 'Protein_Other Water_and_ions',
            'Define': '-DPOSRES -D' + posresifdef
        },
        'simulation_type':'nvt'
    }

    # Create and launch bb
    grompp(input_gro_path=output_min_gro,
        input_top_zip_path=output_genion_top_zip,
        input_ndx_path=output_complex_ndx,
        output_tpr_path=output_gppnvt_tpr,
        properties=prop, restart=True)


    # <a id="eqNVTStep3"></a>
    # ### Step 3: Running NVT equilibration

    # Mdrun: Running NVT System Equilibration

    # Create prop dict and inputs/outputs
    output_nvt_trr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_nvt.trr'
    output_nvt_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_nvt.gro'
    output_nvt_edr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_nvt.edr'
    output_nvt_log = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_nvt.log'
    output_nvt_cpt = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_nvt.cpt'

    # Create and launch bb
    # NOTE: Using mdrun_env (workaround). Revert to mdrun when upstream biobb handles OMP scoping/resume.
    mdrun_env(input_tpr_path=output_gppnvt_tpr,
        output_trr_path=output_nvt_trr,
        output_gro_path=output_nvt_gro,
        output_edr_path=output_nvt_edr,
        output_log_path=output_nvt_log,
        output_cpt_path=output_nvt_cpt,
        use_gpu=usegpu,
        num_threads_omp=nprocs,
        num_threads_mpi=1,
        gpu_id=gpuid)

    # <a id="eqNVTStep4"></a>
    # ### Step 4: Checking NVT Equilibration results
    # Checking **NVT Equilibration** results. Plotting **system temperature** by time during the NVT equilibration process. 

    # GMXEnergy: Getting system temperature by time during NVT Equilibration  

    # Create prop dict and inputs/outputs
    output_nvt_temp_xvg = proteinFile.removesuffix('.pdb') +'_'+ligandCode+'_nvt_temp.xvg'
    prop = {
        'terms':  ["Temperature"]
    }

    # Create and launch bb
    gmx_energy(input_energy_path=output_nvt_edr, 
            output_xvg_path=output_nvt_temp_xvg, 
            properties=prop)

    # <a id="npt"></a>
    # ***
    # ## Equilibrate the system (NPT)
    # Equilibrate the **protein-ligand complex** system in NPT ensemble (constant Number of particles, Pressure and Temperature) .
    # - [Step 1](#eqNPTStep1): Creating portable binary run file for system equilibration
    # - [Step 2](#eqNPTStep2): Equilibrate the **protein-ligand complex** with NPT ensemble.
    # - [Step 3](#eqNPTStep3): Checking **NPT Equilibration** results. Plotting **system pressure and density** by time during the **NPT equilibration** process.
    # ***
    # **Building Blocks** used:
    #  - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp** 
    #  - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun** 
    #  - [GMXEnergy](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_energy) from **biobb_analysis.gromacs.gmx_energy** 
    # ***


    # <a id="eqNPTStep1"></a>
    # ### Step 1: Creating portable binary run file for system equilibration (NPT)


    # Grompp: Creating portable binary run file for (NPT) System Equilibration

    # Create prop dict and inputs/outputs
    output_gppnpt_tpr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_gppnpt.tpr'
    prop = {
        'mdp':{
            'type': 'npt',
            'nsteps':complex['npt_steps'],
            'tc-grps': 'Protein_Other Water_and_ions',
            'Define': '-DPOSRES -D' + posresifdef
        },
        'simulation_type':'npt'
    }

    # Create and launch bb
    grompp(input_gro_path=output_nvt_gro,
        input_top_zip_path=output_genion_top_zip,
        input_ndx_path=output_complex_ndx,
        output_tpr_path=output_gppnpt_tpr,
        input_cpt_path=output_nvt_cpt,
        properties=prop, restart=True)


    # <a id="eqNPTStep2"></a>
    # ### Step 2: Running NPT equilibration


    # Mdrun: Running NPT System Equilibration

    # Create prop dict and inputs/outputs
    output_npt_trr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_npt.trr'
    output_npt_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_npt.gro'
    output_npt_edr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_npt.edr'
    output_npt_log = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_npt.log'
    output_npt_cpt = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_npt.cpt'

    # Create and launch bb
    # NOTE: Using mdrun_env (workaround). Revert to mdrun when upstream biobb handles OMP scoping/resume.
    mdrun_env(input_tpr_path=output_gppnpt_tpr,
        output_trr_path=output_npt_trr,
        output_gro_path=output_npt_gro,
        output_edr_path=output_npt_edr,
        output_log_path=output_npt_log,
        output_cpt_path=output_npt_cpt,
        use_gpu=usegpu,
        num_threads_omp=nprocs,
        num_threads_mpi=1,
        gpu_id=gpuid)


    # <a id="eqNPTStep3"></a>
    # ### Step 3: Checking NPT Equilibration results
    # Checking **NPT Equilibration** results. Plotting **system pressure and density** by time during the **NPT equilibration** process. 

    # GMXEnergy: Getting system pressure and density by time during NPT Equilibration  

    # Create prop dict and inputs/outputs
    output_npt_pd_xvg = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_npt_PD.xvg'
    prop = {
        'terms':  ["Pressure","Density"]
    }

    # Create and launch bb
    gmx_energy(input_energy_path=output_npt_edr, 
            output_xvg_path=output_npt_pd_xvg, 
            properties=prop)


    # <a id="free"></a>
    # ***
    # ## Free Molecular Dynamics Simulation
    # Upon completion of the **two equilibration phases (NVT and NPT)**, the system is now well-equilibrated at the desired temperature and pressure. The **position restraints** can now be released. The last step of the **protein-ligand complex** MD setup is a short, **free MD simulation**, to ensure the robustness of the system. 
    # - [Step 1](#mdStep1): Creating portable binary run file to run a **free MD simulation**.
    # - [Step 2](#mdStep2): Run short MD simulation of the **protein-ligand complex**.
    # - [Step 3](#mdStep3): Checking results for the final step of the setup process, the **free MD run**. Plotting **Root Mean Square deviation (RMSd)** and **Radius of Gyration (Rgyr)** by time during the **free MD run** step. 
    # ***
    # **Building Blocks** used:
    #  - [Grompp](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.grompp) from **biobb_gromacs.gromacs.grompp** 
    #  - [Mdrun](https://biobb-md.readthedocs.io/en/latest/gromacs.html#module-gromacs.mdrun) from **biobb_gromacs.gromacs.mdrun** 
    #  - [GMXRms](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rms) from **biobb_analysis.gromacs.gmx_rms** 
    #  - [GMXRgyr](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_rgyr) from **biobb_analysis.gromacs.gmx_rgyr** 
    # ***


    # <a id="mdStep1"></a>
    # ### Step 1: Creating portable binary run file to run a free MD simulation

    # Grompp: Creating portable binary run file for mdrun
    # Create prop dict and inputs/outputs
    prop = {
        'mdp':{
            #'nsteps':'500000' # 1 ns (500,000 steps x 2fs per step)
            #'nsteps':'5000' # 10 ps (5,000 steps x 2fs per step)
            'nsteps':complex['md_steps'] # 10 ns (5000000 steps x 2fs per step)
        },
        'simulation_type':'free'
    }
    output_gppmd_tpr = proteinFile.removesuffix('.pdb')+'_'+ligandCode + '_gppmd.tpr'

    # Create and launch bb
    grompp(input_gro_path=output_npt_gro,
        input_top_zip_path=output_genion_top_zip,
        output_tpr_path=output_gppmd_tpr,
        input_cpt_path=output_npt_cpt,
        properties=prop)

    # <a id="mdStep2"></a>
    # ### Step 2: Running short free MD simulation

    # Mdrun: Running free dynamics

    # Create prop dict and inputs/outputs
    output_md_trr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md.trr'
    output_md_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md.gro'
    output_md_edr = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md.edr'
    output_md_log = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md.log'
    output_md_cpt = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md.cpt'
    output_md_xtc = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md.xtc'

    # Create and launch bb
    mdrun_env(input_tpr_path=output_gppmd_tpr,
        output_trr_path=output_md_trr,
        output_xtc_path=output_md_xtc,
        output_gro_path=output_md_gro,
        output_edr_path=output_md_edr,
        output_log_path=output_md_log,
        output_cpt_path=output_md_cpt,
        use_gpu=usegpu,
        num_threads_omp=nprocs,
        num_threads_mpi=1,
        gpu_id=gpuid)

    # <a id="post"></a>
    # ***
    # ## Post-processing and Visualizing resulting 3D trajectory
    # Post-processing and Visualizing the **protein-ligand complex system** MD setup **resulting trajectory** using **NGL**
    # - [Step 1](#ppStep1): *Imaging* the resulting trajectory, **stripping out water molecules and ions** and **correcting periodicity issues**.
    # - [Step 2](#ppStep2): Generating a *dry* structure, **removing water molecules and ions** from the final snapshot of the MD setup pipeline.
    # - [Step 3](#ppStep3): Visualizing the *imaged* trajectory using the *dry* structure as a **topology**. 
    # ***
    # **Building Blocks** used:
    #  - [GMXImage](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_image) from **biobb_analysis.gromacs.gmx_image** 
    #  - [GMXTrjConvStr](https://biobb-analysis.readthedocs.io/en/latest/gromacs.html#module-gromacs.gmx_trjconv_str) from **biobb_analysis.gromacs.gmx_trjconv_str** 
    # ***


    # <a id="ppStep1"></a>
    # ### Step 1: *Imaging* the resulting trajectory.
    # Stripping out **water molecules and ions** and **correcting periodicity issues**  


    # GMXImage: "Imaging" the resulting trajectory
    #           Removing water molecules and ions from the resulting structure

    # Create prop dict and inputs/outputs (their way, which didn't work too well)
    #output_imaged_traj = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_imaged_traj.xtc'
    #prop = {
    #    'center_selection':  'Protein_Other',
    #    'output_selection': 'Protein_Other',
    #    'pbc' : 'mol',
    #    'center' : True
    #}
    
    # Create prop dict and inputs/outputs (my way, which worked better)
    output_cluster_traj = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_cluster_traj.xtc'
    prop = {
        'cluster_selection':  'Protein_Other',
        'output_selection': 'Protein_Other',
        'center' : False,
        'center_selection': 'Protein_Other', #This is due to a bug in biobb_gromacs     
        'pbc' : 'cluster'
    }

    # Create and launch bb 
    gmx_image(input_traj_path=output_md_xtc,
            input_top_path=output_gppmd_tpr,
            input_index_path=output_complex_ndx,
            output_traj_path=output_cluster_traj, 
            properties=prop)
    
    output_center_traj = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_cluster_center_traj.xtc'
    prop = {
        'center_selection':  'Protein_Other',
        'output_selection': 'Protein_Other',
        'pbc' : 'none',
        'center' : True
    }
    # Create and launch bb 
    gmx_image(input_traj_path=output_cluster_traj,
            input_top_path=output_gppmd_tpr,
            input_index_path=output_complex_ndx,
            output_traj_path=output_center_traj, 
            properties=prop)


    # <a id="ppStep2"></a>
    # ### Step 2: Generating the output *dry* structure.
    # **Removing water molecules and ions** from the resulting structure

    # GMXTrjConvStr: Converting and/or manipulating a structure
    #                Removing water molecules and ions from the resulting structure
    #                The "dry" structure will be used as a topology to visualize 
    #                the "imaged dry" trajectory generated in the previous step.

    # Create prop dict and inputs/outputs
    output_dry_gro = proteinFile.removesuffix('.pdb')+'_'+ligandCode+'_md_dry.gro'
    prop = {
        'selection':  'Protein_Other'
    }

    # Create and launch bb
    gmx_trjconv_str(input_structure_path=output_md_gro,
                    input_top_path=output_gppmd_tpr,
                    input_index_path=output_complex_ndx,
                    output_str_path=output_dry_gro, 
                    properties=prop) 
    
    os.chdir(current_working_directory)
    # Reset standard output and standard error back to default
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    print("=" * 50)
    print("End of Protein-Ligand Dynamics Simulation")
    print("=" * 50)
    return True

if __name__ == '__main__':
    # Execute main function
    parser = argparse.ArgumentParser(description="Protein-Ligand Molecular Dynamics Simulation Setup")
    parser.add_argument("--input_structure", type=str, required=True, help="Path to the input PDB structure file")
    parser.add_argument("--ligand_code", type=str, required=True, help="Ligand code (e.g., STM)")
    parser.add_argument("--ligand_charge", type=int, default=0, help="Ligand charge (default: 0)")
    parser.add_argument("--outdir", type=str, default="./outputs", help="Output directory (default: ./outputs)")
    parser.add_argument("--nprocs", type=int, default=int(os.environ.get('PBS_NCPUS', '12')), help="Number of processors (default: environment variable PBS_NCPUS or 12)")
    parser.add_argument("--usegpu", action='store_true', help="Use GPU for simulation (default: False)")
    parser.add_argument("--gpuid", type=str, default="0", help="GPU ID to use (default: '0')")
    parser.add_argument("--protonated", action='store_true', help="Set this option if the complex is protonated (default: False)")
    parser.add_argument("--em_steps", type=str, default='15000', help="Max number of Energy minimization steps (default: 15000)")
    parser.add_argument("--npt_steps", type=str, default='50000', help="Number of time steps for NPT equilibration (default: 5000 @ 2 fs per step)")
    parser.add_argument("--nvt_steps", type=str, default='50000', help="Number of time steps for  NVT equilibration (default: 50000 @ 2 fs per step")
    parser.add_argument("--md_steps", type=str, default='5000000', help="Number of time steps for production MD simulation (default: 5000000 @ 2 fs per step")

    args = parser.parse_args()
    
    complex = {
        'input_structure': args.input_structure,
        'ligand_code': args.ligand_code,
        'ligand_charge': args.ligand_charge,
        'outdir': args.outdir,
        'nprocs': args.nprocs,
        'usegpu': args.usegpu,
        'gpuid': args.gpuid,
        'em_steps': args.em_steps,
        'npt_steps': args.npt_steps,
        'nvt_steps': args.nvt_steps,
        'md_steps': args.md_steps
    }

    molecular_dynamics(complex, protonated=args.protonated)


# <a id="output"></a>
# ## Output files
# 
# Important **Output files** generated:
#  - **output_md_gro** : **Final structure** (snapshot) of the MD setup protocol.
#  - **output_md_trr** : **Final trajectory** of the MD setup protocol.
#  - **output_md_xtc** : **Final trajectory** of the MD setup protocol (single precision positions only).
#  - **output_md_cpt** : **Final checkpoint file**, with information about the state of the simulation. It can be used to **restart** or **continue** a MD simulation.
#  - **output_gppmd_tpr** : **Final tpr file**, GROMACS portable binary run input file. This file contains the starting structure of the **MD setup free MD simulation step**, together with the molecular topology and all the simulation parameters. It can be used to **extend** the simulation.
#  - **output_genion_top_zip** : **Final topology** of the MD system. It is a compressed zip file including a **topology file** (.top) and a set of auxiliary **include topology** files (.itp).
