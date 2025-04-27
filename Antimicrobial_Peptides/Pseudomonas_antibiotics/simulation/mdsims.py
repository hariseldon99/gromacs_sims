#!/usr/bin/env python

import importmonkey
import os
from mpi4py import MPI
from Bio.PDB import PDBParser
from math import ceil
importmonkey.add_path("/home/arunima/gitrepos/gromacs_sims/scripts")
import biobb_protein_ligand_simulation as bb

def chunk_into_n(lst, n):
  size = ceil(len(lst) / n)
  return list(
    map(lambda x: lst[x * size:x * size + size],
    list(range(n)))
  )

complexes_dir = "./simulation_complexes"

complex_pdbs = [
     "pbppen1999.pdb", "aminostep2016.pdb", "levgyra2018.pdb"
    ,"pbpmez2018.pdb", "pbppen2018.pdb"
]

complexes = []
# MPI initialization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

if rank == 0:
    for pdb in complex_pdbs:
        pdb_file = os.path.join(complexes_dir, pdb)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('complex', pdb_file)

        ligand_code = None
        for model in structure:
            for chain in model:
                for residue in chain:
                    if not residue.id[0].startswith(' '):  # Non-standard residues
                        ligand_code = residue.resname
                        break
                if ligand_code is not None:
                    break
            if ligand_code is not None:
                break
        print(f"Found ligand code: {ligand_code} in {pdb}")
        if ligand_code is None:
            print(f"No ligand code found in {pdb}")
            raise ValueError(f"Ligand code could not be determined for {pdb}. Exiting.")  

        complex = {
            'input_structure': pdb_file,
            'ligand_code': ligand_code,
            'ligand_charge': '0',
            'outdir': pdb.removesuffix('.pdb'),
            'nprocs': os.environ.get('OMP_NUM_THREADS', '24'),
            'mpithreads': '1',
            'usegpu': True,
            'gpuid': '0',
            'em_steps': '15000',
            'npt_steps': '50000',
            'nvt_steps': '50000',
            'md_steps': '5000000'
        }
        complexes.append(complex)
    print(f"Total complexes: {len(complexes)}")

if rank == 0:
    print(f"Running on {size} processes.")
# Distributing the complexes among the MPI processes

# Scatter the receptor_pdbqts list across the communicator
if rank == 0:
    complexes_loc = chunk_into_n(complexes, size)
else:
    complexes_loc = None

# Scatter the complexes list
complexes_loc = comm.scatter(complexes_loc, root=0)
# Each process will now have its own subset of complexes to work with
for complex in complexes_loc:
    bb.molecular_dynamics(complex, protonated=True)
