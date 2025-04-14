#!/usr/bin/env python3
import importmonkey 
importmonkey.add_path("/host_pwd")
receptors_dir = "/host_pwd/downloaded_models"
ligands_dir = "/host_pwd/ligands"
import blind_docking as bd
from multiprocessing import cpu_count
import os
nprocs_loc = os.environ.get("OMP_NUM_THREADS", cpu_count())
from mpi4py import MPI

from math import ceil
def chunk_into_n(lst, n):
  size = ceil(len(lst) / n)
  return list(
    map(lambda x: lst[x * size:x * size + size],
    list(range(n)))
  )

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if rank == 0:
    receptors = [os.path.join(receptors_dir, f) for f in os.listdir(receptors_dir) if f.endswith('.pdb')]
    print(f"Found {len(receptors)} receptors in {receptors_dir}.")

    ligands = [os.path.join(ligands_dir, f) for f in os.listdir(ligands_dir) if f.endswith('.pdb')]
    print(f"Found {len(ligands)} ligands in {ligands_dir}.")

    print("Preparing receptors and ligands for docking...")
    receptor_pdbqts = []
    #prepare receptors for docking
    for r in receptors:
        name = os.path.basename(r).split(".")[0]
        print(f"Preparing receptor {name} for docking...")
        receptor_pdbqts.append(bd.prepare_receptor(r))
        print(f"{name} prepared.")

    #prepare ligands for docking
    ligand_pdbqts = []
    for ligand in ligands:
        if ligand.endswith('.pdb'):
            ligand_name = os.path.basename(ligand).split(".")[0]
            print(f"Preparing {ligand_name} for docking...")
            ligand_pdbqts.append(bd.prepare_ligand(ligand))
            print(f"{ligand_name} prepared.")
        else:
            print(f"{ligand} is not a PDB file, skipping.")

    print("All receptors and ligands prepared for docking.")
    print("Starting batch docking...")
else:
    receptor_pdbqts = None
    ligands = None
# Broadcast the ligands list to all processes
ligands = comm.bcast(ligands, root=0)

# Scatter the receptor_pdbqts list across the communicator
if rank == 0:
    receptor_pdbqts_loc = chunk_into_n(receptor_pdbqts, size)
else:
    receptor_pdbqts_loc = None


# Scatter the receptor_pdbqts list
receptor_pdbqts_loc = comm.scatter(receptor_pdbqts_loc, root=0)

# Each process will now have a portion of the receptor_pdbqts list
for i, receptor in enumerate(receptor_pdbqts_loc):
    center_x, center_y, center_z = bd.calculate_center_of_mass(receptor)
    box_size = bd.estimate_box_size(receptor)
    size_x, size_y, size_z = box_size
    bd.batch_docking(receptor, ligands_dir, 
                  center_x, center_y, center_z, 
                  size_x, size_y, size_z,
                  nprocs_loc)
    print(f"Process {rank} completed docking for receptor {i+1}/{len(receptor_pdbqts_loc)}: {receptor}")
