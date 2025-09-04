#!/usr/bin/env python
import time
import os, sys
import glob
import subprocess

ext_time_ns = 40  # Extension time in nanoseconds

# Default values
numprocs = 24     # Number of OMP threads
mpithreads = 1    # Number of MPI ranks
gpuid = '0'       # Default GPU id

# Override with SLURM environment variables if present
numprocs = int(os.environ.get("SLURM_CPUS_PER_TASK", numprocs))
mpithreads = int(os.environ.get("SLURM_NTASKS", mpithreads))

# Detect GPU allocation from SLURM or CUDA_VISIBLE_DEVICES
if "CUDA_VISIBLE_DEVICES" in os.environ:
    visible = os.environ["CUDA_VISIBLE_DEVICES"].split(",")[0]
    gpuid = visible.strip()

usegpu = True if gpuid else False

simulation_dirnames = [
    "aminostep2012C", "aminotob2016C", "levparc2018C_chainB",
    "ndxgyra2018C", "ndxparc2018C_chainA", "pbppen1999C",
    "aminostep2016C", "levgyra2018C", "ndxgyra1996C",
    "ndxparc1999C_chainB", "pbpmez2018C", "pbppen2018C"
]

# Unset OMP_NUM_THREADS if set
omp_num_threads_backup = os.environ.get('OMP_NUM_THREADS', None)
if omp_num_threads_backup is not None:
    del os.environ['OMP_NUM_THREADS']

for simdir in simulation_dirnames:
    cwd = os.getcwd()
    target_dir = os.path.join(cwd, simdir)
    if not os.path.isdir(target_dir):
        raise FileNotFoundError(f"Directory '{simdir}' does not exist in {cwd}")
    os.chdir(target_dir)
    print(f"Entering directory {simdir}")

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    log_filename = f"extended_simulation_log_{timestamp}.log"
    sys.stdout = open(log_filename, "w")
    sys.stderr = sys.stdout

    print("="*50)
    print("Beginning of Extended Protein-Ligand Dynamics Simulation")
    print("="*50)

    input_tpr = glob.glob(f"*{simdir}*_gppmd.tpr")[0]
    input_cpt = glob.glob(f"*{simdir}*_md.cpt")[0]
    input_ndx = glob.glob(f"*{simdir}*_index.ndx")[0]
    output_gro = glob.glob(f"*{simdir}*_md.gro")[0].replace('.gro', '_ext.gro')
    output_trr = glob.glob(f"*{simdir}*_md.trr")[0].replace('.trr', '_ext.trr')
    output_edr = glob.glob(f"*{simdir}*_md.edr")[0].replace('.edr', '_ext.edr')
    output_log = glob.glob(f"*{simdir}*_md.log")[0].replace('.log', '_ext.log')
    output_xtc = glob.glob(f"*{simdir}*_md.xtc")[0].replace('.xtc', '_ext.xtc')
    output_cpt = glob.glob(f"*{simdir}*_md.cpt")[0].replace('.cpt', '_ext.cpt')
    output_dry_gro = glob.glob(f"*{simdir}*_md_dry.gro")[0].replace('_dry.gro', '_ext_dry.gro')

    # Step 1: Extend simulation by 40 ns (40000 ps)
    extended_tpr = input_tpr.replace('.tpr', '_ext.tpr')
    subprocess.run([
        "gmx", "convert-tpr",
        "-s", input_tpr,
        "-extend", str(ext_time_ns * 1000),
        "-o", extended_tpr
    ], check=True)

    # Step 2: Run mdrun with GPU (nonbonded + PME) and -noappend
    mdrun_cmd = [
        "gmx", "mdrun",
        "-s", extended_tpr,
        "-cpi", input_cpt,
        "-deffnm", simdir + "_ext",
        "-ntomp", str(numprocs),
        "-ntmpi", str(mpithreads),
        "-noappend"
    ]
    if usegpu:
        mdrun_cmd += ["-gpu_id", gpuid, "-nb", "gpu", "-pme", "gpu", "-update", "gpu"]
    subprocess.run(mdrun_cmd, check=True)

    # Step 3: Post-process with trjconv (cluster PBC)
    #Note: These seem to fail execution with no error messages. Needs debugging
    cluster_traj = output_xtc.replace('.xtc', '_cluster_traj.xtc')
    with subprocess.Popen(
        ["gmx", "trjconv", "-s", extended_tpr, "-n", input_ndx ,"-f", output_xtc, "-o", cluster_traj, "-pbc", "cluster"],
        stdin=subprocess.PIPE, text=True
    ) as proc:
        proc.communicate("Protein_Other\nProtein_Other\n")  # select group twice

    # Step 4: Center trajectory
    #Note: These seem to fail execution with no error messages. Needs debugging
    center_traj = output_xtc.replace('.xtc', '_center_traj.xtc')
    with subprocess.Popen(
        ["gmx", "trjconv", "-s", extended_tpr, "-n", input_ndx, "-f", cluster_traj, "-o", center_traj, "-center"],
        stdin=subprocess.PIPE, text=True
    ) as proc:
        proc.communicate("Protein_Other\n")

    print(f"Extended MD simulation completed for {simdir}")
    os.chdir(cwd)
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    print("="*50)
    print("End of Extended Protein-Ligand Dynamics Simulation")
    print("="*50)

# Restore OMP_NUM_THREADS
if omp_num_threads_backup is not None:
    os.environ['OMP_NUM_THREADS'] = omp_num_threads_backup