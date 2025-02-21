from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import MDAnalysis as mda
import numpy as np
import mdtraj as md
import pickle

pdbfile='Colpm_pyocyanin/colpm_pyc_complex/sys_prod/start.pdb'
xtcfile='Colpm_pyocyanin/colpm_pyc_complex/sys_prod/prod_colpm_pyc_complex_compact_compact.xtc'
tprfile="Colpm_pyocyanin/colpm_pyc_complex/sys_prod/prod_colpm_pyc_complex.tpr"
pickle_filename = 'colpm_pyc_hbdata.pkl'
lig_resname = 'PYC'


def do_hb(ligname,pdbfile, xtcfile, tprfile, stride = 3, verbose=True):
    # Load the trajectory using MDAnalysis to get the total number of frames
    u = mda.Universe(tprfile, xtcfile)
    total_frames = len(u.trajectory)
    total_time_ns = u.trajectory.totaltime / 1000  # Convert from ps to ns
    timestep = (stride / total_frames) * total_time_ns
    nframes = int(np.ceil(total_frames / stride))
    full_traj = md.load_xtc(xtcfile, stride=stride, top=pdbfile)
    lig_atoms = full_traj.topology.select(f'resname {ligname}')

    hbonds = md.baker_hubbard(full_traj, periodic=False, freq=0.01)
    label = lambda hbond : '%s - %s ... %s' % ( full_traj.topology.atom(hbond[0]),\
                                                full_traj.topology.atom(hbond[1]),\
                                                full_traj.topology.atom(hbond[2]))

    lig_hbonds = []
    for hbond in hbonds:
        d,h,a = hbond
        if d in lig_atoms or a in lig_atoms:
            lig_hbonds.append(hbond)
            if verbose:
                print(label(hbond))

    lig_hbonds = np.array(lig_hbonds)
    #Giving full tpr file to include bonding information
    complex_hbonds = []

    for hbond in lig_hbonds:
        d_ix, h_ix, a_ix = hbond
        complex_hbonds.append(HBA(
            universe=u,
            donors_sel=f"index {d_ix}",
            hydrogens_sel=f"index {h_ix}",
            acceptors_sel=f"index {a_ix}",
            update_selections=False
        ).run(verbose=verbose))
    return complex_hbonds


prot_lig_hb = do_hb(lig_resname),pdbfile, xtcfile, tprfile, stride=6, verbose=False)

with open(pickle_filename, 'wb') as file:
    pickle.dump(prot_lig_hb, file)