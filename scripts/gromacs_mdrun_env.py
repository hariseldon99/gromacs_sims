# Temporary workaround for GROMACS thread handling under Slurm and Biobb mdrun.
# Why:
# - On some Slurm clusters, resumed jobs inherit/restore inconsistent OMP_NUM_THREADS.
# - GROMACS 2024+ with GPUs requires explicit -ntmpi when -ntomp/GPU is used.
# - Setting OMP_NUM_THREADS globally can crash unrelated steps; we need per-run scoping.
# What this does:
# - Replicates Biobb's Mdrun staging/sandboxing/logging but launches gmx via subprocess
#   with OMP_NUM_THREADS scoped to the child process only.
# - Passes -ntmpi/-ntomp explicitly; avoids -nt to prevent conflicts.
# TODO (remove this file) when:
# - biobb_gromacs Mdrun natively scopes OMP env per run and/or handles Slurm resume
#   correctly with GROMACS 2024+, so the stock mdrun works reliably again.

import os
import subprocess
from typing import Optional
from biobb_common.tools import file_utils as fu
from biobb_gromacs.gromacs.mdrun import Mdrun as _BiobbMdrun

class MdrunOMPEnv(_BiobbMdrun):
    """Subclass of biobb Mdrun that runs gmx mdrun via subprocess with a scoped OMP_NUM_THREADS."""
    def launch(self) -> int:
        # Honour Biobb restart/skip semantics
        if self.check_restart():
            return 0

        # Ensure outputs staged like Biobb
        if not self.stage_io_dict["out"].get("output_trr_path"):
            self.stage_io_dict["out"]["output_trr_path"] = fu.create_name(prefix=self.prefix, step=self.step, name='trajectory.trr')
            self.tmp_files.append(self.stage_io_dict["out"]["output_trr_path"])

        # Stage files into a unique sandbox
        self.stage_files()

        # Build base command
        self.cmd = [self.binary_path, 'mdrun',
                    '-s', self.stage_io_dict["in"]["input_tpr_path"],
                    '-c', self.stage_io_dict["out"]["output_gro_path"],
                    '-e', self.stage_io_dict["out"]["output_edr_path"],
                    '-g', self.stage_io_dict["out"]["output_log_path"],
                    '-o', self.stage_io_dict["out"]["output_trr_path"]]

        # Optional inputs/outputs
        if self.stage_io_dict["in"].get("input_cpt_path"):
            self.cmd += ['-cpi', self.stage_io_dict["in"]["input_cpt_path"]]
        if self.stage_io_dict["out"].get("output_xtc_path"):
            self.cmd += ['-x', self.stage_io_dict["out"]["output_xtc_path"]]
        if self.stage_io_dict["out"].get("output_cpt_path"):
            self.cmd += ['-cpo', self.stage_io_dict["out"]["output_cpt_path"]]
            if getattr(self, 'checkpoint_time', None):
                self.cmd += ['-cpt', str(self.checkpoint_time)]
        if self.stage_io_dict["out"].get("output_dhdl_path"):
            self.cmd += ['-dhdl', self.stage_io_dict["out"]["output_dhdl_path"]]

        # MPI launcher passthrough (rare with thread-MPI, but keep parity)
        if self.mpi_bin:
            mpi_cmd = [self.mpi_bin]
            if self.mpi_np:
                mpi_cmd += ['-n', str(self.mpi_np)]
            if self.mpi_flags:
                mpi_cmd.extend(self.mpi_flags)
            self.cmd = mpi_cmd + self.cmd

        # Threading: prefer explicit -ntmpi/-ntomp; avoid -nt (total) to prevent conflicts
        # -ntmpi: enforce >=1 when GPU is used (GROMACS 2024+ requirement)
        mpi_ranks = None
        if getattr(self, 'num_threads_mpi', None):
            mpi_ranks = int(self.num_threads_mpi)
        if getattr(self, 'use_gpu', False):
            mpi_ranks = max(1, mpi_ranks or 1)
        if mpi_ranks:
            fu.log(f'User added number of gmx mpi threads: {mpi_ranks}', self.out_log)
            self.cmd += ['-ntmpi', str(mpi_ranks)]

        if getattr(self, 'num_threads_omp', None):
            fu.log(f'User added number of gmx omp threads: {self.num_threads_omp}', self.out_log)
            self.cmd += ['-ntomp', str(self.num_threads_omp)]
        if getattr(self, 'num_threads_omp_pme', None):
            fu.log(f'User added number of gmx omp_pme threads: {self.num_threads_omp_pme}', self.out_log)
            self.cmd += ['-ntomp_pme', str(self.num_threads_omp_pme)]

        # GPU flags
        if getattr(self, 'use_gpu', False):
            fu.log('Adding GPU specific settings: -nb gpu -pme gpu', self.out_log)
            self.cmd += ["-nb", "gpu", "-pme", "gpu"]
        if getattr(self, 'gpu_id', None):
            fu.log(f'List of unique GPU device IDs available to use: {self.gpu_id}', self.out_log)
            self.cmd += ['-gpu_id', str(self.gpu_id)]
        if getattr(self, 'gpu_tasks', None):
            fu.log(f'List of GPU device IDs mapping PP tasks to devices: {self.gpu_tasks}', self.out_log)
            self.cmd += ['-gputasks', str(self.gpu_tasks)]

        # Extra GMX env
        if getattr(self, 'gmx_lib', None):
            self.env_vars_dict['GMXLIB'] = self.gmx_lib

        # Child-only environment; do not leak Slurm/global OMP settings
        env = os.environ.copy()
        ntomp = getattr(self, 'num_threads_omp', None)
        if ntomp:
            env['OMP_NUM_THREADS'] = str(ntomp)
        else:
            env.pop('OMP_NUM_THREADS', None)
        for k, v in getattr(self, 'env_vars_dict', {}).items():
            env[str(k)] = str(v)

        # Log command and run inside sandbox
        fu.log('Running (subprocess): ' + ' '.join(self.cmd), self.out_log)
        workdir = self.stage_io_dict.get("unique_dir") or os.getcwd()
        proc = subprocess.run(self.cmd, cwd=workdir, env=env, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        # Stream output into Biobb logs
        if proc.stdout:
            for line in proc.stdout.splitlines():
                fu.log(line, self.out_log)

        self.return_code = proc.returncode
        if self.return_code != 0:
            fu.log(f'mdrun failed with return code {self.return_code}', self.err_log)
            try:
                self.copy_to_host()
            finally:
                self.tmp_files.append(self.stage_io_dict.get("unique_dir"))
                self.remove_tmp_files()
            raise RuntimeError('gmx mdrun (subprocess) failed')

        # Normal teardown
        self.copy_to_host()
        self.tmp_files.append(self.stage_io_dict.get("unique_dir"))
        self.remove_tmp_files()
        self.check_arguments(output_files_created=True, raise_exception=False)
        return self.return_code


def mdrun_env(input_tpr_path: str, output_gro_path: str, output_edr_path: str,
              output_log_path: str, output_trr_path: Optional[str] = None, input_cpt_path: Optional[str] = None,
              output_xtc_path: Optional[str] = None, output_cpt_path: Optional[str] = None,
              output_dhdl_path: Optional[str] = None, properties: Optional[dict] = None, **kwargs) -> int:
    """Function-style wrapper matching biobb.mdrun signature, using MdrunOMPEnv."""
    # Ensure a dict so we can pass/update log paths if needed
    if properties is None:
        properties = {}
    # Let BiobbObject manage logs (uses properties['out_log']/['err_log'] if provided)
    return MdrunOMPEnv(input_tpr_path=input_tpr_path, output_trr_path=output_trr_path,
                       output_gro_path=output_gro_path, output_edr_path=output_edr_path,
                       output_log_path=output_log_path, input_cpt_path=input_cpt_path,
                       output_xtc_path=output_xtc_path, output_cpt_path=output_cpt_path,
                       output_dhdl_path=output_dhdl_path, properties=properties, **kwargs).launch()