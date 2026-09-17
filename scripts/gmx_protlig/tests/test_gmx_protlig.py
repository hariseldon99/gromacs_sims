"""
Comprehensive unit tests for gmx_protlig package.
"""

import os
import sys
import tempfile
import unittest
import pickle
from pathlib import Path

# Ensure gmx_protlig is on python path
PACKAGE_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PACKAGE_ROOT not in sys.path:
    sys.path.insert(0, PACKAGE_ROOT)

import gmx_protlig
from gmx_protlig.utils import (
    read_gro,
    write_combined_gro,
    extract_atomtypes_block,
    extract_nonbonded_params_block,
    strip_global_directives,
    rewrite_include_paths,
    save_checkpoint,
    load_latest_checkpoint,
    relocate_checkpoint_paths,
    enable_gromacs_py_gpu_hook,
    disable_gromacs_py_gpu_hook,
    is_gpu_hook_enabled,
)
from gmx_protlig.batch import (
    BatchPipelineRunner,
    _parse_float,
    _parse_int,
    _parse_bool,
)
from gmx_protlig.pipeline import ProteinLigandPipeline


class TestGmxProtligImports(unittest.TestCase):
    def test_imports(self):
        self.assertIsNotNone(gmx_protlig.__version__)
        self.assertTrue(hasattr(gmx_protlig, "ProteinLigandPipeline"))
        self.assertTrue(hasattr(gmx_protlig, "BatchPipelineRunner"))
        self.assertTrue(hasattr(gmx_protlig, "step1_extract_and_prepare_ligand"))
        self.assertTrue(hasattr(gmx_protlig, "step8_hpc_equi_prod"))

    def test_protlig_api_compatibility(self):
        import protlig_api
        self.assertTrue(hasattr(protlig_api, "ProteinLigandPipeline"))
        self.assertTrue(hasattr(protlig_api, "BatchPipelineRunner"))
        self.assertTrue(hasattr(protlig_api, "step1_extract_and_prepare_ligand"))


class TestGroUtils(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.prot_gro = os.path.join(self.temp_dir.name, "prot.gro")
        self.lig_gro = os.path.join(self.temp_dir.name, "lig.gro")
        self.out_gro = os.path.join(self.temp_dir.name, "combined.gro")

        # Create mock protein GRO (2 atoms)
        with open(self.prot_gro, "w") as fh:
            fh.write("Protein System\n")
            fh.write("    2\n")
            fh.write("    1ALA      N    1   1.000   1.000   1.000\n")
            fh.write("    1ALA     CA    2   1.100   1.100   1.100\n")
            fh.write("   5.00000   5.00000   5.00000\n")

        # Create mock ligand GRO (1 atom)
        with open(self.lig_gro, "w") as fh:
            fh.write("Ligand System\n")
            fh.write("    1\n")
            fh.write("    1UNL     C1    1   2.000   2.000   2.000\n")
            fh.write("   5.00000   5.00000   5.00000\n")

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_read_gro(self):
        title, atoms, box = read_gro(self.prot_gro)
        self.assertEqual(title.strip(), "Protein System")
        self.assertEqual(len(atoms), 2)
        self.assertIn("5.00000", box)

    def test_write_combined_gro(self):
        total = write_combined_gro(self.prot_gro, self.lig_gro, self.out_gro, title="MergedComplex")
        self.assertEqual(total, 3)

        title, atoms, box = read_gro(self.out_gro)
        self.assertEqual(title.strip(), "MergedComplex")
        self.assertEqual(len(atoms), 3)
        # Verify atom renumbering: atom 3 should have index 3
        self.assertIn("3", atoms[2][15:20])


class TestTopologyUtils(unittest.TestCase):
    def test_extract_atomtypes_and_strip(self):
        itp_content = """\
[ moleculetype ]
; Name   nrexcl
UNL      3

[ atomtypes ]
; name  at.num  mass     charge ptype  sigma      epsilon
  c3    6       12.011   0.0000 A      0.339967   0.457730

[ atoms ]
; nr type resnr resid atom cgnr charge mass
  1  c3   1     UNL   C1   1    0.100  12.011

[ nonbond_params ]
; i j func sigma epsilon
  c3 c3 1 0.339967 0.457730
"""
        atomtypes = extract_atomtypes_block(itp_content)
        self.assertIn("[ atomtypes ]", atomtypes)
        self.assertIn("c3    6", atomtypes)

        nonbonded = extract_nonbonded_params_block(itp_content)
        self.assertIn("[ nonbond_params ]", nonbonded)

        stripped = strip_global_directives(itp_content)
        self.assertNotIn("[ atomtypes ]", stripped)
        self.assertNotIn("[ nonbond_params ]", stripped)
        self.assertIn("[ moleculetype ]", stripped)
        self.assertIn("[ atoms ]", stripped)

    def test_rewrite_include_paths(self):
        top_text = """\
#include "amber99sb-ildn.ff/forcefield.itp"
#include "posre.itp"
#include "amber99sb-ildn.ff/tip3p.itp"
"""
        with tempfile.TemporaryDirectory() as td:
            from_dir = os.path.join(td, "protein", "topol.top")
            to_dir = os.path.join(td, "complex")
            os.makedirs(os.path.join(td, "protein"), exist_ok=True)
            os.makedirs(to_dir, exist_ok=True)

            res = rewrite_include_paths(top_text, from_dir, to_dir)
            # System forcefield should NOT be rewritten
            self.assertIn('#include "amber99sb-ildn.ff/forcefield.itp"', res)
            # Local posre.itp should be rewritten to relative path ../protein/posre.itp
            self.assertIn('../protein/posre.itp', res)


class TestCheckpointUtils(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_save_and_load_checkpoint(self):
        obj = {"name": "test_system", "step": 5}
        path = save_checkpoint(obj, "test_step", search_dir=self.temp_dir.name)
        self.assertTrue(os.path.exists(path))

        loaded_obj, loaded_file = load_latest_checkpoint("test_step", search_dir=self.temp_dir.name)
        self.assertEqual(loaded_obj["name"], "test_system")
        self.assertEqual(loaded_file, path)

    def test_relocate_checkpoint_paths(self):
        class DummySys:
            def __init__(self, coor, top):
                self.coor_file = coor
                self.top_file = top

        # Create dummy existing files in new dir
        new_dir = os.path.join(self.temp_dir.name, "new_run")
        os.makedirs(new_dir, exist_ok=True)
        new_coor = os.path.join(new_dir, "complex.gro")
        with open(new_coor, "w") as f:
            f.write("gro")

        # System with stale old path from another host
        sys_obj = DummySys(coor="/old/cluster/path/complex.gro", top="non_path_string")
        relocate_checkpoint_paths(sys_obj, new_dir)
        self.assertEqual(sys_obj.coor_file, new_coor)
        self.assertEqual(sys_obj.top_file, "non_path_string")


class TestBatchAndManifests(unittest.TestCase):
    def test_manifest_parsers(self):
        self.assertEqual(_parse_float("1.5"), 1.5)
        self.assertEqual(_parse_float("", default=0.15), 0.15)
        self.assertIsNone(_parse_float(None))

        self.assertEqual(_parse_int("4"), 4)
        self.assertEqual(_parse_int("", default=1), 1)

        self.assertTrue(_parse_bool("True"))
        self.assertTrue(_parse_bool("yes"))
        self.assertTrue(_parse_bool("1"))
        self.assertFalse(_parse_bool("False"))
        self.assertFalse(_parse_bool("", default=False))

    def test_load_example_manifests(self):
        csv_path = os.path.join(PACKAGE_ROOT, "manifest_example.csv")
        if os.path.exists(csv_path):
            runner = BatchPipelineRunner(manifest=csv_path)
            self.assertEqual(len(runner.complexes), 2)
            self.assertEqual(runner.complexes[0]["ligand_name"], "Ligand1")
            self.assertEqual(runner.complexes[0]["residue_name"], "UNL")

        json_path = os.path.join(PACKAGE_ROOT, "manifest_example.json")
        if os.path.exists(json_path):
            runner_json = BatchPipelineRunner(manifest=json_path)
            self.assertEqual(len(runner_json.complexes), 2)
            self.assertEqual(runner_json.complexes[1]["ligand_name"], "Ligand2")
            self.assertEqual(runner_json.complexes[1]["residue_name"], "UNL")


class TestPdbqtParsing(unittest.TestCase):
    def test_multi_model_and_missing_endmdl(self):
        sample_pdbqt = """\
MODEL 1
REMARK VINA RESULT:    -7.8      0.000      0.000
ATOM      1  C1  UNL     1       1.000   1.000   1.000  0.00  0.00    +0.000 C
ENDMDL
MODEL 2
REMARK VINA RESULT:    -9.2      1.500      2.100
ATOM      1  C1  UNL     1       2.000   2.000   2.000  0.00  0.00    +0.000 C
MODEL 3
REMARK VINA RESULT:    -6.4      2.500      3.100
ATOM      1  C1  UNL     1       3.000   3.000   3.000  0.00  0.00    +0.000 C
"""
        with tempfile.TemporaryDirectory() as td:
            pdbqt_file = os.path.join(td, "docked.pdbqt")
            with open(pdbqt_file, "w") as f:
                f.write(sample_pdbqt)

            # Test parsing via helper logic
            with open(pdbqt_file, "r") as fh:
                raw = fh.read()

            blocks = {}
            current_model = None
            current_lines = []
            for line in raw.splitlines(keepends=True):
                if line.startswith("MODEL "):
                    if current_model is not None and current_lines:
                        blocks[current_model] = current_lines
                    current_model = int(line.split()[1])
                    current_lines = [line]
                elif line.startswith("ENDMDL") and current_model is not None:
                    current_lines.append(line)
                    blocks[current_model] = current_lines
                    current_model = None
                    current_lines = []
                elif current_model is not None:
                    current_lines.append(line)

            if current_model is not None and current_lines:
                blocks[current_model] = current_lines

            # All 3 models must be captured despite Model 3 lacking ENDMDL
            self.assertEqual(len(blocks), 3)

            # Check best model is Model 2 (-9.2)
            scores = {}
            for mid, lines in blocks.items():
                for l in lines:
                    if "VINA RESULT:" in l:
                        parts = l.split()
                        idx_res = next(i for i, p in enumerate(parts) if "RESULT" in p)
                        scores[mid] = float(parts[idx_res + 1])
                        break
            best_model = min(scores, key=scores.get)
            self.assertEqual(best_model, 2)
            self.assertEqual(scores[best_model], -9.2)


class TestGromacsPyGpuHook(unittest.TestCase):
    def setUp(self):
        disable_gromacs_py_gpu_hook()

    def tearDown(self):
        disable_gromacs_py_gpu_hook()

    def test_hook_toggle(self):
        self.assertFalse(is_gpu_hook_enabled())
        ok = enable_gromacs_py_gpu_hook()
        if not ok:
            # os_command_py not in current Python environment
            self.skipTest("os_command_py not available in environment")
        self.assertTrue(is_gpu_hook_enabled())
        disable_gromacs_py_gpu_hook()
        self.assertFalse(is_gpu_hook_enabled())

    def test_hook_em_omits_bonded_and_update(self):
        ok = enable_gromacs_py_gpu_hook(pin="on", pme="gpu", bonded="gpu", update="gpu")
        if not ok:
            self.skipTest("os_command_py not available in environment")

        import os_command_py.os_command as osc

        # EM command with "Init_em"
        cmd = osc.Command(["/usr/local/bin/gmx", "mdrun", "-s", "Init_em_cuedc2.tpr", "-deffnm", "Init_em_cuedc2"])
        cmd_str = " ".join(cmd.cmd)
        self.assertIn("-pin", cmd.cmd)
        self.assertIn("on", cmd.cmd)
        self.assertNotIn("-bonded", cmd.cmd)
        self.assertNotIn("-update", cmd.cmd)

    def test_hook_dynamical_md_includes_full_gpu(self):
        ok = enable_gromacs_py_gpu_hook(pin="on", pme="gpu", bonded="gpu", update="gpu")
        if not ok:
            self.skipTest("os_command_py not available in environment")

        import os_command_py.os_command as osc

        # Production MD command
        cmd = osc.Command(["/usr/local/bin/gmx", "mdrun", "-s", "prod_cuedc2.tpr", "-deffnm", "prod_cuedc2"])
        self.assertIn("-pin", cmd.cmd)
        self.assertIn("-pme", cmd.cmd)
        self.assertIn("-bonded", cmd.cmd)
        self.assertIn("-update", cmd.cmd)

    def test_hook_with_steep_mdp_file(self):
        ok = enable_gromacs_py_gpu_hook()
        if not ok:
            self.skipTest("os_command_py not available in environment")

        import os_command_py.os_command as osc

        with tempfile.TemporaryDirectory() as td:
            prev = os.getcwd()
            os.chdir(td)
            try:
                # Name doesn't contain 'em' or 'mini', but mdp has integrator = steep
                with open("custom_step.mdp", "w") as f:
                    f.write("integrator = steep\nnsteps = 1000\n")

                cmd = osc.Command(["/usr/local/bin/gmx", "mdrun", "-s", "custom_step.tpr", "-deffnm", "custom_step"])
                self.assertNotIn("-bonded", cmd.cmd)
                self.assertNotIn("-update", cmd.cmd)
                self.assertIn("-pin", cmd.cmd)
            finally:
                os.chdir(prev)


if __name__ == "__main__":
    unittest.main()

