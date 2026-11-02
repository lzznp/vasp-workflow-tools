"""Synthetic regression tests; no project or licensed VASP files are used."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import vasp_audit as va

POSCAR = """Synthetic three atom cell
1
10 0 0
0 10 0
0 0 10
Si H
1 2
Selective dynamics
Direct
0 0 0 F F F
0.1 0 0 T T T
0.2 0 0 F F F
"""
SETTINGS = {key: "0" for key in va.KEYS}
SETTINGS.update(ENCUT="400", PREC="ACCURATE", ISMEAR="0", SIGMA="0.05",
                LREAL=".FALSE.", ISPIN="2", EDIFF="1E-6", NSW="0", IBRION="-1",
                NUPDOWN="2")

def block(energy=-10.0, status="converged", shift=0.0, mobile_force=0.01):
    flag = ("---------------- aborting loop because EDIFF is reached ----------------"
            if status == "converged" else
            "---------------- aborting loop EDIFF was not reached (unconverged) ----------------")
    return f"""{flag}
POSITION                                       TOTAL-FORCE (eV/Angst)
-----------------------------------------------------------------------------------
0 0 0 99 0 0
{1+shift:.8f} 0 0 {mobile_force} 0 0
2 0 0 0 0 0
 free energy TOTEN = -99999
 free  energy   TOTEN = {energy} eV
"""

def make_run(path, text=None, settings=None, poscar=POSCAR):
    path.mkdir()
    (path/"POSCAR").write_text(poscar)
    values = dict(SETTINGS); values.update(settings or {})
    (path/"INCAR").write_text("\n".join(f"{k} = {v}" for k, v in values.items()))
    (path/"OUTCAR").write_text(text or block()+"General timing and accounting\n")
    (path/"KPOINTS").write_text("synthetic gamma\n0\nGamma\n1 1 1\n0 0 0\n")
    # Deliberately invalid placeholder, used ONLY to test equality checks.
    (path/"POTCAR").write_text("SYNTHETIC HASH PLACEHOLDER, NOT PSEUDOPOTENTIAL DATA")
    return path

class RegressionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
    def tearDown(self):
        self.temp.cleanup()
    def test_frozen_atoms_do_not_pollute_force_maximum(self):
        run = make_run(self.root/"a")
        result = va.audit(run)
        self.assertAlmostEqual(result["force_metric"]["fmax_projected_eVA"], 0.01)
        self.assertAlmostEqual(result["force_metric"]["fmax_all_atoms_eVA"], 99)
        self.assertTrue(result["selected_subspace_below_threshold"])
    def test_exit_zero_does_not_establish_scf(self):
        run = make_run(self.root/"a", block(status="unconverged")+"General timing and accounting\n")
        (run/"EXIT_CODE").write_text("0")
        result = va.audit(run)
        self.assertEqual(result["last_frame_scf"], "unconverged")
        self.assertIsNone(result["selected_subspace_below_threshold"])
        self.assertFalse(result["geometry_convergence_message"])
    def test_static_scf_does_not_require_geometry_message(self):
        result = va.audit(make_run(self.root/"a"))
        self.assertEqual(result["kind"], "static")
        self.assertEqual(result["scf_converged_frames"], 1)
        self.assertFalse(result["geometry_convergence_message"])
    def test_unfinished_output_is_not_accepted_as_complete_frame(self):
        text = "POSITION TOTAL-FORCE (eV/Angst)\n---\n0 0 0 1 1 1\n"
        result = va.audit(make_run(self.root/"a", text))
        self.assertEqual(result["scf_complete_force_frames"], 0)
        self.assertIsNone(result["force_metric"])
        self.assertIn("Incomplete force block ignored", result["warnings"])
    def test_exact_double_space_toten_ignores_electronic_iteration_energy(self):
        result = va.audit(make_run(self.root/"a"))
        self.assertEqual(result["force_metric"]["energy_eV"], -10)
    def test_phva_uses_center_and_excludes_imaginary_mode_from_zpe(self):
        modes = """ 1 f  = 10 THz 60 2PiTHz 1000 cm-1 100 meV
 2 f  = 20 THz 120 2PiTHz 2000 cm-1 200 meV
 3 f/i= 3 THz 18 2PiTHz 300 cm-1 40 meV
"""
        text = block(-10)+block(-9, mobile_force=3)+modes+"General timing and accounting\n"
        result = va.audit(make_run(self.root/"a", text, {"IBRION": "5", "NSW": "1"}))
        self.assertEqual(result["force_metric"]["energy_eV"], -10)
        self.assertEqual(result["imaginary_count"], 1)
        self.assertEqual(result["frequencies_cm1"], [1000, 2000, -300])
        self.assertAlmostEqual(result["zpe_real_modes_eV"], 0.15)
        self.assertAlmostEqual(result["force_metric"]["fmax_projected_eVA"], 0.01)
    def test_periodic_distance_crossing_boundary(self):
        cell = np.eye(3)*10
        vector = va.minimum_image(np.array([9.8, 0, 0]), cell)
        self.assertAlmostEqual(np.linalg.norm(vector), 0.2)
    def test_skew_cell_requires_more_than_fractional_rounding(self):
        cell = np.array([[1,0,0],[2,.3,0],[0,0,5.]])
        v = np.array([.01,.147,0])
        got = va.minimum_image(v, cell)
        brute = min(np.linalg.norm(v-np.array([i,j,k])@cell)
                    for i in range(-5,6) for j in range(-5,6) for k in range(-1,2))
        self.assertAlmostEqual(np.linalg.norm(got), brute)
        rounded = (v@np.linalg.inv(cell)-np.rint(v@np.linalg.inv(cell)))@cell
        self.assertLess(np.linalg.norm(got), np.linalg.norm(rounded))
    def test_mixed_flags_follow_direct_lattice_directions(self):
        cell = np.array([[1,1,0],[0,1,0],[0,0,1.]])
        result = va.projected_forces([[0,1,0]], np.array([[True,False,False]]), cell)
        np.testing.assert_allclose(result, [[.5,.5,0]])
    def test_negative_volume_and_cartesian_scale(self):
        p = self.root/"POSCAR"
        p.write_text("demo\n-8\n1 0 0\n0 1 0\n0 0 1\nH\n1\nCartesian\n.5 .5 .5\n")
        s = va.read_poscar(p)
        self.assertAlmostEqual(np.linalg.det(s.cell), 8)
        np.testing.assert_allclose(s.positions, [[1,1,1]])
    def test_three_component_scale(self):
        p = self.root/"POSCAR"
        p.write_text("demo\n2 3 4\n1 0 0\n0 1 0\n0 0 1\nH\n1\nCartesian\n1 1 1\n")
        np.testing.assert_allclose(va.read_poscar(p).positions, [[2,3,4]])
    def test_truncated_poscar_rejected(self):
        p=self.root/"POSCAR";p.write_text(POSCAR.rsplit("\n",2)[0]+"\n")
        with self.assertRaises(ValueError): va.read_poscar(p)
    def test_compare_uses_evaluated_frame_not_proposed_contcar(self):
        left=make_run(self.root/"a", block(-10)+block(-11,shift=.1))
        right=make_run(self.root/"b", block(-8,shift=.1), {"NUPDOWN":"4"})
        (left/"CONTCAR").write_text("this unrelated proposed geometry is not used")
        result=va.compare(left,right)
        self.assertTrue(result["comparable"])
        self.assertAlmostEqual(result["energy_difference_eV"],3)
    def test_geometry_mismatch_blocks_energy_subtraction(self):
        left=make_run(self.root/"a");right=make_run(self.root/"b",block(-8,shift=.5))
        result=va.compare(left,right)
        self.assertFalse(result["comparable"])
        self.assertIsNone(result["energy_difference_eV"])
    def test_reference_mismatch_blocks_comparison(self):
        left=make_run(self.root/"a");right=make_run(self.root/"b")
        (right/"POTCAR").write_text("ANOTHER SYNTHETIC PLACEHOLDER")
        self.assertFalse(va.compare(left,right)["comparable"])
    def test_missing_numerical_metadata_is_unknown(self):
        left=make_run(self.root/"a");right=make_run(self.root/"b")
        (right/"INCAR").write_text("ENCUT=400")
        self.assertFalse(va.compare(left,right)["comparable"])
    def test_nul_corruption_blocks_threshold_acceptance(self):
        run=make_run(self.root/"a",block()+"NUL marker: \0\nGeneral timing and accounting\n")
        result=va.audit(run)
        self.assertIsNone(result["selected_subspace_below_threshold"])
        self.assertTrue(any("NUL" in w for w in result["warnings"]))
    def test_frame_uses_updated_outcar_lattice(self):
        lattice=" direct lattice vectors                 reciprocal lattice vectors\n20 0 0 0 0 0\n0 20 0 0 0 0\n0 0 20 0 0 0\n"
        run=make_run(self.root/"a",lattice+block()+"General timing and accounting\n")
        structure=va.read_poscar(run/"POSCAR")
        output=va.read_outcar(run/"OUTCAR",3,structure.cell)
        np.testing.assert_allclose(output["frames"][0]["cell"],np.eye(3)*20)
    def test_process_inventory_reports_overlapping_affinities(self):
        proc=self.root/"proc";proc.mkdir()
        for pid,cwd in [(101,"job-a"),(102,"job-b")]:
            d=proc/str(pid);d.mkdir()
            (d/"exe").write_text("synthetic proc executable")
            (d/"cwd").write_text(cwd)
            (d/"status").write_text("Name: vasp\nCpus_allowed_list:\t0-3\n")
        with patch.object(va.Path, "resolve", return_value=Path("/fake/vasp_std")), patch.object(va.os, "readlink", side_effect=lambda p: "/fake/job-"+p.parent.name):
            result=va.processes(proc)
        self.assertEqual(len(result["jobs"]),2)
        self.assertEqual(result["affinity_overlaps"][0]["cpus"],[0,1,2,3])

if __name__ == "__main__":
    unittest.main()
