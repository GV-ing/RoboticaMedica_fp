import tempfile
import unittest
from pathlib import Path

import numpy as np

from modules.config import NOMINAL_PARAMS, UNCERTAIN_GROUPS
from modules.contact_force import read_contact_force
from modules.geometry_gen import check_tet_mesh, create_hollow_tube_meshes
from modules.monte_carlo import analyze, build_tasks, regression_tornado, sample_realizations


class MonteCarloTests(unittest.TestCase):
    def test_grouped_perturbations_reproducible_and_bounded(self):
        first = sample_realizations(seed=17, count=20)
        self.assertEqual(first, sample_realizations(seed=17, count=20))
        self.assertEqual(len(first), 20)
        for r in first:
            self.assertEqual(set(r["deltas"]), set(UNCERTAIN_GROUPS))
            for group, keys in UNCERTAIN_GROUPS.items():
                d = r["deltas"][group]
                self.assertGreaterEqual(abs(d), 0.10)
                self.assertLessEqual(abs(d), 0.20)
                for key in keys:  # stesso fattore per tutte le chiavi del gruppo
                    self.assertAlmostEqual(r["overrides"][key], NOMINAL_PARAMS[key] * (1 + d))
            ratio = r["overrides"]["E_load"] / r["overrides"]["E_unload"]
            self.assertAlmostEqual(ratio, NOMINAL_PARAMS["E_load"] / NOMINAL_PARAMS["E_unload"])

    def test_task_list(self):
        tasks = build_tasks(sample_realizations(count=20), "/tmp/series")
        self.assertEqual(len(tasks), 21 * 2 * 3 + 1)
        self.assertEqual(len({t["case_id"] for t in tasks}), len(tasks))
        self.assertEqual(tasks[0]["test"], "soft")  # i casi piu' lunghi per primi

    def test_tornado_recovers_linear_model(self):
        reals = sample_realizations(seed=3, count=20)
        D = [r["deltas"] for r in reals]
        H = [15 + 4 * d["yt"] / 0.2 - 2 * d["E"] / 0.2 + 0.5 * d["k"] / 0.2 for d in D]
        t = regression_tornado(D, H)
        coef = {c["parameter"]: c["coef_pct_per_20pct"] for c in t["coefficients"]}
        self.assertAlmostEqual(coef["yt"], 4.0, places=6)
        self.assertAlmostEqual(coef["E"], -2.0, places=6)
        self.assertEqual([c["parameter"] for c in t["coefficients"][:2]], ["yt", "E"])
        self.assertAlmostEqual(t["r2"], 1.0, places=9)

    def test_analyze_statistics_and_target(self):
        reals = sample_realizations(seed=5, count=20)
        results = [{"case_id": "x", "realization": "nominal", "config": "direct",
                    "test": "free", "metrics": {"H_pct": 30.0}}]
        for cfg, base in (("fixed", 14.0), ("switched", 16.0)):
            results.append({"case_id": "x", "realization": "nominal", "config": cfg,
                            "test": "free", "metrics": {"H_pct": base}})
            for r in reals:
                H = base * (1 + 0.5 * r["deltas"]["yt"])
                results.append({"case_id": "x", "realization": r["realization"], "config": cfg,
                                "test": "free", "metrics": {"H_pct": H}})
                for obj in ("rigid", "soft"):
                    stable = not (cfg == "switched" and obj == "rigid" and r["realization"] == "r03")
                    results.append({"case_id": "x", "realization": r["realization"],
                                    "config": cfg, "test": obj,
                                    "metrics": {"stable": stable, "in_range_10_50mN": True}})
        s = analyze(results, reals)
        f, w = s["configs"]["fixed"], s["configs"]["switched"]
        self.assertEqual(f["n_realizations"], 20)
        self.assertAlmostEqual(f["robustness_index"], f["H_worst_pct"] / 14.0)
        self.assertLessEqual(f["robustness_index"], 1.3)
        self.assertTrue(f["meets_target"])
        self.assertEqual(w["unstable_total"], 1)
        self.assertFalse(w["meets_target"])
        self.assertEqual(s["most_robust_config"], "fixed")
        self.assertEqual(f["dominant_parameters"][0], "yt")
        self.assertAlmostEqual(f["H_nominal_over_direct"], 14.0 / 30.0)


class MeshTests(unittest.TestCase):
    def test_tube_surface_is_closed_and_tetrahedra_are_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            obj_path = Path(directory) / "tube.obj"
            nodes, tetrahedra = create_hollow_tube_meshes(str(obj_path))
            bad_count, min_volume, total_volume = check_tet_mesh(nodes, tetrahedra)
            self.assertEqual(bad_count, 0)
            self.assertGreater(min_volume, 0.0)
            self.assertGreater(total_volume, 0.0)

            faces = [
                [int(index) for index in line.split()[1:]]
                for line in obj_path.read_text(encoding="utf-8").splitlines()
                if line.startswith("f ")
            ]
            edge_counts = {}
            for face in faces:
                for start, end in zip(face, face[1:] + face[:1]):
                    edge = tuple(sorted((start, end)))
                    edge_counts[edge] = edge_counts.get(edge, 0) + 1
            self.assertEqual(set(edge_counts.values()), {2})


class GradedMeshTests(unittest.TestCase):
    def test_graded_mesh_is_valid_and_surfaces_lie_on_tube(self):
        from modules.geometry_gen import graded_x_vals, write_surface_obj

        x_vals = graded_x_vals(100.0, 40.0, h_flex=0.5, n_shaft=4)
        self.assertAlmostEqual(x_vals[0], 0.0)
        self.assertAlmostEqual(x_vals[-1], 100.0)
        self.assertIn(40.0, np.round(x_vals, 9))
        self.assertTrue(np.all(np.diff(x_vals[5:]) <= 0.5 + 1e-9))

        with tempfile.TemporaryDirectory() as directory:
            nodes, tetrahedra = create_hollow_tube_meshes(
                str(Path(directory) / "fem.obj"), x_vals=x_vals)
            bad, _, volume = check_tet_mesh(nodes, tetrahedra)
            self.assertEqual(bad, 0)
            r_o, r_i = 1.93 / 2, 1.93 / 2 - 0.22
            # volume del prisma poligonale a 12 lati ~ volume del tubo (entro 5%)
            self.assertAlmostEqual(volume / (np.pi * (r_o**2 - r_i**2) * 100.0), 1.0, delta=0.05)

            obj = Path(directory) / "coll.obj"
            n_vertices = write_surface_obj(str(obj), 1.93, 0.22, 75.0, 100.0, 1.0)
            verts = np.array([[float(v) for v in line.split()[1:]]
                              for line in obj.read_text().splitlines()
                              if line.startswith("v ")])
            self.assertEqual(len(verts), n_vertices)
            self.assertGreaterEqual(verts[:, 0].min(), 75.0 - 1e-9)
            radii = np.hypot(verts[:, 1], verts[:, 2])
            self.assertTrue(np.all(np.isclose(radii, r_o) | np.isclose(radii, r_i)))


class GradedBoxMeshTests(unittest.TestCase):
    def test_graded_axis_and_box_mesh(self):
        from modules.geometry_gen import box_hexa_mesh, box_surface, graded_axis

        zs = graded_axis(-10.0, 10.0, 0.0, 0.4, 5.0, 1.4)
        self.assertAlmostEqual(zs[0], -10.0)
        self.assertAlmostEqual(zs[-1], 10.0)
        self.assertIn(0.0, zs)
        d = np.diff(zs)
        self.assertTrue(np.all(d > 0))
        self.assertAlmostEqual(d.min(), 0.4, places=9)
        self.assertLessEqual(d.max(), 5.0 + 1e-9)

        xs = graded_axis(80.0, 105.0, 100.0, 0.4, 5.0)
        ys = graded_axis(1.0, 21.0, 1.0, 0.4, 5.0)
        pos, hexa = box_hexa_mesh(xs, ys, zs)
        p, h = np.array(pos), np.array(hexa)
        self.assertEqual(len(p), len(xs) * len(ys) * len(zs))
        vol = np.einsum("ij,ij->i", np.cross(p[h[:, 1]] - p[h[:, 0]], p[h[:, 3]] - p[h[:, 0]]),
                        p[h[:, 4]] - p[h[:, 0]])
        self.assertTrue(np.all(vol > 0))                       # orientamento positivo
        self.assertAlmostEqual(vol.sum(), 25.0 * 20.0 * 20.0)  # volume del blocco

        sp, st = box_surface(xs, ys, zs)
        sp = np.array(sp)
        on_face = (np.isclose(sp[:, 0], 80) | np.isclose(sp[:, 0], 105) |
                   np.isclose(sp[:, 1], 1) | np.isclose(sp[:, 1], 21) |
                   np.isclose(sp[:, 2], -10) | np.isclose(sp[:, 2], 10))
        self.assertTrue(np.all(on_face))
        edges = {}
        for t in st:
            for e in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
                key = tuple(sorted(e))
                edges[key] = edges.get(key, 0) + 1
        self.assertEqual(set(edges.values()), {2})             # superficie chiusa


class ContactForceProviderTests(unittest.TestCase):
    """read_contact_force legge constraintForces del solver ed esclude la riga del cavo."""

    @staticmethod
    def _fake_root(forces, cable_impulse, dt=0.01):
        class Data:
            def __init__(self, value):
                self.value = value

        class Obj:
            def __init__(self, **data):
                for k, v in data.items():
                    setattr(self, k, Data(v))

        class Node:
            def __init__(self, children=None, objects=None):
                self.children = children or {}
                self.objects = objects or {}

            def getChild(self, name):
                return self.children.get(name)

            def getObject(self, name):
                return self.objects.get(name)

        tendon = Node(objects={"tendonCable": Obj(value=[cable_impulse])})
        root = Node(
            children={"NitinolTube": Node(children={"Tendon": tendon})},
            objects={"ConstraintSolver": Obj(constraintForces=forces)},
        )
        root.dt = Data(dt)
        return root

    def test_sums_normal_rows_and_excludes_cable(self):
        dt = 0.01
        cable = 4.0 * dt  # T = 4 N imposta come impulso
        # due contatti (n, t1, t2) + cavo in ultima riga
        forces = [0.0003, 1e-5, -2e-5, 0.0002, 3e-5, 0.0, cable]
        root = self._fake_root(forces, cable, dt)
        self.assertAlmostEqual(read_contact_force(root), 0.05)  # (3e-4 + 2e-4) / 0.01

    def test_cable_row_found_by_value_anywhere(self):
        dt = 0.01
        cable = 2.0 * dt
        forces = [cable, 0.0001, 0.0, 0.0]
        self.assertAlmostEqual(read_contact_force(self._fake_root(forces, cable, dt)), 0.01)

    def test_no_contact_returns_zero_and_bad_rows_raise(self):
        dt = 0.01
        self.assertEqual(read_contact_force(self._fake_root([0.03], 0.03, dt)), 0.0)
        with self.assertRaises(RuntimeError):
            read_contact_force(self._fake_root([0.001, 0.0, 0.03], 0.03, dt))


if __name__ == "__main__":
    unittest.main()
