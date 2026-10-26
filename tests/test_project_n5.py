import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path

import numpy as np

from modules.config import NOMINAL_PARAMS, UNCERTAIN_PARAMS
from modules.contact_force import read_contact_force
from modules.geometry_gen import check_tet_mesh, create_hollow_tube_meshes
from modules.monte_carlo import (
    contact_metrics,
    hysteresis_area_Nmm,
    run_monte_carlo,
    sample_realizations,
    summarize_hysteresis,
    tornado_coefficients,
)


class MonteCarloTests(unittest.TestCase):
    def test_sampled_perturbations_are_reproducible_and_bounded(self):
        first = sample_realizations(seed=17, count=20)
        second = sample_realizations(seed=17, count=20)
        self.assertEqual(first, second)
        self.assertEqual(len(first), 20)

        for realization in first:
            self.assertGreater(
                realization["params"]["E_load"],
                realization["params"]["E_unload"],
            )
            for key in UNCERTAIN_PARAMS:
                delta = realization["relative_perturbations"][key]
                self.assertGreaterEqual(abs(delta), 0.10)
                self.assertLessEqual(abs(delta), 0.20)
                self.assertAlmostEqual(
                    realization["params"][key],
                    NOMINAL_PARAMS[key] * (1.0 + delta),
                )

    def test_metric_helpers(self):
        self.assertAlmostEqual(
            hysteresis_area_Nmm([0.0, 1.0, 0.0], [1.0, 1.0, 0.0]),
            0.5,
        )
        contact = contact_metrics([0.0, 0.020, 0.021], 0.005, 0.020, 0.005)
        self.assertAlmostEqual(contact["contact_fraction"], 2.0 / 3.0)
        self.assertTrue(contact["stable"])
        summary = summarize_hysteresis([1.0, 2.0], 0.5)
        self.assertEqual(summary["worst_Nmm"], 2.0)
        self.assertEqual(summary["robustness_index"], 4.0)

        perturbations = [
            {"p": -0.1, "q": 0.1},
            {"p": 0.1, "q": -0.1},
            {"p": -0.2, "q": -0.2},
            {"p": 0.2, "q": 0.2},
        ]
        tornado = tornado_coefficients(perturbations, [0.8, 1.2, 0.6, 1.4])
        self.assertEqual({row["parameter"] for row in tornado}, {"p", "q"})

    def test_monte_carlo_runs_all_contact_and_spring_scenarios(self):
        def simulate(params, steps):
            phase = np.linspace(0.0, 2.0 * np.pi, 101)
            gain = (
                params["E_load"] / NOMINAL_PARAMS["E_load"]
                + params["k_spring"] / NOMINAL_PARAMS["k_spring"]
            )
            return {
                "tip_displacement_mm": np.sin(phase),
                "tendon_force_N": 2.0 + gain * np.cos(phase),
                "contact_force_N": np.full(
                    steps, params["contact_target_force"]
                ),
            }

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "results.json"
            progress = io.StringIO()
            with redirect_stderr(progress):
                results = run_monte_carlo(
                    simulate,
                    seed=17,
                    count=20,
                    output_path=output,
                )
            self.assertEqual(len(results), 4)
            self.assertTrue(all(
                len(result["realizations"]) == 20
                for result in results.values()
            ))
            self.assertEqual(set(json.loads(output.read_text(encoding="utf-8"))), set(results))
            self.assertIn("0/84", progress.getvalue())
            self.assertIn("84/84", progress.getvalue())


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


class ContactForceProviderTests(unittest.TestCase):
    def test_reads_resultant_normal_reaction_from_contact_surface(self):
        class Data:
            def __init__(self, value):
                self.value = value

        class MechanicalObject:
            def getData(self, name):
                self.asserted_name = name
                return Data([[1.0, 0.003, 2.0], [-1.0, -0.001, 1.0]])

        class Node:
            def __init__(self, children=None, objects=None):
                self.children = children or {}
                self.objects = objects or {}

            def getChild(self, name):
                return self.children.get(name)

            def getObject(self, name):
                return self.objects.get(name)

        dofs = MechanicalObject()
        contact = Node(children={
            "CollisionModel": Node(objects={"dofs": dofs}),
        })
        root = Node(children={"SoftContactObject": contact})
        self.assertAlmostEqual(read_contact_force(root), 0.002)
        self.assertEqual(dofs.asserted_name, "lambda")


if __name__ == "__main__":
    unittest.main()
