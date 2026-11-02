import unittest

import numpy as np

from modules.config import NOMINAL_PARAMS
from modules.metrics import contact_stability, hysteresis_H
from modules.spring_model import SeriesSpring

DT = 0.01


def run_cycle(params, rate=1.0, u1_of_T=lambda T: 0.0, contact_of_T=lambda T: 0.0):
    """Ciclo 0 -> T_max -> 0 con il solo modello della molla (tubo idealmente rigido)."""
    s = SeriesSpring(params)
    u, T, out = 0.0, 0.0, []
    phase = "load"
    while True:
        r = s.max_motor_rate(rate) * DT
        if phase == "load":
            u += r
            if s.T >= s.T_max * (1 - 1e-6) or u >= params["u_max"]:
                phase = "unload"
        else:
            u = max(0.0, u - r)
        T = s.step(u, u1_of_T(T), contact_of_T(T), DT)
        out.append((u, T, s.k, s.z, phase))
        if phase == "unload" and u <= 0.0:
            return s, np.array(out, dtype=object)


class SpringModelTests(unittest.TestCase):
    def setUp(self):
        self.p = NOMINAL_PARAMS.copy()

    def test_dahl_state_saturates_and_reverses_continuously(self):
        s = SeriesSpring(self.p)
        for _ in range(100):
            s.update_dahl(0.01)          # 1 mm in carico: z -> ~1 (sigma0 = 10/mm)
        self.assertGreater(s.z, 0.9999)
        z0 = s.z
        s.update_dahl(-0.1)              # 0.1 mm in scarico: z = -1 + (1 + z0)*exp(-1)
        self.assertAlmostEqual(s.z, -1 + (1 + z0) * np.exp(-1), places=12)

    def test_tension_is_continuous_saturated_and_rate_limited(self):
        _, out = run_cycle(self.p)
        T = out[:, 1].astype(float)
        self.assertAlmostEqual(T.max(), self.p["T_max"], places=6)
        self.assertGreaterEqual(T.min(), 0.0)
        self.assertLessEqual(np.max(np.abs(np.diff(T))), self.p["max_force_rate"] * DT + 1e-12)

    def test_loop_width_matches_glossary_deadband(self):
        """Glossario: dY_deadband = 2*F_hyst,max / k = larghezza in u del ciclo T-u a pari T."""
        s, out = run_cycle(self.p)
        u, T = out[:, 0].astype(float), out[:, 1].astype(float)
        i_rev = int(np.argmax(u))
        T_level = 5.0  # a meta' ciclo, con z gia' saturato in entrambi i rami
        u_load = np.interp(T_level, T[:i_rev], u[:i_rev])
        u_unload = np.interp(T_level, T[i_rev:][::-1], u[i_rev:][::-1])
        expected = 2 * s.F_hyst_max / self.p["k_spring"]
        self.assertAlmostEqual(u_unload - u_load, expected, delta=0.03 * expected)

    def test_switched_stiffness_is_slew_limited(self):
        p = dict(self.p, switched_mode=True)
        s, out = run_cycle(p, contact_of_T=lambda T: 0.02 if T > 1.0 else 0.0)
        k = out[:, 2].astype(float)
        self.assertAlmostEqual(k.min(), p["k_low"])
        self.assertAlmostEqual(k.max(), p["k_high"])
        self.assertLessEqual(np.max(np.abs(np.diff(k))), p["dk_dt_max"] * DT + 1e-12)

    def test_slack_spring_transmits_no_force(self):
        s = SeriesSpring(self.p)
        self.assertEqual(s.step(0.5, 0.6, 0.0, DT), 0.0)  # u1 > u: molla scarica


class MetricsTests(unittest.TestCase):
    def test_H_on_synthetic_loop(self):
        u_up = np.linspace(0, 8, 200)
        u_dn = u_up[::-1]
        y_up = 0.03 * u_up
        y_dn = 0.03 * u_dn + 0.012 * np.sin(np.pi * u_dn / 8)  # scarto massimo 0.012 a u = 4
        u = np.concatenate([u_up, u_dn])
        y = np.concatenate([y_up, y_dn])
        ph = np.array(["load"] * 200 + ["unload"] * 200)
        r = hysteresis_H(u, y, ph)
        self.assertAlmostEqual(r["H_pct"], 100 * 0.012 / 0.24, places=2)
        self.assertAlmostEqual(r["u_at_gap_max_mm"], 4.0, delta=0.05)

    def test_H_zero_without_hysteresis(self):
        u = np.concatenate([np.linspace(0, 5, 50), np.linspace(5, 0, 50)])
        ph = np.array(["load"] * 50 + ["unload"] * 50)
        self.assertAlmostEqual(hysteresis_H(u, 0.1 * u, ph)["H_pct"], 0.0, places=9)

    def test_contact_stability_criteria(self):
        t = np.arange(0, 6, 0.01)
        ph = np.array(["hold"] * t.size)
        ok = contact_stability(t, 0.02 + 0.002 * np.sin(5 * t), ph)
        self.assertTrue(ok["stable"])
        self.assertTrue(ok["in_range_10_50mN"])
        bad = contact_stability(t, 0.02 + 0.006 * np.sin(5 * t), ph)  # p-p 60%
        self.assertFalse(bad["stable"])
        F = np.full(t.size, 0.02)
        F[-10] = 0.0                                                   # perdita di contatto
        self.assertFalse(contact_stability(t, F, ph)["stable"])
        short = contact_stability(t[:200], np.full(200, 0.02), ph[:200])  # solo 2 s
        self.assertFalse(short["stable"])


if __name__ == "__main__":
    unittest.main()
