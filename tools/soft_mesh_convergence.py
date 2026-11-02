#!/usr/bin/env python3
# tools/soft_mesh_convergence.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Convergenza della mesh GRADUATA dell'oggetto soffice (passo minimo soft_h_min al contatto).

Per ogni h_min e ogni tensione T (controller disattivato, T*dt sul cavo):
  F      forza di contatto letta dal solver (modules/contact_force.py, metodo validato)
  F_kin  forza dalla riduzione di freccia rispetto al moto libero (verifica indipendente)
  costo  ms/passo e stima della durata delle prove sul soffice nel Monte Carlo

Riferimento uniforme gia' misurato: griglia 5 / 9 / 13 -> 15.4 / 12.1 / 10.4 mN a T = 10 N.

Esecuzione (cartella SOFA_Project, source env.sh), ~30-60 min:
    python3.12 -u tools/soft_mesh_convergence.py
    python3.12 -u tools/soft_mesh_convergence.py --hmin 1.6 0.8 0.4 0.2 --levels 2 6 10
"""

import argparse
import math
import os
import sys
import time

import numpy as np

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

import Sofa.Core  # noqa: E402
import Sofa.Simulation  # noqa: E402

import main_scene  # noqa: E402
from modules.config import NOMINAL_PARAMS  # noqa: E402
from modules.contact_force import read_contact_force  # noqa: E402

MC_SOFT_TESTS = 42          # 21 casi x 2 configurazioni
STEPS_PER_TEST = 1300       # ~8 s di carico + 5 s di mantenimento a dt = 0.01


def build(**overrides):
    p = NOMINAL_PARAMS.copy()
    p.update(switched_mode=False, live_plot=False)
    p.update(overrides)
    root = Sofa.Core.Node("root")
    main_scene.createScene(root, p)
    ctrl = root.getObject("ActuationController")
    ctrl.update_actuation_physics = lambda _dt: None
    cable = root.getChild("NitinolTube").getChild("Tendon").getObject("tendonCable")
    cable.value.value = [0.0]
    Sofa.Simulation.init(root)
    return root, ctrl, cable


def centroid_x(root):
    """Centro del contatto lungo x (pesi = |lambda_y| sulla superficie dell'oggetto)."""
    mo = root.getChild("SoftContactObject").getChild("CollisionModel").getObject("dofs")
    lam = np.asarray(mo.getData("lambda").value, dtype=float)
    if lam.size == 0:
        return float("nan")
    w = np.abs(lam.reshape((-1, 3))[:, 1])
    pos = np.asarray(mo.position.value, dtype=float)
    return float(np.dot(w, pos[:, 0]) / w.sum()) if w.sum() > 0 else float("nan")


def sweep(root, ctrl, cable, levels, a, contact):
    dt = float(root.dt.value)
    n_ramp = max(1, int(round(a.ramp / dt)))
    n_tot = n_ramp + max(1, int(round(a.hold / dt)))
    n_avg = max(1, int(round(a.avg / dt)))
    y0 = float(ctrl.get_tip_position()[1])
    rows, T_prev = [], 0.0
    for T in levels:
        dY, F, xc = [], [], []
        t0 = time.perf_counter()
        for i in range(n_tot):
            Ti = T_prev + (T - T_prev) * min(1.0, (i + 1) / n_ramp)
            cable.value.value = [float(Ti) * dt]
            Sofa.Simulation.animate(root, dt)
            if i >= n_tot - n_avg:
                dY.append(float(ctrl.get_tip_position()[1]) - y0)
                if contact:
                    F.append(read_contact_force(root))
                    xc.append(centroid_x(root))
        rows.append({"T": T, "dY": float(np.mean(dY)),
                     "F": float(np.mean(F)) if F else 0.0,
                     "F_ptp": float(np.ptp(F)) if F else 0.0,
                     "xc": float(np.nanmean(xc)) if xc else float("nan"),
                     "ms": (time.perf_counter() - t0) / n_tot * 1e3})
        T_prev = T
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hmin", type=float, nargs="+", default=[1.6, 0.8, 0.4, 0.2])
    ap.add_argument("--levels", type=float, nargs="+", default=[2.0, 6.0, 10.0])
    ap.add_argument("--growth", type=float, default=NOMINAL_PARAMS["soft_growth"])
    ap.add_argument("--hmax", type=float, default=NOMINAL_PARAMS["soft_h_max"])
    ap.add_argument("--ramp", type=float, default=0.3)
    ap.add_argument("--hold", type=float, default=1.0)
    ap.add_argument("--avg", type=float, default=0.3)
    a = ap.parse_args()

    P = NOMINAL_PARAMS
    L, Ls = P["length_flex"], P["length_shaft"]
    r_o = P["od"] / 2.0
    r_i = r_o - P["thickness"]
    EI = P["E_load"] * math.pi / 4.0 * (r_o ** 4 - r_i ** 4)
    c_beam = P["yt_offset"] * L ** 2 / (2.0 * EI)

    print("\n[0] Riferimento in moto libero")
    root, ctrl, cable = build(contact_object="none")
    free = {r["T"]: r["dY"] for r in sweep(root, ctrl, cable, a.levels, a, contact=False)}
    del root, ctrl, cable
    c = float(np.dot(list(free), list(free.values())) / np.dot(list(free), list(free)))
    EI_eff = EI * c_beam / c
    print(f"    c = {c:.5f} mm/N")

    summary = []
    for h in sorted(a.hmin, reverse=True):
        root, ctrl, cable = build(contact_object="soft", soft_mesh="graded", soft_h_min=h,
                                  soft_h_max=a.hmax, soft_growth=a.growth)
        n_nodes = len(root.getChild("SoftContactObject").getObject("state").position.value)
        print(f"\n[soffice graduato] h_min = {h:g} mm, {n_nodes} nodi")
        rows = sweep(root, ctrl, cable, a.levels, a, contact=True)
        del root, ctrl, cable
        print("     T |  F [mN] | F_kin [mN] | F/F_kin | p-p [mN] | x_c [mm] | ms/passo")
        for r in rows:
            arm = (r["xc"] if np.isfinite(r["xc"]) else Ls + L) - Ls
            F_kin = max(0.0, (free[r["T"]] - r["dY"]) * 6 * EI_eff / (arm ** 2 * (3 * L - arm)))
            rat = r["F"] / F_kin if F_kin > 1e-6 else float("nan")
            print(f"  {r['T']:4.1f} | {r['F'] * 1e3:7.2f} | {F_kin * 1e3:10.2f} | {rat:7.3f} | "
                  f"{r['F_ptp'] * 1e3:8.3f} | {r['xc']:8.2f} | {r['ms']:8.0f}")
        ms = float(np.mean([r["ms"] for r in rows]))
        summary.append((h, n_nodes, rows[-1]["F"], ms))

    print("\n=== Riepilogo (forza all'ultima tensione) ===")
    print("   h_min | nodi | F [mN] | variazione | ms/passo | stima 42 prove soffice [h]")
    prev = None
    for h, n, F, ms in summary:
        var = "" if prev is None else f"{(F - prev) / prev:+.1%}"
        hours = ms * 1e-3 * STEPS_PER_TEST * MC_SOFT_TESTS / 3600
        print(f"  {h:6.2f} | {n:4d} | {F * 1e3:6.2f} | {var:>10} | {ms:8.0f} | {hours:8.1f}")
        prev = F
    print("  Riferimento uniforme: 5/9/13 nodi per lato -> 15.4 / 12.1 / 10.4 mN a T = 10 N")


if __name__ == "__main__":
    main()
