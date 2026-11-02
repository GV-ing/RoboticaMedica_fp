#!/usr/bin/env python3
# tools/check_soft_contact.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Diagnosi del contatto con l'oggetto soffice.

Per ogni livello di tensione T (controller disattivato, T*dt sul cavo) misura:
  F_obj   reazione letta sul modello di collisione dell'OGGETTO (come contact_force.py)
  F_tube  reazione letta sul modello di collisione del TUBO (deve coincidere con F_obj)
  x_c     punto medio di contatto lungo x (media pesata con la forza) e estensione
  F_kin   forza stimata dalla riduzione di freccia rispetto al moto libero,
          con il carico applicato in x_c:  delta = F * a^2 (3L - a) / (6 EI_eff),
          a = x_c - length_shaft, EI_eff = EI_trave * c_beam / c_misurata
e ripete per diverse risoluzioni della mesh dell'oggetto soffice (obj_grid_n).

Interpretazione:
  F_obj != F_tube       -> lettura di lambda non affidabile sul modello mappato
  F_obj ~ F_tube ~ F_kin -> forza corretta; eventuali salti sono fisici/di mesh
  forte dipendenza da obj_grid_n -> mesh dell'oggetto soffice non convergente

Esecuzione (cartella SOFA_Project, source env.sh), ~20-30 min con due griglie:
    python3.12 -u tools/check_soft_contact.py
    python3.12 -u tools/check_soft_contact.py --grids 5 9 13
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


def reaction(mo, dt):
    """Risultante Y [N], punto medio x [mm] e estensione x [mm] delle reazioni su un MO."""
    lam = np.asarray(mo.getData("lambda").value, dtype=float)
    if lam.size == 0:
        return 0.0, float("nan"), float("nan")
    f = lam.reshape((-1, 3)) / dt
    pos = np.asarray(mo.position.value, dtype=float)
    w = np.abs(f[:, 1])
    if w.sum() <= 0:
        return 0.0, float("nan"), float("nan")
    xs = pos[w > 1e-3 * w.max(), 0]
    return float(abs(f[:, 1].sum())), float(np.dot(w, pos[:, 0]) / w.sum()), float(np.ptp(xs))


def run(root, ctrl, cable, levels, a, contact):
    dt = float(root.dt.value)
    n_ramp = max(1, int(round(a.ramp / dt)))
    n_tot = n_ramp + max(1, int(round(a.hold / dt)))
    n_avg = max(1, int(round(a.avg / dt)))
    y0 = float(ctrl.get_tip_position()[1])
    if contact:
        obj_mo = root.getChild("SoftContactObject").getChild("CollisionModel").getObject("dofs")
        tube_mo = root.getChild("NitinolTube").getChild("CollisionModel").getObject("dofs")
    rows, T_prev = [], 0.0
    for T in levels:
        acc = {"dY": [], "Fo": [], "Ft": [], "xc": [], "ext": []}
        t0 = time.perf_counter()
        for i in range(n_tot):
            Ti = T_prev + (T - T_prev) * min(1.0, (i + 1) / n_ramp)
            cable.value.value = [float(Ti) * dt]
            Sofa.Simulation.animate(root, dt)
            if i >= n_tot - n_avg:
                acc["dY"].append(float(ctrl.get_tip_position()[1]) - y0)
                if contact:
                    fo, xc, ext = reaction(obj_mo, dt)
                    ft, _, _ = reaction(tube_mo, dt)
                    acc["Fo"].append(fo)
                    acc["Ft"].append(ft)
                    acc["xc"].append(xc)
                    acc["ext"].append(ext)
        row = {"T": T, "t_step": (time.perf_counter() - t0) / n_tot}
        for k, v in acc.items():
            row[k] = float(np.nanmean(v)) if v else float("nan")
        rows.append(row)
        T_prev = T
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grids", type=int, nargs="+", default=[5, 9],
                    help="Nodi per lato della griglia FEM dell'oggetto soffice")
    ap.add_argument("--levels", type=float, nargs="+", default=[1, 2, 4, 6, 8, 10])
    ap.add_argument("--ramp", type=float, default=0.3)
    ap.add_argument("--hold", type=float, default=1.0)
    ap.add_argument("--avg", type=float, default=0.3)
    a = ap.parse_args()

    P = NOMINAL_PARAMS
    L, Ls, yt = P["length_flex"], P["length_shaft"], P["yt_offset"]
    r_o = P["od"] / 2.0
    r_i = r_o - P["thickness"]
    EI = P["E_load"] * math.pi / 4.0 * (r_o ** 4 - r_i ** 4)
    c_beam = yt * L ** 2 / (2.0 * EI)

    print("\n[0] Riferimento in moto libero")
    root, ctrl, cable = build(contact_object="none")
    free = {r["T"]: r["dY"] for r in run(root, ctrl, cable, a.levels, a, contact=False)}
    del root, ctrl, cable
    c = float(np.dot(list(free), list(free.values())) / np.dot(list(free), list(free)))
    EI_eff = EI * c_beam / c
    print(f"    c = {c:.5f} mm/N, EI_eff = {EI_eff:.4g} N*mm^2")

    for n in a.grids:
        print(f"\n[soffice] obj_grid_n = {n} (elementi da {P['obj_size'] / (n - 1):.2f} mm)")
        root, ctrl, cable = build(contact_object="soft", obj_grid_n=n)
        rows = run(root, ctrl, cable, a.levels, a, contact=True)
        del root, ctrl, cable
        print("     T |   dY [mm] | F_obj [mN] | F_tube [mN] | x_c [mm] | estens. [mm] | F_kin [mN] | ms/passo")
        for r in rows:
            arm = (r["xc"] if np.isfinite(r["xc"]) else Ls + L) - Ls
            k_inv = arm ** 2 * (3 * L - arm) / (6 * EI_eff)  # cedevolezza punta per forza in x_c
            F_kin = max(0.0, (free[r["T"]] - r["dY"]) / k_inv)
            print(f"  {r['T']:4.1f} | {r['dY']:9.4f} | {r['Fo'] * 1e3:10.2f} | {r['Ft'] * 1e3:11.2f} | "
                  f"{r['xc']:8.2f} | {r['ext']:12.2f} | {F_kin * 1e3:10.2f} | {r['t_step'] * 1e3:8.1f}")


if __name__ == "__main__":
    main()
