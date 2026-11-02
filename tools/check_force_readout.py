#!/usr/bin/env python3
# tools/check_force_readout.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Confronto di metodi di lettura della forza di contatto, per scegliere quello affidabile.

Contesto (tools/check_soft_contact.py): sul soffice il lambda letto sul modello di
collisione MAPPATO vale esattamente 2x oppure 1x la forza ricavata dalla cinematica
(doppio conteggio intermittente); sul rigido (modello non mappato) la lettura e' esatta.

Metodi confrontati a ogni livello di T (controller disattivato, T*dt sul cavo):
  A  coll_obj   lambda sul MechanicalObject di collisione dell'oggetto (metodo attuale)
  B  fem_obj    lambda sul MechanicalObject FEM dell'oggetto soffice ('state', non mappato)
  C  solver     constraintForces del constraint solver: somma delle righe normali di
                contatto (FrictionContact: blocchi di 3 righe n,t1,t2), escludendo la riga
                del cavo; riportate due ipotesi sulla posizione della riga del cavo
Riferimenti:
  rigido   F_teo = 1.5*(yt/L)*(T - T_c)        (equilibrio, gia' verificato entro 0.5%)
  soffice  F_kin dalla riduzione di freccia con carico in x_c (cinematica, senza lambda)

Esecuzione (cartella SOFA_Project, source env.sh), ~15-20 min:
    python3.12 -u tools/check_force_readout.py
    python3.12 -u tools/check_force_readout.py --objects soft --grid 9
"""

import argparse
import math
import os
import sys

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


def lam_y(mo, dt):
    """Risultante Y [N] e centroide x [mm] delle reazioni lambda su un MechanicalObject."""
    if mo is None:
        return float("nan"), float("nan")
    try:
        lam = np.asarray(mo.getData("lambda").value, dtype=float)
    except Exception:
        return float("nan"), float("nan")
    if lam.size == 0 or lam.size % 3:
        return 0.0, float("nan")
    f = lam.reshape((-1, 3)) / dt
    w = np.abs(f[:, 1])
    pos = np.asarray(mo.position.value, dtype=float)
    xc = float(np.dot(w, pos[:, 0]) / w.sum()) if w.sum() > 0 else float("nan")
    return float(abs(f[:, 1].sum())), xc


def solver_normals(solver, dt, T):
    """Somma delle righe normali di contatto da constraintForces, due ipotesi sul cavo."""
    try:
        cf = np.asarray(solver.constraintForces.value, dtype=float).ravel()
    except Exception:
        return float("nan"), float("nan"), 0
    n = cf.size
    if n == 0:
        return 0.0, 0.0, 0
    out = []
    for cable_row in (0, n - 1):           # cavo come prima o come ultima riga
        rest = np.delete(cf, cable_row)
        if rest.size % 3:
            out.append(float("nan"))
            continue
        out.append(float(abs(rest[0::3].sum()) / dt))
    cable_check = cf[0] / dt, cf[-1] / dt  # una delle due deve valere ~T
    return out[0], out[1], n, cable_check


def run(obj, levels, a, grid):
    root, ctrl, cable = build(contact_object=obj, obj_grid_n=grid)
    dt = float(root.dt.value)
    node = root.getChild("SoftContactObject" if obj == "soft" else "RigidContactObject")
    coll_mo = node.getChild("CollisionModel").getObject("dofs")
    fem_mo = node.getObject("state") if obj == "soft" else None
    solver = root.getObject("ConstraintSolver")
    n_ramp = max(1, int(round(a.ramp / dt)))
    n_tot = n_ramp + max(1, int(round(a.hold / dt)))
    n_avg = max(1, int(round(a.avg / dt)))
    y0 = float(ctrl.get_tip_position()[1])
    rows, T_prev = [], 0.0
    for T in levels:
        acc = {k: [] for k in ("dY", "A", "B", "C0", "Cn", "xc", "nrows", "cab0", "cabn")}
        for i in range(n_tot):
            Ti = T_prev + (T - T_prev) * min(1.0, (i + 1) / n_ramp)
            cable.value.value = [float(Ti) * dt]
            Sofa.Simulation.animate(root, dt)
            if i >= n_tot - n_avg:
                acc["dY"].append(float(ctrl.get_tip_position()[1]) - y0)
                fa, xc = lam_y(coll_mo, dt)
                fb, _ = lam_y(fem_mo, dt)
                res = solver_normals(solver, dt, T)
                acc["A"].append(fa)
                acc["xc"].append(xc)
                acc["B"].append(fb)
                if len(res) == 4:
                    c0, cn, nr, (k0, kn) = res
                    acc["C0"].append(c0)
                    acc["Cn"].append(cn)
                    acc["nrows"].append(nr)
                    acc["cab0"].append(k0)
                    acc["cabn"].append(kn)
        row = {"T": T}
        for k, v in acc.items():
            row[k] = float(np.nanmean(v)) if len(v) and not np.all(np.isnan(v)) else float("nan")
        rows.append(row)
        T_prev = T
    del root, ctrl, cable
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--objects", nargs="+", default=["rigid", "soft"], choices=["rigid", "soft"])
    ap.add_argument("--grid", type=int, default=NOMINAL_PARAMS["obj_grid_n"])
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
    dt = float(root.dt.value)
    y0 = float(ctrl.get_tip_position()[1])
    free, T_prev = {}, 0.0
    n_ramp = max(1, int(round(a.ramp / dt)))
    n_tot = n_ramp + max(1, int(round(a.hold / dt)))
    for T in a.levels:
        for i in range(n_tot):
            cable.value.value = [float(T_prev + (T - T_prev) * min(1.0, (i + 1) / n_ramp)) * dt]
            Sofa.Simulation.animate(root, dt)
        free[T] = float(ctrl.get_tip_position()[1]) - y0
        T_prev = T
    del root, ctrl, cable
    c = float(np.dot(list(free), list(free.values())) / np.dot(list(free), list(free)))
    EI_eff = EI * c_beam / c
    T_c = max(0.0, P["obj_gap"] - P["contact_distance"]) / c
    print(f"    c = {c:.5f} mm/N, inizio contatto atteso a T ~ {T_c:.2f} N")

    for obj in a.objects:
        print(f"\n[{obj}] obj_grid_n = {a.grid}" if obj == "soft" else f"\n[{obj}]")
        rows = run(obj, a.levels, a, a.grid)
        print("     T | rif [mN] | A coll_obj | B fem_obj | C cavo=0 | C cavo=fine | righe | "
              "T dal solver (0 / fine) [N]")
        for r in rows:
            if obj == "rigid":
                ref = max(0.0, 1.5 * yt / L * (r["T"] - T_c))
            else:
                arm = (r["xc"] if np.isfinite(r["xc"]) else Ls + L) - Ls
                ref = max(0.0, (free[r["T"]] - r["dY"]) * 6 * EI_eff / (arm ** 2 * (3 * L - arm)))
            print(f"  {r['T']:4.1f} | {ref * 1e3:8.2f} | {r['A'] * 1e3:10.2f} | "
                  f"{r['B'] * 1e3:9.2f} | {r['C0'] * 1e3:8.2f} | {r['Cn'] * 1e3:11.2f} | "
                  f"{r['nrows']:5.0f} | {r['cab0']:6.3f} / {r['cabn']:6.3f}")
        print("  rif = equilibrio (rigido) oppure cinematica F_kin (soffice). "
              "Il metodo affidabile coincide con rif a tutti i livelli.")


if __name__ == "__main__":
    main()
