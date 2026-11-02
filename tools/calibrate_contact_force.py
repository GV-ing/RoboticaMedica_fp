#!/usr/bin/env python3
# tools/calibrate_contact_force.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Calibrazione di `force_scale`: confronta la reazione di contatto letta da SOFA
(lambda del MechanicalObject di collisione) con la tensione T imposta al tendine.

Due verifiche indipendenti:

  TEST A - dipendenza da dt (decisivo)
      Stessa scena (oggetto rigido), stessa T, due passi temporali diversi.
      Se il valore grezzo letto si dimezza dimezzando dt, e' un IMPULSO [N*s]
      -> force_scale = 1/dt. Se resta uguale, e' gia' una FORZA [N] -> force_scale = 1.

  TEST B - equilibrio meccanico (ordine di grandezza)
      Il tendine con offset yt applica alla punta un momento M = T*yt.
      Trave a sbalzo (L = parte flessibile) con momento in punta e appoggio in
      punta dopo aver chiuso il gap g:
          F_teorica = 1.5 * (yt / L) * (T - g / c)
      dove c = dY_libero / T e' la cedevolezza misurata in moto libero (stessa
      mesh, quindi include l'irrigidimento numerico dei tetraedri lineari).

La legge della molla del controller viene disattivata: T e' imposta direttamente
al CableConstraint con una rampa, poi mantenuta finche' la scena e' quasi statica.

Esecuzione (dalla cartella SOFA_Project, con le variabili d'ambiente del README):
    /usr/bin/python3 -u tools/calibrate_contact_force.py
    /usr/bin/python3 -u tools/calibrate_contact_force.py --levels 3 5 7 9 --hold 2.0
"""

import argparse
import csv
import math
import os
import sys

import numpy as np

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

import Sofa.Core  # noqa: E402
import Sofa.Simulation  # noqa: E402

from main_scene import createScene  # noqa: E402
from modules.config import NOMINAL_PARAMS  # noqa: E402
from modules.contact_force import read_contact_force  # noqa: E402


def build_scene(contact_object, dt):
    """Crea e inizializza la scena con la molla virtuale disattivata."""
    params = NOMINAL_PARAMS.copy()
    params.update(contact_object=contact_object, dt=dt,
                  switched_mode=False, live_plot=False)
    root = Sofa.Core.Node("root")
    createScene(root, params)
    Sofa.Simulation.init(root)

    ctrl = root.getObject("ActuationController")
    if ctrl is None:
        raise RuntimeError("ActuationController non trovato")
    # Disattiva legge molla e commutazione di E: T la impone questo script.
    ctrl.update_actuation_physics = lambda _dt: None

    cable = root.getChild("NitinolTube").getChild("Tendon").getObject("tendonCable")
    if cable is None:
        raise RuntimeError("CableConstraint 'tendonCable' non trovato")
    return root, ctrl, cable, params


def run_levels(root, ctrl, cable, levels, t_ramp, t_hold, t_avg):
    """Per ogni livello di T: rampa, mantenimento, media sull'ultima finestra."""
    dt = float(root.dt.value)
    n_ramp = max(1, int(round(t_ramp / dt)))
    n_hold = max(1, int(round(t_hold / dt)))
    n_avg = min(n_hold, max(1, int(round(t_avg / dt))))

    tip0_y = float(ctrl.get_tip_position()[1])
    rows = []
    T_prev = 0.0
    for T in levels:
        raw, dy = [], []
        n_tot = n_ramp + n_hold
        for i in range(n_tot):
            T_i = T_prev + (T - T_prev) * min(1.0, (i + 1) / n_ramp)
            cable.value.value = [float(T_i) * dt]  # impulso: forza T_i
            Sofa.Simulation.animate(root, dt)
            if i >= n_tot - n_avg:
                raw.append(read_contact_force(root))  # gia' in N (lambda/dt)
                dy.append(float(ctrl.get_tip_position()[1]) - tip0_y)
        raw = np.asarray(raw)
        rows.append({
            "T_N": float(T),
            "raw_mean": float(raw.mean()),
            "raw_ptp": float(np.ptp(raw)),
            "dY_mm": float(np.mean(dy)),
        })
        T_prev = T
        print(f"    dt={dt:g}  T={T:5.2f} N  dY={rows[-1]['dY_mm']:+.4f} mm  "
              f"F_contatto={rows[-1]['raw_mean'] * 1e3:.3f} mN  (p-p {rows[-1]['raw_ptp']:.2e})")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--levels", type=float, nargs="+", default=[2, 4, 6, 8, 10],
                    help="Livelli di tensione T [N]")
    ap.add_argument("--dt", type=float, default=NOMINAL_PARAMS["dt"])
    ap.add_argument("--ramp", type=float, default=0.5, help="Durata rampa [s]")
    ap.add_argument("--hold", type=float, default=1.5, help="Mantenimento [s]")
    ap.add_argument("--avg", type=float, default=0.5, help="Finestra di media [s]")
    ap.add_argument("--csv", default="calibration_force.csv")
    a = ap.parse_args()

    dt1, dt2 = a.dt, a.dt / 2.0
    p = NOMINAL_PARAMS
    L, yt = p["length_flex"], p["yt_offset"]
    g_geom = p["obj_gap"]
    g_eff = max(0.0, g_geom - p["contact_distance"])  # il contatto scatta a contactDistance

    # Rigidezza analitica (solo per confronto con la mesh)
    r_o = p["od"] / 2.0
    r_i = r_o - p["thickness"]
    I = math.pi / 4.0 * (r_o ** 4 - r_i ** 4)
    EI = p["E_load"] * I
    c_analytic = yt * L ** 2 / (2.0 * EI)

    print("\n[1/3] Moto libero (nessun oggetto), dt =", dt1)
    root, ctrl, cable, _ = build_scene("none", dt1)
    free = run_levels(root, ctrl, cable, a.levels, a.ramp, a.hold, a.avg)
    del root, ctrl, cable

    T_arr = np.array([r["T_N"] for r in free])
    dY_arr = np.array([r["dY_mm"] for r in free])
    c_meas = float(np.dot(T_arr, dY_arr) / np.dot(T_arr, T_arr))  # retta per l'origine

    print("\n[2/3] Oggetto rigido, dt =", dt1)
    root, ctrl, cable, _ = build_scene("rigid", dt1)
    rig1 = run_levels(root, ctrl, cable, a.levels, a.ramp, a.hold, a.avg)
    del root, ctrl, cable

    print("\n[3/3] Oggetto rigido, dt =", dt2)
    root, ctrl, cable, _ = build_scene("rigid", dt2)
    rig2 = run_levels(root, ctrl, cable, a.levels, a.ramp, a.hold, a.avg)
    del root, ctrl, cable

    # ---------------- Report ----------------
    print("\n=== Rigidezza flessionale ===")
    print(f"  cedevolezza analitica  c = {c_analytic:.4e} mm/N  (EI = {EI:.4g} N*mm^2)")
    print(f"  cedevolezza misurata   c = {c_meas:.4e} mm/N  "
          f"(rapporto mesh/analitica = {c_meas / c_analytic:.3f})")
    if c_meas <= 0:
        print("  ATTENZIONE: cedevolezza non positiva, controllare verso del tendine.")
    T_contact = g_eff / c_meas if c_meas > 0 else float("nan")
    print(f"  T necessaria per chiudere il gap ({g_eff:.3f} mm) ~ {T_contact:.2f} N")

    print("\n=== TEST A: dipendenza da dt (atteso ~1 sia per F che per dY) ===")
    ratios = []
    for r1, r2 in zip(rig1, rig2):
        if r1["raw_mean"] > 1e-12 and r2["raw_mean"] > 1e-12:
            ratios.append(r1["raw_mean"] / r2["raw_mean"])
            print(f"  T={r1['T_N']:5.2f} N  F(dt)/F(dt/2) = {ratios[-1]:.3f}   "
                  f"dY(dt)/dY(dt/2) = {r1['dY_mm'] / r2['dY_mm']:.3f}")
    verdict_scale = None
    if ratios:
        rho = float(np.median(ratios))
        if abs(rho - 2.0) < 0.3:
            print("  -> rapporto ~2: la forza letta dipende ancora da dt. Correzione NON applicata.")
        elif abs(rho - 1.0) < 0.3:
            verdict_scale = 1.0
            print("  -> rapporto ~1: forza letta indipendente da dt. OK (force_scale = 1)")
        else:
            print(f"  -> rapporto {rho:.2f} non conclusivo: allungare --hold o controllare il contatto")
    else:
        print("  Nessun livello in contatto: aumentare --levels.")

    print("\n=== TEST B: confronto con l'equilibrio (dt = %g) ===" % dt1)
    print("     T [N] | F_teorica [mN] | F_letta [mN] | rapporto")
    out_rows = []
    for r in rig1:
        F_th = max(0.0, 1.5 * yt / L * (r["T_N"] - T_contact))
        F_read = r["raw_mean"]
        ratio = F_read / F_th if F_th > 0 else float("nan")
        print(f"    {r['T_N']:6.2f} | {F_th * 1e3:14.2f} | {F_read * 1e3:12.2f} | {ratio:8.3f}")
        out_rows.append({**r, "F_theory_N": F_th, "F_read_N": F_read, "dt": dt1})
    print("  Atteso rapporto ~1 (entro ~30%): il modello di trave e' approssimato.")

    path = os.path.join(PROJECT_DIR, a.csv)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
        w.writeheader()
        w.writerows(out_rows)
    print(f"\nDati salvati in {path}")
    if verdict_scale is not None:
        print(f"Valore suggerito in modules/config.py: \"force_scale\": {verdict_scale:g}"
              + "  (la divisione per dt e' gia' in modules/contact_force.py)")


if __name__ == "__main__":
    main()
