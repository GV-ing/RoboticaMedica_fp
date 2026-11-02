#!/usr/bin/env python3
# tools/check_scene.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Verifica del blocco A (fisica della scena) in un'unica esecuzione.

  [1] Masse        massa del tubo e dell'oggetto soffice calcolate dalla densita'
  [2] Modulo E     il cambio di youngModulus A RUNTIME (come fa il controller tra
                   carico e scarico) modifica davvero la rigidezza? Con e senza reinit().
                   Atteso: dY(E_unload)/dY(E_load) = 77000/65000 = 1.185
  [3] Cedevolezza  dY/T in moto libero -> tensione di inizio contatto con il gap comune
  [4] Rigido       con il gap comune: tensione di inizio contatto, forza misurata vs
                   equilibrio F = 1.5*(yt/L)*(T - T_c), tempo per passo; confronto tra
                   modelli di collisione completi (Triangle+Line+Point) e ridotti
  [5] Soffice      con lo stesso gap: forza massima raggiungibile fino a --T-max e
                   confronto con la soglia di commutazione (5 mN) e il range 10-50 mN

Il controller e' disattivato: la tensione T e' imposta direttamente (T*dt sul cavo,
convenzione verificata con tools/check_cable.py).

Esecuzione (cartella SOFA_Project, variabili d'ambiente del README), ~15-20 min:
    python3.12 -u tools/check_scene.py
    python3.12 -u tools/check_scene.py --skip-full            # gap di config.py
    python3.12 -u tools/check_scene.py --gap 0.05 --skip-full # gap comune imposto
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

FULL_MODELS = {"tube_collision_models": ("Triangle", "Line", "Point"),
               "object_collision_models": ("Triangle", "Line", "Point")}


def build(**overrides):
    p = NOMINAL_PARAMS.copy()
    p.update(switched_mode=False, live_plot=False)
    p.update(overrides)
    root = Sofa.Core.Node("root")
    main_scene.createScene(root, p)
    ctrl = root.getObject("ActuationController")
    ctrl.update_actuation_physics = lambda _dt: None  # T e E imposti dallo script
    tube = root.getChild("NitinolTube")
    cable = tube.getChild("Tendon").getObject("tendonCable")
    cable.value.value = [0.0]
    Sofa.Simulation.init(root)
    return root, ctrl, tube, cable, p


def hold(root, ctrl, cable, T_from, T_to, t_ramp, t_hold, t_avg, with_force=False):
    """Rampa T_from -> T_to, mantenimento, medie finali di dY e forza di contatto."""
    dt = float(root.dt.value)
    n_ramp = max(1, int(round(t_ramp / dt)))
    n_tot = n_ramp + max(1, int(round(t_hold / dt)))
    n_avg = max(1, int(round(t_avg / dt)))
    y0 = ctrl._y0
    dys, fs = [], []
    t0 = time.perf_counter()
    for i in range(n_tot):
        T = T_from + (T_to - T_from) * min(1.0, (i + 1) / n_ramp)
        cable.value.value = [float(T) * dt]
        Sofa.Simulation.animate(root, dt)
        if i >= n_tot - n_avg:
            dys.append(float(ctrl.get_tip_position()[1]) - y0)
            if with_force:
                fs.append(read_contact_force(root))
    t_step = (time.perf_counter() - t0) / n_tot
    out = {"dY": float(np.mean(dys)), "dY_ptp": float(np.ptp(dys)), "t_step": t_step}
    if with_force:
        out.update(F=float(np.mean(fs)), F_ptp=float(np.ptp(fs)))
    return out


def sweep(root, ctrl, cable, levels, a, with_force):
    ctrl._y0 = float(ctrl.get_tip_position()[1])
    rows, T_prev = [], 0.0
    for T in levels:
        r = hold(root, ctrl, cable, T_prev, T, a.ramp, a.hold, a.avg, with_force)
        r["T"] = T
        rows.append(r)
        T_prev = T
        extra = (f"  F={r['F'] * 1e3:7.2f} mN (p-p {r['F_ptp'] * 1e3:.2f})"
                 if with_force else "")
        print(f"    T={T:5.2f} N  dY={r['dY']:+.4f} mm{extra}  {r['t_step'] * 1e3:6.1f} ms/passo")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--T-contact", type=float, default=5.0,
                    help="Tensione a cui deve iniziare il contatto [N] (meta' corsa)")
    ap.add_argument("--T-max", type=float, default=10.0, help="Tensione massima [N]")
    ap.add_argument("--ramp", type=float, default=0.3)
    ap.add_argument("--hold", type=float, default=1.0)
    ap.add_argument("--avg", type=float, default=0.3)
    ap.add_argument("--skip-full", action="store_true",
                    help="Salta il confronto con i modelli di collisione ridotti/completi")
    ap.add_argument("--gap", type=float, default=None,
                    help="Gap comune per rigido e soffice [mm] (default: obj_gap di config.py; "
                         "con --from-T-contact viene ricalcolato da --T-contact)")
    ap.add_argument("--from-T-contact", action="store_true",
                    help="Calcola il gap dalla tensione di inizio contatto --T-contact")
    ap.add_argument("--only-soft", action="store_true",
                    help="Esegue solo le sezioni [3] e [5]")
    a = ap.parse_args()

    P = NOMINAL_PARAMS
    L, yt = P["length_flex"], P["yt_offset"]
    r_o = P["od"] / 2.0
    r_i = r_o - P["thickness"]
    EI = P["E_load"] * math.pi / 4.0 * (r_o ** 4 - r_i ** 4)
    c_beam = yt * L ** 2 / (2.0 * EI)
    summary = {}

    # ---------------- [1] masse + [2] modulo E a runtime ----------------
    T_E = 8.0
    if not a.only_soft:
        print("\n[1-2] Masse e cambio di E a runtime (moto libero, T = 8 N)")
    for use_reinit in (() if a.only_soft else (False, True)):
        root, ctrl, tube, cable, p = build(contact_object="none")
        if not use_reinit:
            m = float(np.ravel(tube.getObject("mass").totalMass.value)[0])
            print(f"    massa tubo = {m:.3e} t = {m * 1e6:.3f} g")
        fem = tube.getObject("FEM")
        ctrl._y0 = float(ctrl.get_tip_position()[1])
        r1 = hold(root, ctrl, cable, 0.0, T_E, a.ramp, a.hold, a.avg)
        fem.youngModulus.value = [P["E_unload"]]
        if use_reinit:
            fem.reinit()
        r2 = hold(root, ctrl, cable, T_E, T_E, a.ramp, a.hold, a.avg)
        ratio = r2["dY"] / r1["dY"]
        tag = "con reinit()" if use_reinit else "senza reinit()"
        print(f"    {tag:15s}: dY(E_load)={r1['dY']:.4f}  dY(E_unload)={r2['dY']:.4f}  "
              f"rapporto={ratio:.3f}  (atteso {P['E_load'] / P['E_unload']:.3f})")
        summary["E_ratio_reinit" if use_reinit else "E_ratio_plain"] = ratio
        del root, ctrl, tube, cable, fem

    # ---------------- [3] cedevolezza e gap ----------------
    print("\n[3] Cedevolezza in moto libero")
    root, ctrl, tube, cable, p = build(contact_object="none")
    levels = [2.0, 4.0, 6.0, 8.0, a.T_max]
    free = sweep(root, ctrl, cable, levels, a, with_force=False)
    del root, ctrl, tube, cable
    T_arr = np.array([r["T"] for r in free])
    dY_arr = np.array([r["dY"] for r in free])
    c = float(np.dot(T_arr, dY_arr) / np.dot(T_arr, T_arr))
    gap_from_T = c * a.T_contact + P["contact_distance"]
    if a.from_T_contact:
        gap = gap_from_T
    else:
        gap = P["obj_gap"] if a.gap is None else a.gap
    T_c = max(0.0, gap - P["contact_distance"]) / c
    print(f"    c = {c:.5f} mm/N (trave: {c_beam:.5f}, rapporto {c / c_beam:.3f}); "
          f"dY a T_max = {c * a.T_max:.3f} mm")
    print(f"    gap usato = {gap:.4f} mm -> contatto atteso a T ~ {T_c:.2f} N "
          f"(per contatto a {a.T_contact:g} N servirebbe {gap_from_T:.4f} mm)")

    # ---------------- [4] rigido ----------------
    contact_levels = [1.0, 2.0, 4.0, 6.0, 8.0, a.T_max]
    REDUCED = {"tube_collision_models": ("Point",), "object_collision_models": ("Triangle",)}
    configs = [] if a.only_soft else (
        [("completi", FULL_MODELS)] + ([] if a.skip_full else [("ridotti", REDUCED)]))
    for label, models in configs:
        print(f"\n[4] Oggetto rigido, modelli di collisione {label}, obj_gap = {gap:.4f} mm")
        root, ctrl, tube, cable, p = build(contact_object="rigid", obj_gap=gap, **models)
        rows = sweep(root, ctrl, cable, contact_levels, a, with_force=True)
        del root, ctrl, tube, cable
        g_eff = gap - P["contact_distance"]
        print("     T [N] | F_teorica [mN] | F_letta [mN] | rapporto")
        for r in rows:
            F_th = max(0.0, 1.5 * yt / L * (r["T"] - g_eff / c))
            rat = r["F"] / F_th if F_th > 1e-6 else float("nan")
            print(f"    {r['T']:6.2f} | {F_th * 1e3:14.2f} | {r['F'] * 1e3:12.2f} | {rat:8.3f}")
        in_contact = [r for r in rows if r["F"] > 1e-4]
        if in_contact:
            summary[f"t_contact_{label}"] = float(np.mean([r["t_step"] for r in in_contact]))
            summary[f"F_rigid_{label}"] = [r["F"] for r in rows]

    # ---------------- [5] soffice ----------------
    print(f"\n[5] Oggetto soffice (E = {P['E_soft'] * 1e3:g} kPa), obj_gap = {gap:.4f} mm")
    root, ctrl, tube, cable, p = build(contact_object="soft", obj_gap=gap)
    soft = root.getChild("SoftContactObject")
    m_soft = float(np.ravel(soft.getObject("mass").totalMass.value)[0])
    print(f"    massa oggetto soffice = {m_soft:.3e} t = {m_soft * 1e6:.2f} g")
    rows_soft = sweep(root, ctrl, cable, contact_levels, a, with_force=True)
    del root, ctrl, tube, cable, soft
    r_soft_max = max(rows_soft, key=lambda r: r["F"])
    F_soft_max = r_soft_max["F"]

    # ---------------- Riepilogo ----------------
    print("\n=== Riepilogo blocco A ===")
    er, ei = summary.get("E_ratio_plain", float("nan")), summary.get("E_ratio_reinit", float("nan"))
    target = P["E_load"] / P["E_unload"]
    if a.only_soft:
        pass
    elif abs(er - target) < 0.03:
        print(f"  E a runtime: FUNZIONA senza reinit (rapporto {er:.3f})")
    elif abs(ei - target) < 0.03:
        print(f"  E a runtime: serve reinit() (senza {er:.3f}, con {ei:.3f}); "
              f"il controller lo chiama dopo ogni cambio di E")
    else:
        print(f"  E a runtime: NON verificato (senza {er:.3f}, con {ei:.3f}): da indagare")
    print(f"  Gap comune usato: {gap:.4f} mm (contatto atteso a T ~ {T_c:.2f} N)")
    tr = summary.get("t_contact_ridotti")
    tf = summary.get("t_contact_completi")
    if tr is not None and tf is not None:
        dev = max(abs(x - y) / max(y, 1e-9) for x, y in
                  zip(summary["F_rigid_ridotti"], summary["F_rigid_completi"]) if y > 1e-4)
        print(f"  Contatto rigido: {tf * 1e3:.0f} ms/passo (completi) -> {tr * 1e3:.0f} ms/passo "
              f"(ridotti), scarto massimo sulla forza {dev:.1%}")
    elif tr is not None:
        print(f"  Contatto rigido (modelli ridotti): {tr * 1e3:.0f} ms/passo")
    th = P["F_contact_th"]
    stable_soft = all(r["F_ptp"] <= 0.3 * r["F"] for r in rows_soft if r["F"] > 1e-4)
    print(f"  Soffice: oscillazione p-p < 30% della media a tutti i livelli in contatto: "
          f"{'si' if stable_soft else 'NO'}")
    print(f"  Soffice: forza massima {F_soft_max * 1e3:.2f} mN a T = {r_soft_max['T']:g} N "
          f"(soglia {th * 1e3:g} mN: {'raggiunta' if F_soft_max >= th else 'NON raggiunta'}; "
          f"range 10-50 mN: {'raggiunto' if F_soft_max >= 0.010 else 'NON raggiunto'})")


if __name__ == "__main__":
    main()
