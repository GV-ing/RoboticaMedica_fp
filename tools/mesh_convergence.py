#!/usr/bin/env python3
# tools/mesh_convergence.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Studio di convergenza della mesh del tubo (TetrahedronFEMForceField).

Per ogni passo assiale h_flex della parte flessibile (mesh graduata: pochi elementi
sull'albero vincolato, superfici visiva/collisione separate e piu' grossolane):
  1. CAVO  : tensione T sul CableConstraint (carico reale del progetto, momento T*yt)
             teoria: delta = T*yt*L^2 / (2*E*I)
  2. FORZA : forza laterale F sui nodi di punta (ConstantForceField)
             teoria: delta = F*L^3 / (3*E*I)
Riporta il rapporto SOFA/teoria, l'estrapolazione di Richardson sulle tre mesh
piu' fitte, il tempo medio per passo (moto libero e, se richiesto, con contatto
rigido) e una stima della durata della campagna Monte Carlo.

Convenzione verificata durante lo sviluppo: il valore del CableConstraint e'
un impulso, quindi si impone T*dt.

Esecuzione (cartella SOFA_Project, variabili d'ambiente del README):
    python3.12 -u tools/mesh_convergence.py --contact-timing
    python3.12 -u tools/mesh_convergence.py --hflex 0.8 0.4 0.3 0.2 --contact-timing
"""

import argparse
import csv
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

# Campagna attuale: 4 scenari x (1 nominale + 20 realizzazioni), 1600 passi a caso
MC_CASES = 4 * 21
MC_STEPS_PER_CASE = 1600


def build(h_flex, n_theta, contact_object, dt, load):
    p = NOMINAL_PARAMS.copy()
    p.update(contact_object=contact_object, dt=dt, switched_mode=False, live_plot=False,
             mesh_h_flex=h_flex, mesh_n_theta=n_theta)
    root = Sofa.Core.Node("root")
    main_scene.createScene(root, p)

    ctrl = root.getObject("ActuationController")
    ctrl.update_actuation_physics = lambda _dt: None  # carico imposto dallo script
    tube = root.getChild("NitinolTube")
    cable = tube.getChild("Tendon").getObject("tendonCable")
    cable.value.value = [0.0]

    cff = None
    if load == "force":
        root.addObject("RequiredPlugin", pluginName="Sofa.Component.MechanicalLoad")
        cff = tube.addObject("ConstantForceField", name="tipLoad",
                             indices=list(ctrl.tip_node_indices),
                             totalForce=[0.0, 0.0, 0.0])
    Sofa.Simulation.init(root)
    info = {
        "n_nodes": len(tube.getObject("dofs").position.value),
        "n_tets": len(tube.getObject("topology").tetrahedra.value),
    }
    return root, ctrl, cable, cff, info


def run_static(root, ctrl, cable, cff, magnitude, t_ramp, t_hold, t_avg):
    """Rampa del carico, mantenimento, media di dY; misura il tempo per passo."""
    dt = float(root.dt.value)
    n_ramp = max(1, int(round(t_ramp / dt)))
    n_hold = max(1, int(round(t_hold / dt)))
    n_avg = min(n_hold, max(1, int(round(t_avg / dt))))
    y0 = float(ctrl.get_tip_position()[1])
    dys = []
    t0 = time.perf_counter()
    for i in range(n_ramp + n_hold):
        s = min(1.0, (i + 1) / n_ramp)
        if cff is None:
            cable.value.value = [float(s * magnitude) * dt]  # impulso -> forza s*T
        else:
            cff.totalForce.value = [0.0, float(s * magnitude), 0.0]
        Sofa.Simulation.animate(root, dt)
        if i >= n_ramp + n_hold - n_avg:
            dys.append(float(ctrl.get_tip_position()[1]) - y0)
    t_step = (time.perf_counter() - t0) / (n_ramp + n_hold)
    dys = np.asarray(dys)
    return float(dys.mean()), float(np.ptp(dys)), t_step


def time_contact(h_flex, n_theta, dt, T, n_steps):
    """Tempo per passo con oggetto rigido e contatto attivo (tendine a T)."""
    root, ctrl, cable, _, _ = build(h_flex, n_theta, "rigid", dt, "cable")
    for _ in range(int(round(1.0 / dt))):  # porta la punta in contatto
        cable.value.value = [T * dt]
        Sofa.Simulation.animate(root, dt)
    t0 = time.perf_counter()
    for _ in range(n_steps):
        cable.value.value = [T * dt]
        Sofa.Simulation.animate(root, dt)
    return (time.perf_counter() - t0) / n_steps


def richardson(h, q):
    """Estrapola q(h->0) dalle ultime tre mesh: q(h) = q_inf - C*h^p (h anche non uniforme)."""
    if len(q) < 3:
        return None, None
    (h1, h2, h3), (q1, q2, q3) = h[-3:], q[-3:]
    if not (h1 > h2 > h3 > 0) or (q2 - q1) == 0 or (q3 - q2) == 0:
        return None, None
    target = (q2 - q1) / (q3 - q2)
    if target <= 0:
        return None, None  # convergenza non monotona

    def g(p):
        return (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - target

    lo, hi = 0.05, 8.0
    if g(lo) * g(hi) > 0:
        return None, None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    order = 0.5 * (lo + hi)
    C = (q3 - q2) / (h2 ** order - h3 ** order)
    return q3 + C * h3 ** order, order


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hflex", type=float, nargs="+", default=[0.8, 0.4, 0.3, 0.2],
                    help="Passi assiali della parte flessibile [mm]")
    ap.add_argument("--ntheta", type=int, default=NOMINAL_PARAMS["mesh_n_theta"])
    ap.add_argument("--T", type=float, default=10.0, help="Tensione cavo [N]")
    ap.add_argument("--F", type=float, default=0.1, help="Forza laterale in punta [N]")
    ap.add_argument("--dt", type=float, default=NOMINAL_PARAMS["dt"])
    ap.add_argument("--ramp", type=float, default=0.5)
    ap.add_argument("--hold", type=float, default=2.0)
    ap.add_argument("--avg", type=float, default=0.5)
    ap.add_argument("--contact-timing", action="store_true",
                    help="Misura anche il tempo per passo con contatto rigido")
    ap.add_argument("--contact-T", type=float, default=12.0,
                    help="Tensione usata per la prova di tempo con contatto [N]")
    ap.add_argument("--tol", type=float, default=0.05, help="Errore ammesso (5%%)")
    ap.add_argument("--csv", default="mesh_convergence.csv")
    a = ap.parse_args()

    p = NOMINAL_PARAMS
    L, yt = p["length_flex"], p["yt_offset"]
    r_o = p["od"] / 2.0
    r_i = r_o - p["thickness"]
    EI = p["E_load"] * math.pi / 4.0 * (r_o ** 4 - r_i ** 4)
    d_cable = a.T * yt * L ** 2 / (2.0 * EI)
    d_force = a.F * L ** 3 / (3.0 * EI)
    print(f"\nTeoria: cavo T={a.T:g} N -> {d_cable:.4f} mm ;  forza F={a.F:g} N -> {d_force:.4f} mm")

    rows = []
    for h in sorted(a.hflex, reverse=True):
        n_x = int(round(L / h))  # elementi sulla parte flessibile
        print(f"\n--- h_flex={h:g} mm ({n_x} elementi flessibili), n_theta={a.ntheta} ---")
        root, ctrl, cable, cff, info = build(h, a.ntheta, "none", a.dt, "cable")
        dy_c, ptp_c, ts_free = run_static(root, ctrl, cable, None, a.T,
                                          a.ramp, a.hold, a.avg)
        del root, ctrl, cable, cff
        root, ctrl, cable, cff, _ = build(h, a.ntheta, "none", a.dt, "force")
        dy_f, ptp_f, _ = run_static(root, ctrl, cable, cff, a.F, a.ramp, a.hold, a.avg)
        del root, ctrl, cable, cff
        ts_contact = (time_contact(h, a.ntheta, a.dt, a.contact_T, 100)
                      if a.contact_timing else None)
        row = {
            "n_x": n_x, "h_mm": h, **info,
            "dY_cable_mm": dy_c, "ratio_cable": dy_c / d_cable,
            "dY_force_mm": dy_f, "ratio_force": dy_f / d_force,
            "ptp_cable": ptp_c, "ptp_force": ptp_f,
            "t_step_free_ms": ts_free * 1e3,
            "t_step_contact_ms": None if ts_contact is None else ts_contact * 1e3,
        }
        rows.append(row)
        msg = (f"  nodi={info['n_nodes']} tet={info['n_tets']} | "
               f"cavo {dy_c:.4f} mm (x{row['ratio_cable']:.3f}) | "
               f"forza {dy_f:.4f} mm (x{row['ratio_force']:.3f}) | "
               f"{row['t_step_free_ms']:.1f} ms/passo")
        if ts_contact is not None:
            msg += f", contatto {row['t_step_contact_ms']:.1f} ms/passo"
        if max(ptp_c / max(abs(dy_c), 1e-12), ptp_f / max(abs(dy_f), 1e-12)) > 1e-3:
            msg += "  [ATTENZIONE: non ancora statico, aumentare --hold]"
        print(msg)

    # ---------------- Riepilogo ----------------
    print("\n=== Riepilogo (rapporto SOFA/teoria, 1 = trave di Eulero-Bernoulli) ===")
    print("  n_flex | h [mm] |  nodi | cavo  | forza | ms/passo libero | ms/passo contatto | stima MC [h]")
    for r in rows:
        t_ms = r["t_step_contact_ms"] or r["t_step_free_ms"]
        mc_h = t_ms * 1e-3 * MC_CASES * MC_STEPS_PER_CASE / 3600.0
        tc = "-" if r["t_step_contact_ms"] is None else f"{r['t_step_contact_ms']:.1f}"
        print(f"  {r['n_x']:6d} | {r['h_mm']:6.3f} | {r['n_nodes']:5d} | "
              f"{r['ratio_cable']:.3f} | {r['ratio_force']:.3f} | "
              f"{r['t_step_free_ms']:15.1f} | {tc:>17} | {mc_h:8.2f}")
    print("  (stima MC: 84 casi x 1600 passi; usa il tempo con contatto se misurato)")

    h_arr = [r["h_mm"] for r in rows]
    for key, label in (("ratio_cable", "cavo"), ("ratio_force", "forza")):
        q_inf, order = richardson(h_arr, [r[key] for r in rows])
        if q_inf is not None:
            print(f"  Richardson {label}: limite h->0 = {q_inf:.3f}, ordine osservato = {order:.2f}")
            if key == "ratio_cable":
                for r in rows:
                    r["err_vs_limit"] = abs(r[key] / q_inf - 1.0)
                print("    errore rispetto al limite: " + ", ".join(
                    f"h={r['h_mm']:g}: {r['err_vs_limit']:.1%}" for r in rows))
        else:
            print(f"  Richardson {label}: non applicabile (servono 3 mesh con raffinamento "
                  f"costante e convergenza monotona)")

    ok = [r for r in rows if abs(r["ratio_cable"] - 1.0) <= a.tol]
    if ok:
        best = ok[0]
        print(f"\n  Mesh piu' grossolana con errore cavo <= {a.tol:.0%}: h_flex = {best['h_mm']:g} mm "
              f"(rapporto {best['ratio_cable']:.3f})")
    else:
        print(f"\n  Nessuna mesh entro {a.tol:.0%}: aumentare n_x o valutare "
              f"HexahedronFEMForceField / mesh rifinita solo sulla parte flessibile.")

    path = os.path.join(PROJECT_DIR, a.csv)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  Dati salvati in {path}")


if __name__ == "__main__":
    main()
