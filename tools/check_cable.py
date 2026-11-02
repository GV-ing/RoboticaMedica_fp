#!/usr/bin/env python3
# tools/check_cable.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Diagnosi del carico del tendine: confronta la freccia di punta in moto libero
prodotta dallo STESSO momento flettente M applicato in quattro modi diversi.

  1. cavo_nominale        CableConstraint come in main_scene (yt=0.6, pullPoint in x=40)
  2. pullpoint_arretrato  come 1, ma pullPoint spostato in x=0 (primo segmento non nullo)
  3. cavo_nella_parete    punti del cavo a meta' parete (y=0.855) invece che nel lume vuoto
  4. momento_puro         nessun cavo: coppia di forze assiali sui nodi di punta
                          (riferimento FEM puro, resultante nulla, momento M)

Per una mensola con momento in punta:  delta_teorica = M * L^2 / (2 * E * I)
Il caso 4 da' la risposta "vera" della mesh (incluso il locking); se i casi 1-3
se ne discostano molto, il problema e' nel modo in cui il cavo carica il tubo.

Esecuzione:
    python3.12 -u tools/check_cable.py
    python3.12 -u tools/check_cable.py --moment 1.2 --nx 30
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

VARIANTS = ("cavo_nominale", "pullpoint_arretrato", "cavo_nella_parete", "momento_puro",
            "forza_puntuale", "compressione_distribuita")
# forza_puntuale: niente CableConstraint; forza -T lungo x sull'ultimo punto del
#     cavo (stessa mappatura baricentrica): stesso carico "statico" del cavo dritto.
# compressione_distribuita: niente cavo; sull'anello di punta compressione
#     uniforme -T piu' il momento T*yt: stesso carico globale, ma distribuito.


def run_variant(variant, moment, n_x, dt, t_ramp, t_hold, t_avg):
    p = NOMINAL_PARAMS.copy()
    if n_x is not None:
        p["mesh_h_flex"] = p["length_flex"] / n_x
    r_o = p["od"] / 2.0
    r_i = r_o - p["thickness"]
    if variant == "cavo_nella_parete":
        p["yt_offset"] = 0.5 * (r_i + r_o)
    p.update(contact_object="none", dt=dt, switched_mode=False, live_plot=False)
    arm = p["yt_offset"]
    T = moment / arm  # tensione che produce il momento M con questo braccio

    root = Sofa.Core.Node("root")
    main_scene.createScene(root, p)
    ctrl = root.getObject("ActuationController")
    ctrl.update_actuation_physics = lambda _dt: None
    tube = root.getChild("NitinolTube")
    cable = tube.getChild("Tendon").getObject("tendonCable")
    cable.value.value = [0.0]

    if variant == "pullpoint_arretrato":
        cable.pullPoint.value = [0.0, float(arm), 0.0]

    cff = None
    unit_forces = None
    tip_idx = list(ctrl.tip_node_indices)
    if variant in ("momento_puro", "compressione_distribuita", "forza_puntuale"):
        root.addObject("RequiredPlugin", pluginName="Sofa.Component.MechanicalLoad")
    if variant in ("momento_puro", "compressione_distribuita"):
        rest = np.asarray(tube.getObject("topology").position.value, dtype=float)
        y = rest[tip_idx, 1]
        # f_x,i = -M * y_i / sum(y^2): risultante nulla, momento Mz = +M (flette verso +Y)
        unit_forces = np.zeros((len(tip_idx), 3))
        unit_forces[:, 0] = -y / np.sum(y ** 2)
        if variant == "compressione_distribuita":
            # + compressione uniforme T = M/yt divisa sui nodi dell'anello (per unita' di M)
            unit_forces[:, 0] += -(1.0 / arm) / len(tip_idx)
        cff = tube.addObject("ConstantForceField", name="tipLoad", indices=tip_idx,
                             forces=(0.0 * unit_forces).tolist())
    elif variant == "forza_puntuale":
        tendon = tube.getChild("Tendon")
        n_pts = len(tendon.getObject("cablePts").position.value)
        unit_forces = np.array([[-1.0 / arm, 0.0, 0.0]])  # per unita' di M: T = M/yt
        cff = tendon.addObject("ConstantForceField", name="tipPointLoad",
                               indices=[n_pts - 1], forces=[[0.0, 0.0, 0.0]])

    Sofa.Simulation.init(root)

    y0 = float(ctrl.get_tip_position()[1])
    n_ramp = max(1, int(round(t_ramp / dt)))
    n_tot = n_ramp + max(1, int(round(t_hold / dt)))
    n_avg = max(1, int(round(t_avg / dt)))
    dys = []
    for i in range(n_tot):
        s = min(1.0, (i + 1) / n_ramp)
        if cff is None:
            cable.value.value = [float(s * T) * dt]  # impulso: forza s*T
        else:
            cff.forces.value = (s * moment * unit_forces).tolist()
        Sofa.Simulation.animate(root, dt)
        if i >= n_tot - n_avg:
            dys.append(float(ctrl.get_tip_position()[1]) - y0)

    disp = None
    if cff is None:
        try:
            disp = float(np.ravel(cable.displacement.value)[0])
        except Exception:
            disp = None
    dys = np.asarray(dys)
    return {"T": T if cff is None else 0.0, "arm": arm, "dY": float(dys.mean()),
            "dY_ptp": float(np.ptp(dys)), "cable_disp": disp}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--moment", type=float, default=1.2, help="Momento in punta [N*mm]")
    ap.add_argument("--nx", type=int, default=None,
                    help="Elementi sulla parte flessibile (default: mesh_h_flex di config.py)")
    ap.add_argument("--dt", type=float, default=NOMINAL_PARAMS["dt"])
    ap.add_argument("--ramp", type=float, default=0.5)
    ap.add_argument("--hold", type=float, default=2.0)
    ap.add_argument("--avg", type=float, default=0.5)
    ap.add_argument("--variants", nargs="+", default=list(VARIANTS), choices=VARIANTS)
    a = ap.parse_args()

    p = NOMINAL_PARAMS
    L = p["length_flex"]
    r_o = p["od"] / 2.0
    r_i = r_o - p["thickness"]
    EI = p["E_load"] * math.pi / 4.0 * (r_o ** 4 - r_i ** 4)
    d_th = a.moment * L ** 2 / (2.0 * EI)
    print(f"\nMomento M={a.moment:g} N*mm, n_flex={a.nx or round(p['length_flex'] / p['mesh_h_flex'])}: freccia teorica = {d_th:.4f} mm")

    res = {}
    try:
        for v in a.variants:
            print(f"\n--- {v} ---")
            r = run_variant(v, a.moment, a.nx, a.dt, a.ramp, a.hold, a.avg)
            res[v] = r
            extra = (f"  T={r['T']:.3f} N (braccio {r['arm']:.3f} mm)"
                     f"  spostamento cavo={r['cable_disp']}" if r["T"] > 0 else "")
            print(f"  dY={r['dY']:+.4f} mm (p-p {r['dY_ptp']:.1e})  "
                  f"SOFA/teoria={r['dY'] / d_th:.3f}{extra}")
    finally:
        pass

    print("\n=== Riepilogo: freccia / freccia del momento puro ===")
    ref = res.get("momento_puro", {}).get("dY")
    for v, r in res.items():
        rel = f"{r['dY'] / ref:.2f}" if ref else "n/d"
        print(f"  {v:>20}: dY={r['dY']:+.4f} mm  rispetto al momento puro = {rel}")
    print("  Atteso ~1 per un cavo che applica correttamente il momento T*yt.")


if __name__ == "__main__":
    main()
