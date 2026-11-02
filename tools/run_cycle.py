#!/usr/bin/env python3
# tools/run_cycle.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Esegue i casi NOMINALI con il controller completo (molla Dahl, limiti del glossario) e
calcola le metriche del Progetto 5. Serve a verificare modello e metriche prima del
Monte Carlo.

Casi:
  free     ciclo libero 0 -> T_max -> 0 senza oggetto        -> H [%]
           configurazioni: fixed (k = 1.2), switched (0.6/4.8), direct (k = 50)
  contact  carico (contact_rate) fino a contact_stop_force o a fine corsa (T_max), poi u fermo:
           contact_settle_s di assestamento + hold_duration s di valutazione
           su oggetto rigido e soffice                         -> stabilita' del contatto

Ogni caso salva la cronologia in results/<caso>.csv.

Esecuzione (cartella SOFA_Project, source env.sh):
    python3.12 -u tools/run_cycle.py --free                  # 3 cicli liberi (~5-10 min)
    python3.12 -u tools/run_cycle.py --contact               # 4 prove di contatto (~30-40 min)
    python3.12 -u tools/run_cycle.py --free --contact
"""

import argparse
import os
import sys
import time

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

import Sofa.Core  # noqa: E402
import Sofa.Simulation  # noqa: E402

from main_scene import createScene  # noqa: E402
from modules.config import NOMINAL_PARAMS  # noqa: E402
from modules.metrics import contact_stability, hysteresis_H  # noqa: E402

CONFIGS = {
    "fixed": {"switched_mode": False},
    "switched": {"switched_mode": True},
    "direct": {"switched_mode": False, "k_spring": 50.0},
}


def simulate(overrides, hold_s, stop_after_hold, max_time=120.0):
    p = NOMINAL_PARAMS.copy()
    p.update(live_plot=False)
    p.update(overrides)
    root = Sofa.Core.Node("root")
    createScene(root, p)
    Sofa.Simulation.init(root)
    ctrl = root.getObject("ActuationController")
    dt = float(root.dt.value)
    ctrl.start_autonomous_cycle(hold_s=hold_s, repeat=False)
    seen_hold = False
    n = 0
    t0 = time.perf_counter()
    while ctrl.phase != "done" and n * dt < max_time:
        Sofa.Simulation.animate(root, dt)
        n += 1
        seen_hold |= ctrl.phase == "hold"
        if stop_after_hold and seen_hold and ctrl.phase == "unload":
            break
    elapsed = time.perf_counter() - t0
    if n * dt >= max_time:
        print(f"    ATTENZIONE: raggiunto max_time = {max_time} s (fase {ctrl.phase})")
    if ctrl.contact_read_errors:
        print(f"    ATTENZIONE: {ctrl.contact_read_errors} errori di lettura della forza")
    return root, ctrl, n, elapsed


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--free", action="store_true", help="Cicli liberi per H")
    ap.add_argument("--contact", action="store_true", help="Prove di stabilita' del contatto")
    ap.add_argument("--configs", nargs="+", default=list(CONFIGS), choices=list(CONFIGS))
    ap.add_argument("--objects", nargs="+", default=["rigid", "soft"], choices=["rigid", "soft"])
    ap.add_argument("--outdir", default=os.path.join(PROJECT_DIR, "results"))
    a = ap.parse_args()
    if not (a.free or a.contact):
        ap.error("specificare --free e/o --contact")
    P = NOMINAL_PARAMS
    H = {}

    if a.free:
        print("\n=== Cicli liberi (H) ===")
        for name in a.configs:
            root, ctrl, n, el = simulate({"contact_object": "none", **CONFIGS[name]},
                                         hold_s=0.0, stop_after_hold=False)
            path = ctrl.export_csv(os.path.join(a.outdir, f"free_{name}.csv"))
            h = ctrl.history
            r = hysteresis_H(h["u"], h["dY"], h["phase"])
            H[name] = r["H_pct"]
            print(f"  {name:8s}: H = {r['H_pct']:6.2f} %  dY_max = {r['dY_max_mm']:.4f} mm  "
                  f"u_max = {r['u_max_mm']:.3f} mm  scarto max {r['gap_max_mm']:.4f} mm "
                  f"a u = {r['u_at_gap_max_mm']:.3f} mm  T_max raggiunta = {max(h['T']):.3f} N")
            print(f"            {n} passi, {el:.0f} s ({el / n * 1e3:.0f} ms/passo) -> {path}")
            del root, ctrl
        if "direct" in H:
            for name in ("fixed", "switched"):
                if name in H:
                    print(f"  H_{name} / H_direct = {H[name] / H['direct']:.3f}")

    if a.contact:
        print(f"\n=== Prove di contatto (arresto a {P['contact_stop_force'] * 1e3:g} mN o a fine corsa,"
              f" {P['contact_settle_s']:g} s + {P['hold_duration']:g} s fermo) ===")
        for obj in a.objects:
            for name in [c for c in a.configs if c != "direct"]:
                overrides = {"contact_object": obj, "auto_stop_force": P["contact_stop_force"],
                             "auto_rate": P["contact_rate"], **CONFIGS[name]}
                root, ctrl, n, el = simulate(overrides,
                                             hold_s=P["contact_settle_s"] + P["hold_duration"],
                                             stop_after_hold=True)
                path = ctrl.export_csv(os.path.join(a.outdir, f"contact_{obj}_{name}.csv"))
                h = ctrl.history
                s = contact_stability(h["t"], h["F_contact"], h["phase"],
                                      duration=P["hold_duration"],
                                      ptp_max=P["stability_ptp_max"],
                                      loss_force=P["contact_loss_force"])
                hold = [i for i, ph in enumerate(h["phase"]) if ph == "hold"]
                F_stop = h["F_contact"][hold[0]] if hold else float("nan")
                F_hold_max = max(h["F_contact"][i] for i in hold) if hold else float("nan")
                print(f"  {obj:5s} {name:8s}: arresto a F = {F_stop * 1e3:.2f} mN, T = "
                      f"{h['T'][hold[0]] if hold else float('nan'):.2f} N; max nel mantenimento "
                      f"{F_hold_max * 1e3:.2f} mN")
                print(f"  {'':5s} {'':8s}  F media = {s['F_mean_N'] * 1e3:7.2f} mN  "
                      f"p-p/media = {s['ptp_over_mean']:.3f}  perso = {s['contact_lost']}  "
                      f"durata = {s['hold_duration_s']:.2f} s  -> stabile = {s['stable']}  "
                      f"(range 10-50 mN: {s['in_range_10_50mN']}, k finale {h['k'][-1]:.2f})")
                print(f"            {n} passi, {el:.0f} s ({el / n * 1e3:.0f} ms/passo) -> {path}")
                del root, ctrl


if __name__ == "__main__":
    main()
