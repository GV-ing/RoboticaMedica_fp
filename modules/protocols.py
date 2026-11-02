# modules/protocols.py
"""Protocolli di prova del Progetto 5 eseguiti in SOFA (senza GUI).

  free_cycle    ciclo libero 0 -> T_max (o u_max) -> 0, senza oggetto  -> H [%]
  contact_test  carico a contact_rate fino a contact_stop_force o a fine corsa, poi u fermo
                per contact_settle_s + hold_duration; stabilita' sulla finestra finale

Usato da tools/run_cycle.py e dal Monte Carlo (modules/mc_worker.py).
"""

import time

from modules.config import NOMINAL_PARAMS
from modules.metrics import contact_stability, hysteresis_H

CONFIGS = {
    "fixed": {"switched_mode": False},
    "switched": {"switched_mode": True},
    "direct": {"switched_mode": False, "k_spring": 50.0},
}


def simulate(overrides, hold_s, stop_after_hold, max_time=120.0):
    """Costruisce la scena, esegue un ciclo autonomo e restituisce (root, ctrl, passi, s)."""
    import Sofa.Core
    import Sofa.Simulation

    from main_scene import createScene

    p = NOMINAL_PARAMS.copy()
    p.update(live_plot=False)
    p.update(overrides)
    root = Sofa.Core.Node("root")
    createScene(root, p)
    Sofa.Simulation.init(root)
    ctrl = root.getObject("ActuationController")
    dt = float(root.dt.value)
    ctrl.start_autonomous_cycle(hold_s=hold_s, repeat=False)
    seen_hold, n = False, 0
    t0 = time.perf_counter()
    while ctrl.phase != "done" and n * dt < max_time:
        Sofa.Simulation.animate(root, dt)
        n += 1
        seen_hold |= ctrl.phase == "hold"
        if stop_after_hold and seen_hold and ctrl.phase == "unload":
            break
    timed_out = n * dt >= max_time
    return root, ctrl, n, time.perf_counter() - t0, timed_out


def free_cycle(overrides):
    """Ciclo libero: restituisce (metriche, cronologia)."""
    root, ctrl, n, elapsed, timed_out = simulate(
        {"contact_object": "none", **overrides}, hold_s=0.0, stop_after_hold=False)
    h = {k: list(v) for k, v in ctrl.history.items()}
    out = hysteresis_H(h["u"], h["dY"], h["phase"])
    out.update(T_peak_N=float(max(h["T"])), steps=n, elapsed_s=elapsed, timed_out=timed_out)
    return out, h


def contact_test(obj, overrides):
    """Prova di stabilita' del contatto su obj ('rigid' | 'soft')."""
    P = {**NOMINAL_PARAMS, **overrides}
    o = {"contact_object": obj, "auto_stop_force": P["contact_stop_force"],
         "auto_rate": P["contact_rate"], **overrides}
    root, ctrl, n, elapsed, timed_out = simulate(
        o, hold_s=P["contact_settle_s"] + P["hold_duration"], stop_after_hold=True)
    h = {k: list(v) for k, v in ctrl.history.items()}
    out = contact_stability(h["t"], h["F_contact"], h["phase"], duration=P["hold_duration"],
                            ptp_max=P["stability_ptp_max"], loss_force=P["contact_loss_force"])
    hold = [i for i, ph in enumerate(h["phase"]) if ph == "hold"]
    out.update(
        F_stop_N=float(h["F_contact"][hold[0]]) if hold else float("nan"),
        T_stop_N=float(h["T"][hold[0]]) if hold else float("nan"),
        F_hold_max_N=float(max(h["F_contact"][i] for i in hold)) if hold else float("nan"),
        k_final=float(h["k"][-1]), steps=n, elapsed_s=elapsed, timed_out=timed_out,
        contact_read_errors=int(ctrl.contact_read_errors))
    return out, h
