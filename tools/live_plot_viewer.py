#!/usr/bin/env python3
# tools/live_plot_viewer.py
r"""
TELEMETRIA E PLOT REAL-TIME - Progetto 5 (tubo a tendine con molla in serie)
============================================================================
Struttura e informazioni allineate al Modulo 06 del corso
("Telemetria e Plot Real-Time per Molla ed Isteresi Meccanica"), adattate al tubo:

COLONNA SINISTRA - LOOP DI ISTERESI
  1. Loop della molla T(DeltaLs): confronto tra la legge elastica ideale di Hooke
     (T = k*DeltaLs, verde) e la traiettoria reale con isteresi di Dahl (magenta);
     banda isteretica T = k*DeltaLs +/- (Fh_molla + Fh_mat) ed energia dissipata [mJ].
  2. Loop della punta DeltaY(u): rami di carico e scarico a pari u e indice H
     (massima differenza di DeltaY a pari u, in % di DeltaYmax), come da traccia.

COLONNA DESTRA - DINAMICA TEMPORALE (grafici obbligatori della traccia)
  3. Posizione della punta: DeltaY = Y - Y0 e DeltaX = X - X0  [mm]
  4. Attuazione: u (motore) e DeltaLs (allungamento molla)    [mm]
  5. Tensione del tendine T                                    [N]
  6. Forza di contatto, soglia di commutazione 5 mN e range 10-50 mN  [mN]
  7. Rigidezza della molla k (commutazione istantanea 0.6 <-> 4.8 N/mm)

DIFFERENZA RISPETTO AL MODULO 06 (voluta): il plot NON gira dentro SOFA. Il Python 3.12 di
SOFA non ha matplotlib e il ridisegno rallenterebbe ogni passo; questo script e' un processo
separato (Python di sistema) che riceve un campione JSON per riga da modules/live_plot.py.

CONTROLLI DA TASTIERA (finestra dei grafici):
  - Tasto 'R': reset dei loop e dell'energia dissipata (come il tasto R del Modulo 06)

Uso manuale / prova senza finestra:
    env -u PYTHONPATH python3 tools/live_plot_viewer.py < campioni.jsonl
    env -u PYTHONPATH python3 tools/live_plot_viewer.py --snapshot out.png < campioni.jsonl
"""

import json
import os
import sys
import threading
from collections import deque

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib  # noqa: E402

if "--snapshot" in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.animation import FuncAnimation  # noqa: E402

from modules.config import NOMINAL_PARAMS  # noqa: E402
from modules.metrics import hysteresis_H  # noqa: E402

# ==============================================================================
# PARAMETRI (dal modello: modules/config.py)
# ==============================================================================
F_HYST_TOT = NOMINAL_PARAMS["Fh_spring"] + NOMINAL_PARAMS["Fh_mat"]  # [N] semi-ampiezza banda
F_TH_mN = NOMINAL_PARAMS["F_contact_th"] * 1e3                       # soglia 5 mN
TIME_WINDOW_S = 10.0          # finestra mobile dei grafici temporali [s simulati]
MAX_POINTS = 20000

# Palette scura del Modulo 06
BG_FIG, BG_AX, GRID, TXT, TXT_2 = "#1E1E24", "#141419", "#2D2D38", "#EAEAEA", "#C0C0C0"
C_HOOKE, C_DAHL, C_POINT, C_ENERGY = "#50FA7B", "#FF79C6", "#FF5555", "#F1FA8C"
C_BAND, C_BLUE, C_ORANGE, C_PURPLE = "#BD93F9", "#8BE9FD", "#FFB86C", "#BD93F9"
LEGEND = dict(facecolor="#282A36", edgecolor="#44475A", labelcolor="#F8F8F2", fontsize=8)

TIME_KEYS = ("t", "X", "Y", "dX", "dY", "u", "dLs", "T", "F", "k")
data = {k: deque(maxlen=MAX_POINTS) for k in TIME_KEYS}
loop = {k: [] for k in ("u", "dY", "dLs", "T", "branch")}
state = {"x0": None, "phase": "", "z": 0.0, "E_diss": 0.0, "u_prev": None,
         "branch": "load", "closed": False}
lock = threading.Lock()


def reset_loops():
    """Tasto R: azzera loop ed energia dissipata (la dinamica temporale resta)."""
    with lock:
        for v in loop.values():
            v.clear()
        state["E_diss"] = 0.0


def add_sample(s):
    if state["x0"] is None:
        state["x0"] = s["tip_x"]
    data["t"].append(s["t"])
    data["X"].append(s["tip_x"])
    data["Y"].append(s.get("tip_y", s["dY"]))
    data["dX"].append(s["tip_x"] - state["x0"])
    data["dY"].append(s["dY"])
    data["u"].append(s["u"])
    data["dLs"].append(s["dLs"])
    data["T"].append(s["T"])
    data["F"].append(s["F_contact"] * 1e3)
    data["k"].append(s["k"])
    state["phase"], state["z"] = s.get("phase", ""), s.get("z", 0.0)

    # Ramo del ciclo: fasi del modo autonomo, altrimenti verso di u (modo manuale)
    if s.get("phase") in ("load", "unload"):
        state["branch"] = s["phase"]
    elif state["u_prev"] is not None and abs(s["u"] - state["u_prev"]) > 1e-9:
        state["branch"] = "load" if s["u"] > state["u_prev"] else "unload"
    state["u_prev"] = s["u"]

    # Energia dissipata: lavoro della parte non elastica (isteresi + smorzamento),
    # analogo di sum |f_hyst * dy| del Modulo 06. N*mm = mJ.
    if loop["dLs"]:
        d = s["dLs"] - loop["dLs"][-1]
        f_ne = 0.5 * ((s["T"] - s["k"] * s["dLs"])
                      + (loop["T"][-1] - s["k"] * loop["dLs"][-1]))
        state["E_diss"] += abs(f_ne * d)
    loop["u"].append(s["u"])
    loop["dY"].append(s["dY"])
    loop["dLs"].append(s["dLs"])
    loop["T"].append(s["T"])
    loop["branch"].append(state["branch"])


def reader():
    for line in sys.stdin:
        try:
            s = json.loads(line)
        except json.JSONDecodeError:
            continue
        with lock:
            add_sample(s)
    state["closed"] = True


def style_axis(ax, title, xlabel, ylabel):
    ax.set_facecolor(BG_AX)
    ax.set_title(title, color=TXT, fontsize=10.5, fontweight="bold", loc="left")
    ax.set_xlabel(xlabel, color=TXT_2, fontsize=8.5)
    ax.set_ylabel(ylabel, color=TXT_2, fontsize=8.5)
    ax.tick_params(colors="#A0A0A0", labelsize=8)
    ax.grid(True, linestyle="--", color=GRID, alpha=0.7)
    for sp in ax.spines.values():
        sp.set_color(GRID)


def main():
    fig = plt.figure(figsize=(15, 9))
    fig.patch.set_facecolor(BG_FIG)
    try:
        fig.canvas.manager.set_window_title(
            "Robotica Medica - Progetto 5: Telemetria Real-Time Tubo, Molla ed Isteresi")
    except AttributeError:
        pass
    gs = fig.add_gridspec(1, 2, width_ratios=(1.0, 1.25), wspace=0.16)
    g_left = gs[0].subgridspec(2, 1, hspace=0.32)
    g_right = gs[1].subgridspec(5, 1, hspace=0.42)

    # --- 1. LOOP MOLLA T(DeltaLs) ---
    a_sp = fig.add_subplot(g_left[0])
    style_axis(a_sp, "Loop di isteresi della molla: T vs DeltaLs",
               "Allungamento molla DeltaLs [mm]", "Tensione tendine T [N]")
    l_hooke, = a_sp.plot([], [], color=C_HOOKE, ls=":", lw=1.6, label="Elastico ideale (Hooke) T = k*DeltaLs")
    l_bup, = a_sp.plot([], [], color=C_BAND, ls="--", lw=1.0,
                       label=f"Banda isteretica +/-{F_HYST_TOT:.2f} N (Fh_molla + Fh_mat)")
    l_blo, = a_sp.plot([], [], color=C_BAND, ls="--", lw=1.0)
    l_dahl, = a_sp.plot([], [], color=C_DAHL, lw=2.2, label="Traiettoria con isteresi (Dahl)")
    p_sp, = a_sp.plot([], [], color=C_POINT, marker="o", ms=7, ls="", label="Punto attuale")
    t_energy = a_sp.text(0.03, 0.93, "", transform=a_sp.transAxes, color=C_ENERGY,
                         fontsize=9, fontweight="bold", va="top",
                         bbox=dict(boxstyle="round,pad=0.3", facecolor="#282A36",
                                   edgecolor="#6272A4", alpha=0.9))
    a_sp.legend(loc="lower right", **LEGEND)

    # --- 2. LOOP PUNTA DeltaY(u) e indice H ---
    a_tl = fig.add_subplot(g_left[1])
    style_axis(a_tl, "Loop di isteresi della punta: DeltaY vs u",
               "Corsa motore u [mm]", "DeltaY punta [mm]")
    l_load, = a_tl.plot([], [], color=C_BLUE, lw=2.0, label="Carico (u crescente)")
    l_unl, = a_tl.plot([], [], color=C_ORANGE, lw=2.0, label="Scarico (u decrescente)")
    p_tl, = a_tl.plot([], [], color=C_POINT, marker="o", ms=7, ls="", label="Punto attuale")
    l_gap, = a_tl.plot([], [], color=C_ENERGY, lw=1.4, ls="-", marker="_", ms=10,
                       label="Max differenza a pari u")
    t_H = a_tl.text(0.03, 0.93, "", transform=a_tl.transAxes, color=C_ENERGY, fontsize=9,
                    fontweight="bold", va="top",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#282A36",
                              edgecolor="#6272A4", alpha=0.9))
    a_tl.legend(loc="lower right", **LEGEND)

    # --- 3-7. DINAMICA TEMPORALE ---
    a_tip = fig.add_subplot(g_right[0])
    a_act = fig.add_subplot(g_right[1], sharex=a_tip)
    a_T = fig.add_subplot(g_right[2], sharex=a_tip)
    a_F = fig.add_subplot(g_right[3], sharex=a_tip)
    a_k = fig.add_subplot(g_right[4], sharex=a_tip)
    style_axis(a_tip, "Posizione della punta (X, Y assoluti nella barra in alto)", "", "[mm]")
    style_axis(a_act, "Attuazione: motore e molla", "", "[mm]")
    style_axis(a_T, "Tensione del tendine", "", "T [N]")
    style_axis(a_F, "Forza di contatto", "", "F [mN]")
    style_axis(a_k, "Rigidezza della molla", "Tempo simulato [s]", "k [N/mm]")

    l_dY, = a_tip.plot([], [], color=C_BLUE, lw=1.8, label="DeltaY = Y - Y0")
    l_dX, = a_tip.plot([], [], color=C_ORANGE, lw=1.8, label="DeltaX = X - X0")
    l_u, = a_act.plot([], [], color=C_BLUE, lw=1.8, label="u (motore)")
    l_dL, = a_act.plot([], [], color=C_ORANGE, lw=1.8, ls="--", label="DeltaLs (molla)")
    l_T, = a_T.plot([], [], color=C_DAHL, lw=1.8, label="T")
    a_F.axhspan(10, 50, color=C_BAND, alpha=0.18, label="Range di interesse 10-50 mN")
    a_F.axhline(F_TH_mN, color=C_ORANGE, ls="--", lw=1.2, label=f"Soglia commutazione {F_TH_mN:.0f} mN")
    l_F, = a_F.plot([], [], color=C_HOOKE, lw=1.6, label="F contatto")
    l_k, = a_k.plot([], [], color=C_PURPLE, lw=1.8, drawstyle="steps-post", label="k")
    for ax in (a_tip, a_act, a_T, a_F, a_k):
        ax.legend(loc="upper left", ncol=3, **LEGEND)
        if ax is not a_k:
            plt.setp(ax.get_xticklabels(), visible=False)

    status = fig.text(0.5, 0.985, "", ha="center", va="top", color=TXT, fontsize=9.5,
                      fontweight="bold")
    fig.text(0.005, 0.005, "Tasto R (in questa finestra): reset loop ed energia",
             ha="left", va="bottom", color=TXT_2, fontsize=8)
    fig.subplots_adjust(left=0.05, right=0.99, top=0.925, bottom=0.075)

    def on_key(event):
        if event.key in ("r", "R"):
            reset_loops()
            print("[LIVE PLOT] loop ed energia dissipata azzerati (tasto R)")
    fig.canvas.mpl_connect("key_press_event", on_key)

    def autoscale(ax, xs, ys, pad=0.08, x0=None):
        if len(xs) < 2:
            return
        x_lo, x_hi = (min(xs), max(xs)) if x0 is None else x0
        y_lo, y_hi = min(ys), max(ys)
        dy = (y_hi - y_lo) or max(abs(y_hi), 1e-3)
        dx = (x_hi - x_lo) or 1.0
        ax.set_xlim(x_lo - (pad * dx if x0 is None else 0), x_hi + (pad * dx if x0 is None else 0))
        ax.set_ylim(y_lo - pad * dy, y_hi + pad * dy)

    def update(_frame):
        with lock:
            t = list(data["t"])
            ts = {k: list(v) for k, v in data.items()}
            lp = {k: list(v) for k, v in loop.items()}
            st = dict(state)
        if not t:
            return ()
        k_now = ts["k"][-1]

        # 1. loop molla
        dls, T = lp["dLs"], lp["T"]
        if dls:
            x_ref = np.linspace(0.0, max(max(dls) * 1.1, 0.5), 50)
            l_hooke.set_data(x_ref, k_now * x_ref)
            l_bup.set_data(x_ref, k_now * x_ref + F_HYST_TOT)
            l_blo.set_data(x_ref, np.maximum(0.0, k_now * x_ref - F_HYST_TOT))
            l_dahl.set_data(dls, T)
            p_sp.set_data([dls[-1]], [T[-1]])
            a_sp.set_xlim(0.0, x_ref[-1])
            a_sp.set_ylim(0.0, max(max(T), k_now * x_ref[-1]) * 1.1 + 0.1)
        t_energy.set_text(f"Energia dissipata: {st['E_diss']:.2f} mJ\n"
                          f"Dahl z: {st['z']:+.2f}   k = {k_now:.2f} N/mm")

        # 2. loop punta e H
        u, dY, br = np.array(lp["u"]), np.array(lp["dY"]), np.array(lp["branch"])
        if u.size:
            ld, ul = br == "load", br == "unload"
            l_load.set_data(u[ld], dY[ld])
            l_unl.set_data(u[ul], dY[ul])
            p_tl.set_data([u[-1]], [dY[-1]])
            autoscale(a_tl, list(u), list(dY))
            try:
                h = hysteresis_H(u, dY, br)
                ug = h["u_at_gap_max_mm"]
                yl = np.interp(ug, *_sorted(u[ld], dY[ld]))
                yu = np.interp(ug, *_sorted(u[ul], dY[ul]))
                l_gap.set_data([ug, ug], [yl, yu])
                t_H.set_text(f"H = {h['H_pct']:.1f} % di DeltaYmax\n"
                             f"gap max = {h['gap_max_mm']*1e3:.1f} um a u = {ug:.2f} mm\n"
                             f"DeltaYmax = {h['dY_max_mm']:.3f} mm")
            except ValueError:
                l_gap.set_data([], [])
                t_H.set_text("H: serve un ciclo carico + scarico")

        # 3-7. dinamica temporale con finestra mobile (come il Modulo 06)
        for line, key in ((l_dY, "dY"), (l_dX, "dX"), (l_u, "u"), (l_dL, "dLs"),
                          (l_T, "T"), (l_F, "F"), (l_k, "k")):
            line.set_data(t, ts[key])
        t_lo = max(0.0, t[-1] - TIME_WINDOW_S)
        t_hi = max(TIME_WINDOW_S, t[-1] + 0.2)
        i0 = next((i for i, tv in enumerate(t) if tv >= t_lo), 0)
        autoscale(a_tip, t[i0:], ts["dY"][i0:] + ts["dX"][i0:], x0=(t_lo, t_hi))
        autoscale(a_act, t[i0:], ts["u"][i0:] + ts["dLs"][i0:], x0=(t_lo, t_hi))
        autoscale(a_T, t[i0:], ts["T"][i0:], x0=(t_lo, t_hi))
        a_F.set_ylim(0, max(55.0, max(ts["F"][i0:]) * 1.1))
        a_k.set_ylim(0, max(ts["k"]) * 1.25)
        a_k.set_xlim(t_lo, t_hi)

        phase = {"load": "CARICO", "unload": "SCARICO", "hold": "MANTENIMENTO",
                 "done": "FINE CICLO"}.get(st["phase"], "MANUALE - " + {"load": "CARICO", "unload": "SCARICO"}[st["branch"]])
        status.set_text(f"t = {t[-1]:.2f} s  |  {phase}  |  u = {ts['u'][-1]:.2f} mm  |  "
                        f"T = {ts['T'][-1]:.2f} N  |  DeltaLs = {ts['dLs'][-1]:.2f} mm  |  "
                        f"punta X = {ts['X'][-1]:.3f}, Y = {ts['Y'][-1]:.4f}, DeltaY = {ts['dY'][-1]:.4f} mm  |  "
                        f"F = {ts['F'][-1]:.1f} mN"
                        + ("   [SOFA chiuso]" if st["closed"] else ""))
        return ()

    if len(sys.argv) > 2 and sys.argv[1] == "--snapshot":   # prova senza finestra
        reader()
        update(0)
        fig.savefig(sys.argv[2], dpi=110, facecolor=BG_FIG)
        return
    threading.Thread(target=reader, daemon=True).start()
    _anim = FuncAnimation(fig, update, interval=200, cache_frame_data=False)  # noqa: F841
    plt.show()


def _sorted(x, y):
    order = np.argsort(x, kind="stable")
    return x[order], y[order]


if __name__ == "__main__":
    main()
