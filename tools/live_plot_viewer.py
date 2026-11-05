#!/usr/bin/env python3
# tools/live_plot_viewer.py
"""Finestra dei grafici in tempo reale (processo separato da SOFA).

Riceve dal controller (modules/live_plot.py) un campione JSON per riga su stdin e aggiorna:
  1. posizione della punta: dY (= Y - Y0) e dX = X - X0  [mm]
  2. attuazione: u (motore) e dLs (allungamento molla)    [mm]
  3. tensione del tendine T                              [N]
  4. forza di contatto, con soglia 5 mN e range 10-50 mN [mN]
  5. rigidezza della molla k (commutazione istantanea 0.6 <-> 4.8 N/mm) [N/mm]
Una grandezza per asse y (niente doppi assi): le unita' diverse stanno in pannelli diversi.
Si chiude da solo quando SOFA chiude lo stdin.

Gira con il Python di sistema (Ubuntu 22.04: python3 3.10 + python3-matplotlib + python3-tk),
non con il python3.12 di SOFA. Avvio manuale per prova:
    env -u PYTHONPATH python3 tools/live_plot_viewer.py < campioni.jsonl
Immagine senza finestra:
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
from matplotlib.animation import FuncAnimation  # noqa: E402

from modules import plot_style as st  # noqa: E402

MAX_POINTS = 20000
data = {k: deque(maxlen=MAX_POINTS) for k in ("t", "dX", "dY", "u", "dLs", "T", "F", "k")}
state = {"x0": None, "phase": "", "closed": False}
lock = threading.Lock()


def reader():
    for line in sys.stdin:
        try:
            s = json.loads(line)
        except json.JSONDecodeError:
            continue
        with lock:
            if state["x0"] is None:
                state["x0"] = s["tip_x"]
            data["t"].append(s["t"])
            data["dX"].append(s["tip_x"] - state["x0"])
            data["dY"].append(s["dY"])
            data["u"].append(s["u"])
            data["dLs"].append(s["dLs"])
            data["T"].append(s["T"])
            data["F"].append(s["F_contact"] * 1e3)
            data["k"].append(s["k"])
            state["phase"] = s.get("phase", "")
    state["closed"] = True


def main():
    st.apply(plt)
    fig, axes = plt.subplots(5, 1, sharex=True, figsize=(8.5, 11),
                             gridspec_kw={"height_ratios": [1, 1, 1, 1, 0.6]})
    fig.canvas.manager.set_window_title("Progetto 5 - grafici in tempo reale")
    a_tip, a_act, a_T, a_F, a_k = axes

    l_dY, = a_tip.plot([], [], color=st.SERIES[0], label="dY = Y - Y0")
    l_dX, = a_tip.plot([], [], color=st.SERIES[1], label="dX = X - X0")
    a_tip.set_ylabel("punta [mm]")
    a_tip.set_title("Posizione della punta", loc="left")

    l_u, = a_act.plot([], [], color=st.SERIES[0], label="u (motore)")
    l_dL, = a_act.plot([], [], color=st.SERIES[1], label="dLs (molla)")
    a_act.set_ylabel("corsa [mm]")
    a_act.set_title("Attuazione", loc="left")

    l_T, = a_T.plot([], [], color=st.SERIES[0], label="T")
    a_T.set_ylabel("tensione [N]")
    a_T.set_title("Tensione del tendine", loc="left")

    a_F.axhspan(10, 50, color=st.BAND, zorder=0, label="range 10-50 mN")
    a_F.axhline(5, color=st.REF, lw=1, ls="--", label="soglia commutazione 5 mN")
    l_F, = a_F.plot([], [], color=st.SERIES[0], label="F contatto")
    a_F.set_ylabel("forza [mN]")
    a_F.set_title("Forza di contatto", loc="left")

    l_k, = a_k.plot([], [], color=st.SERIES[0], label="k", drawstyle="steps-post")
    a_k.set_ylabel("k [N/mm]")
    a_k.set_xlabel("tempo simulato [s]")
    a_k.set_title("Rigidezza della molla", loc="left")

    for ax in axes:  # legenda sopra il pannello, a destra del titolo: non copre le curve
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3, borderaxespad=0.2)
    status = fig.text(0.99, 0.995, "", ha="right", va="top", color=st.TEXT_2, fontsize=8)
    fig.tight_layout(rect=(0, 0, 1, 0.98))

    def update(_frame):
        with lock:
            t = list(data["t"])
            series = {k: list(v) for k, v in data.items()}
            phase = state["phase"]
        if not t:
            return ()
        for line, key in ((l_dY, "dY"), (l_dX, "dX"), (l_u, "u"), (l_dL, "dLs"),
                          (l_T, "T"), (l_F, "F"), (l_k, "k")):
            line.set_data(t, series[key])
        for ax in axes:
            ax.relim()
            ax.autoscale_view()
        a_F.set_ylim(bottom=0, top=max(55.0, max(series["F"]) * 1.1))
        a_k.set_ylim(0, max(series["k"]) * 1.2)
        status.set_text(f"t = {t[-1]:.2f} s | fase: {phase} | k = {series['k'][-1]:.2f} N/mm"
                        f" | T = {series['T'][-1]:.2f} N | dY = {series['dY'][-1]:.4f} mm"
                        + ("  [SOFA chiuso]" if state["closed"] else ""))
        return ()

    if len(sys.argv) > 2 and sys.argv[1] == "--snapshot":   # prova senza finestra
        reader()
        update(0)
        fig.savefig(sys.argv[2], dpi=110)
        return
    threading.Thread(target=reader, daemon=True).start()
    _anim = FuncAnimation(fig, update, interval=200, cache_frame_data=False)  # noqa: F841
    plt.show()


if __name__ == "__main__":
    main()
