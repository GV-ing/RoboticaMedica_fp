#!/usr/bin/env python3
# tools/make_figures.py  -  figure della relazione dai risultati del Monte Carlo
"""
Genera le figure del Progetto 5 da una cartella di campagna (results/mc_finale):

  fig_tornado.png          coefficienti di regressione di H (punti di H per +20%)
  fig_H_realizzazioni.png  H per realizzazione, nominale e soglia 1.3 x nominale
  fig_forze_contatto.png   forza di contatto nella finestra di 5 s, rigido e soffice
  fig_isteresi_nominale.png   cicli dY(u) nominali    (se ci sono le serie temporali)
  fig_contatto_nominale.png   F(t) e k(t) nominali    (se ci sono le serie temporali)

Serie temporali: <campagna>/series/<caso>.csv; in mancanza, per i cicli liberi si usano
results/free_<config>.csv prodotti da tools/run_cycle.py.

Uso (non richiede SOFA):
    python3.12 tools/make_figures.py --campaign results/mc_finale --out results/figures
"""

import argparse
import csv
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Palette categoriale validata (slot 1-2) + inchiostri
C_FIXED, C_SWITCH, C_DIRECT = "#2a78d6", "#eb6834", "#52514e"
INK, INK2, GRID, BAND = "#0b0b0b", "#52514e", "#e4e3df", "#1baf7a"
LABEL = {"fixed": "(a) molla fissa k = 1.2 N/mm", "switched": "(b) due livelli 0.6 / 4.8 N/mm",
         "direct": "attuazione diretta k = 50 N/mm"}
COLOR = {"fixed": C_FIXED, "switched": C_SWITCH, "direct": C_DIRECT}
PARAM = {"E": "E (carico/scarico)", "Fh_mat": "F_h materiale", "Fh_spring": "F_h molla",
         "c": "smorzamento c", "k": "rigidezza k", "yt": "offset tendine y_t"}


def style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 200, "font.size": 9.5,
        "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.titlesize": 10.5,
        "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "xtick.color": INK2, "ytick.color": INK2, "legend.frameon": False,
        "lines.linewidth": 2.0,
    })


def load(campaign):
    with open(os.path.join(campaign, "summary.json"), encoding="utf-8") as f:
        summary = json.load(f)
    with open(os.path.join(campaign, "realizations.json"), encoding="utf-8") as f:
        reals = json.load(f)["realizations"]
    latest = {}
    with open(os.path.join(campaign, "cases.jsonl"), encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                latest[r["case_id"]] = r          # l'ultima riga del caso prevale
    return summary, reals, latest


def metric(latest, real, cfg, test, key):
    r = latest.get(f"{real}_{cfg}_{test}")
    return r["metrics"].get(key) if r and "metrics" in r else None


def fig_tornado(summary, out):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4), sharex=True)
    for ax, cfg in zip(axes, ("fixed", "switched")):
        t = summary["configs"][cfg]["tornado"]
        rows = sorted(t["coefficients"], key=lambda r: abs(r["coef_pct_per_20pct"]))
        names = [PARAM[r["parameter"]] for r in rows]
        vals = [r["coef_pct_per_20pct"] for r in rows]
        ax.barh(names, vals, color=COLOR[cfg], height=0.55, edgecolor="white", linewidth=1)
        ax.axvline(0, color=INK2, linewidth=1)
        for i, v in enumerate(vals):
            if abs(v) >= 0.3:
                ax.text(v + (0.06 if v > 0 else -0.06), i, f"{v:+.2f}", va="center",
                        ha="left" if v > 0 else "right", color=INK, fontsize=8.5)
        ax.set_title(f"{LABEL[cfg]}\nR² = {t['r2']:.3f}", loc="left")
        ax.grid(axis="y", visible=False)
        ax.set_xlabel("variazione di H [punti %] per +20% del parametro")
    lim = max(abs(r["coef_pct_per_20pct"]) for c in ("fixed", "switched")
              for r in summary["configs"][c]["tornado"]["coefficients"]) * 1.35
    axes[0].set_xlim(-lim, lim)
    fig.suptitle("Diagramma tornado: regressione lineare di H sulle perturbazioni normalizzate",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_tornado.png"))
    plt.close(fig)


def fig_H(summary, reals, latest, out):
    ids = [r["realization"] for r in reals]
    x = np.arange(len(ids))
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    for cfg, dx in (("fixed", -0.12), ("switched", 0.12)):
        c = summary["configs"][cfg]
        H = [metric(latest, i, cfg, "free", "H_pct") for i in ids]
        nom = c["H_nominal_pct"]
        ax.axhline(nom, color=COLOR[cfg], linewidth=1.2, alpha=0.9)
        ax.axhline(1.3 * nom, color=COLOR[cfg], linewidth=1.2, linestyle=(0, (4, 3)), alpha=0.9)
        ax.scatter(x + dx, H, s=34, color=COLOR[cfg], edgecolor="white", linewidth=1, zorder=3,
                   label=f"{LABEL[cfg]}: media {c['H_mean_pct']:.1f} ± {c['H_std_pct']:.1f} %, "
                         f"robustezza {c['robustness_index']:.2f}")
        ax.text(len(ids) - 0.4, 1.3 * nom, f"1.3 × nominale {'(a)' if cfg == 'fixed' else '(b)'}", color=INK2,
                fontsize=8, va="bottom", ha="right")
    Hd = summary.get("H_direct_nominal_pct")
    if Hd:
        ax.axhline(Hd, color=C_DIRECT, linewidth=1.2, linestyle=(0, (1, 2)))
        ax.text(len(ids) - 0.4, Hd, f"attuazione diretta (nominale) {Hd:.1f} %", color=INK2,
                fontsize=8, va="bottom", ha="right")
    ax.set_xticks(x, ids, rotation=0, fontsize=8)
    ax.set_xlim(-0.6, len(ids) - 0.4)
    ax.set_ylabel("H [%]")
    ax.set_title("Isteresi H nelle 20 realizzazioni (linea continua: nominale; tratteggio: target)",
                 loc="left")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=1, fontsize=8.5)
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig_H_realizzazioni.png"), bbox_inches="tight")
    plt.close(fig)


def fig_forces(reals, latest, out):
    ids = [r["realization"] for r in reals]
    x = np.arange(len(ids))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    for ax, obj, title in ((axes[0], "rigid", "Oggetto rigido"),
                           (axes[1], "soft", "Oggetto soffice (E = 5 kPa)")):
        ax.axhspan(10, 50, color=BAND, alpha=0.10, linewidth=0)
        ax.axhline(10, color=BAND, linewidth=1)
        ax.axhline(50, color=BAND, linewidth=1)
        for cfg, dx in (("fixed", -0.12), ("switched", 0.12)):
            F = [1e3 * (metric(latest, i, cfg, obj, "F_mean_N") or np.nan) for i in ids]
            st = [metric(latest, i, cfg, obj, "stable") for i in ids]
            ax.scatter(x + dx, F, s=30, color=COLOR[cfg], edgecolor="white", linewidth=1,
                       zorder=3, label=LABEL[cfg])
            bad = [j for j, s in enumerate(st) if s is False]
            if bad:
                ax.scatter(x[bad] + dx, np.array(F)[bad], s=80, facecolor="none",
                           edgecolor=INK, linewidth=1.2, zorder=4)
            Fn = metric(latest, "nominal", cfg, obj, "F_mean_N")
            if Fn is not None:
                ax.axhline(1e3 * Fn, color=COLOR[cfg], linewidth=1, alpha=0.8)
        ax.text(0.98, 38, "range di interesse\n10–50 mN", transform=ax.get_yaxis_transform(),
                color=INK2, fontsize=8, va="center", ha="right")
        ax.set_xticks(x[::2], ids[::2], fontsize=8)
        ax.set_ylim(0, 65)
        ax.set_ylabel("forza media nella finestra di 5 s [mN]")
        ax.set_title(title, loc="left")
        ax.grid(axis="x", visible=False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.suptitle("Prova di contatto: arresto a 30 mN (o a fine corsa), poi motore fermo; "
                 "tutte le prove risultano stabili", x=0.01, ha="left", fontsize=10.5,
                 fontweight="bold")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=8.5)
    fig.savefig(os.path.join(out, "fig_forze_contatto.png"), bbox_inches="tight")
    plt.close(fig)


def read_series(path):
    with open(path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    out = {}
    for k in rows[0]:
        if k == "phase":
            out[k] = np.array([r[k] for r in rows])
        else:
            out[k] = np.array([float(r[k]) for r in rows])
    return out


def find_series(campaign, case_id, fallback=None):
    p = os.path.join(campaign, "series", f"{case_id}.csv")
    if os.path.exists(p):
        return read_series(p)
    if fallback and os.path.exists(fallback):
        return read_series(fallback)
    return None


def fig_hysteresis(summary, campaign, out):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4))
    found = False
    for ax, cfg in zip(axes, ("fixed", "switched", "direct")):
        s = find_series(campaign, f"nominal_{cfg}_free",
                        os.path.join(PROJECT_DIR, "results", f"free_{cfg}.csv"))
        if s is None:
            ax.set_visible(False)
            continue
        found = True
        load_ = s["phase"] == "load"
        unl = s["phase"] == "unload"
        ax.plot(s["u"][load_], 1e3 * s["dY"][load_], color=COLOR[cfg], label="carico")
        ax.plot(s["u"][unl], 1e3 * s["dY"][unl], color=COLOR[cfg], linestyle=(0, (4, 2)),
                label="scarico")
        H = (summary["H_direct_nominal_pct"] if cfg == "direct"
             else summary["configs"][cfg]["H_nominal_pct"])
        ax.set_title(f"{LABEL[cfg]}\nH = {H:.1f} %", loc="left")
        ax.set_xlabel("comando motore u [mm]")
        ax.set_ylabel("deflessione punta ΔY [µm]")
        ax.legend(loc="upper left", fontsize=8)
    if found:
        fig.suptitle("Cicli di carico e scarico nominali in moto libero", x=0.01, ha="left",
                     fontsize=11, fontweight="bold")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "fig_isteresi_nominale.png"))
    plt.close(fig)
    return found


def fig_contact_time(campaign, out):
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 5.2), sharex="col")
    found = False
    for col, obj in enumerate(("rigid", "soft")):
        for cfg in ("fixed", "switched"):
            s = find_series(campaign, f"nominal_{cfg}_{obj}")
            if s is None:
                continue
            found = True
            axes[0, col].plot(s["t"], 1e3 * s["F_contact"], color=COLOR[cfg], label=LABEL[cfg])
            axes[1, col].plot(s["t"], s["k"], color=COLOR[cfg])
        axes[0, col].axhspan(10, 50, color=BAND, alpha=0.10, linewidth=0)
        axes[0, col].set_title("Oggetto rigido" if obj == "rigid" else "Oggetto soffice",
                               loc="left")
        axes[0, col].set_ylabel("forza di contatto [mN]")
        axes[1, col].set_ylabel("rigidezza molla k [N/mm]")
        axes[1, col].set_xlabel("tempo [s]")
    if found:
        axes[0, 0].legend(loc="lower right", fontsize=8)
        fig.suptitle("Prove di contatto nominali: commutazione istantanea di k a 5 mN",
                     x=0.01, ha="left", fontsize=11, fontweight="bold")
        fig.tight_layout()
        fig.savefig(os.path.join(out, "fig_contatto_nominale.png"))
    plt.close(fig)
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--campaign", default=os.path.join(PROJECT_DIR, "results", "mc_finale"))
    ap.add_argument("--out", default=os.path.join(PROJECT_DIR, "results", "figures"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    style()
    summary, reals, latest = load(a.campaign)
    fig_tornado(summary, a.out)
    fig_H(summary, reals, latest, a.out)
    fig_forces(reals, latest, a.out)
    h = fig_hysteresis(summary, a.campaign, a.out)
    c = fig_contact_time(a.campaign, a.out)
    print(f"Figure in {a.out}" + ("" if h else "  (cicli di isteresi: serie mancanti)")
          + ("" if c else "  (contatto nel tempo: serie mancanti)"))


if __name__ == "__main__":
    main()