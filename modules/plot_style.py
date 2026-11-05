# modules/plot_style.py
"""Stile comune per grafici live e figure della relazione (matplotlib, senza SOFA).

Palette categoriale validata (primi tre slot sicuri anche per daltonismo, a coppie):
blu, arancio, acqua. Grigi per assi, griglia e linee di riferimento.
"""

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
REF = "#8a8984"          # linee di riferimento (soglie, target)
BAND = "#f0efec"         # bande di riferimento (range 10-50 mN)

SERIES = ["#2a78d6", "#eb6834", "#1baf7a"]   # slot 1-3
CONFIG_COLOR = {"fixed": SERIES[0], "switched": SERIES[1], "direct": SERIES[2]}
CONFIG_LABEL = {"fixed": "(a) molla fissa 1.2 N/mm",
                "switched": "(b) molla a due livelli 0.6/4.8 N/mm",
                "direct": "attuazione diretta 50 N/mm"}


def apply(plt):
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": GRID, "axes.labelcolor": TEXT_2, "axes.titlecolor": TEXT,
        "axes.titlesize": 11, "axes.titleweight": "bold", "axes.labelsize": 9,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": TEXT_2, "ytick.color": TEXT_2,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "legend.frameon": False, "legend.labelcolor": TEXT,
        "lines.linewidth": 2.0, "lines.solid_capstyle": "round",
        "font.size": 9, "text.color": TEXT,
    })
