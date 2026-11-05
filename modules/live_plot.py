# modules/live_plot.py
"""Listener del controller usati SOLO nella GUI (runSofa): grafici live e registrazione CSV.

  LivePlotListener  avvia tools/live_plot_viewer.py come processo separato (Python di sistema
                    con matplotlib e tkinter, NON il python3.12 di SOFA) e gli
                    invia un campione ogni `plot_every` passi; se la finestra non parte o viene
                    chiusa, si disattiva senza fermare la simulazione.
  CsvRecorder       scrive ogni passo in results/gui_<data_ora>.csv (export CSV richiesto).

Il viewer parte SENZA il PYTHONPATH di SOFA, cosi' usa il proprio matplotlib/numpy.
"""

import atexit
import csv
import datetime as dt
import json
import os
import shutil
import subprocess

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VIEWER = os.path.join(PROJECT_DIR, "tools", "live_plot_viewer.py")


def _viewer_env():
    """Ambiente pulito per il viewer: niente librerie/numpy di SOFA (compilati per 3.12)."""
    env = dict(os.environ)
    for var in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH"):
        env.pop(var, None)
    env["PYTHONPATH"] = PROJECT_DIR          # modules.config e modules.metrics (puri, numpy)
    return env


def find_viewer_python(preferred=None):
    """Primo interprete che importa matplotlib e tkinter (finestra TkAgg) senza SOFA.

    Ordine: preferito (config 'live_plot_python'), Python di sistema (Ubuntu 22.04: 3.10 con
    python3-matplotlib e python3-tk da apt), poi python3.12.
    """
    candidates = [preferred, "/usr/bin/python3", shutil.which("python3"), shutil.which("python3.12")]
    seen = set()
    for exe in candidates:
        if not exe or exe in seen or not os.path.exists(exe):
            continue
        seen.add(exe)
        try:
            ok = subprocess.run([exe, "-c", "import matplotlib, tkinter"], env=_viewer_env(),
                                capture_output=True, timeout=20).returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            ok = False
        if ok:
            return exe
    return None


class LivePlotListener:
    def __init__(self, every=5, python=None):
        self.every = max(1, int(every))
        self.count = 0
        self.proc = None
        exe = find_viewer_python(python)
        if exe is None:
            print("[LIVE PLOT] disattivati: nessun Python con matplotlib e tkinter. "
                  "Installare con: sudo apt install python3-matplotlib python3-tk")
            return
        try:
            self.proc = subprocess.Popen([exe, "-u", VIEWER], stdin=subprocess.PIPE,
                                         text=True, env=_viewer_env())
            print(f"[LIVE PLOT] finestra avviata con {exe}")
        except OSError as exc:
            print(f"[LIVE PLOT] disattivati: {exc}")

    def __call__(self, sample):
        if self.proc is None:
            return
        self.count += 1
        if self.count % self.every:
            return
        if self.proc.poll() is not None:     # finestra chiusa o viewer fallito
            print("[LIVE PLOT] finestra chiusa: grafici disattivati (la simulazione continua)")
            self.proc = None
            return
        try:
            self.proc.stdin.write(json.dumps(sample) + "\n")
            self.proc.stdin.flush()
        except (BrokenPipeError, OSError):
            self.proc = None


class CsvRecorder:
    def __init__(self, fields, outdir=None):
        outdir = outdir or os.path.join(PROJECT_DIR, "results")
        os.makedirs(outdir, exist_ok=True)
        self.path = os.path.join(outdir, "gui_" + dt.datetime.now().strftime("%Y%m%d_%H%M%S")
                                 + ".csv")
        self.fields = list(fields)
        self.file = open(self.path, "w", newline="", encoding="utf-8")
        self.writer = csv.writer(self.file)
        self.writer.writerow(self.fields)
        self.rows = 0
        atexit.register(self.file.close)     # chiude (e svuota) il file all'uscita di runSofa
        print(f"[CSV] registrazione in {self.path}")

    def __call__(self, sample):
        self.writer.writerow([sample[f] for f in self.fields])
        self.rows += 1
        if self.rows % 100 == 0:
            self.file.flush()
