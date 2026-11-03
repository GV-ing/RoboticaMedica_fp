#!/usr/bin/env python3
# tools/run_monte_carlo.py  -  SOFA v26.06 / SofaPython3 / Python 3.12
"""
Campagna Monte Carlo del Progetto 5, eseguita in parallelo (un processo SOFA per caso).

Casi: (nominale + N realizzazioni) x {fixed, switched} x {free, rigid, soft} + diretta nominale.
Con N = 20: 127 casi. Ogni caso completato viene aggiunto subito a cases.jsonl: se la campagna
si interrompe, rilanciando lo stesso comando con --resume riparte dai casi mancanti.

Output in --outdir:
  realizations.json   perturbazioni campionate (seed riproducibile)
  cases.jsonl         un risultato per caso (metriche o errore)
  series/<caso>.csv   cronologie complete (per le figure)
  summary.json        statistiche, tornado, verifica dei target

Esecuzione (cartella SOFA_Project, source env.sh):
    python3.12 -u tools/run_monte_carlo.py --count 1 --tests free      # prova rapida (~5 min)
    python3.12 -u tools/run_monte_carlo.py                              # campagna completa
    python3.12 -u tools/run_monte_carlo.py --resume --outdir results/mc_2026   # riprende
"""

import argparse
import datetime as dt
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)

from modules.config import MC_NUM_REALIZATIONS, NOMINAL_PARAMS  # noqa: E402
from modules.mc_worker import run_task  # noqa: E402
from modules.monte_carlo import analyze, build_tasks, sample_realizations  # noqa: E402


def load_done(path):
    done = []
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    done.append(json.loads(line))
                except json.JSONDecodeError:  # riga troncata (crash durante la scrittura)
                    print(f"  riga non valida in {path} ignorata: il caso verra' rieseguito")
    return done


def print_summary(s):
    print("\n=== Riepilogo Monte Carlo ===")
    if s.get("H_direct_nominal_pct") is not None:
        print(f"  H attuazione diretta (nominale): {s['H_direct_nominal_pct']:.2f} %")
    for cfg, c in s["configs"].items():
        print(f"\n  [{cfg}]  realizzazioni con H: {c['n_realizations']}")
        if "H_mean_pct" in c:
            print(f"    H nominale {c['H_nominal_pct']:.2f} %  media {c['H_mean_pct']:.2f} %  "
                  f"dev.std {c['H_std_pct']:.2f} %  peggiore {c['H_worst_pct']:.2f} % "
                  f"({c['worst_realization']})")
            if "robustness_index" in c:
                print(f"    indice di robustezza (peggiore/nominale) = {c['robustness_index']:.3f}"
                      f"  (target <= 1.3)")
            if c.get("H_nominal_over_direct") is not None:
                print(f"    H nominale / H diretta = {c['H_nominal_over_direct']:.3f}")
        print(f"    contatto instabile: rigido {c['unstable_contact']['rigid']}, soffice "
              f"{c['unstable_contact']['soft']} (su {c['n_contact_tests']} prove)")
        print(f"    forza fuori range 10-50 mN: rigido {c['contact_out_of_range_10_50mN']['rigid']}, "
              f"soffice {c['contact_out_of_range_10_50mN']['soft']}")
        if "tornado" in c:
            t = c["tornado"]
            print(f"    tornado (punti di H per +20%), R^2 = {t['r2']:.3f}:")
            for row in t["coefficients"]:
                print(f"      {row['parameter']:10s} {row['coef_pct_per_20pct']:+7.3f}   "
                      f"(standardizzato {row['standardized']:+.3f})")
            print(f"    parametri dominanti: {', '.join(c['dominant_parameters'])}")
        print(f"    target rispettato: {c['meets_target']}")
    if "most_robust_config" in s:
        print(f"\n  Configurazione piu' robusta: {s['most_robust_config']}")
    if s["failed_cases"]:
        print(f"\n  ATTENZIONE: {len(s['failed_cases'])} casi falliti: {s['failed_cases']}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--count", type=int, default=MC_NUM_REALIZATIONS)
    ap.add_argument("--seed", type=int, default=NOMINAL_PARAMS["random_seed"])
    ap.add_argument("--workers", type=int, default=NOMINAL_PARAMS["mc_workers"])
    ap.add_argument("--tests", nargs="+", default=["free", "rigid", "soft"],
                    choices=["free", "rigid", "soft"])
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--configs", nargs="+", default=["fixed", "switched", "direct"],
                    choices=["fixed", "switched", "direct"])
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--rerun", action="store_true",
                    help="Con --resume: riesegue i casi selezionati (--configs/--tests) anche se "
                         "gia' fatti; nel riepilogo il risultato nuovo sostituisce il vecchio")
    ap.add_argument("--tag", default="",
                    help="Etichetta salvata nei risultati (es. 'k_istantaneo')")
    ap.add_argument("--dry-run", action="store_true",
                    help="Verifica la pipeline con metriche sintetiche, senza SOFA")
    a = ap.parse_args()

    outdir = a.outdir or os.path.join(
        PROJECT_DIR, "results", "mc_" + dt.datetime.now().strftime("%Y%m%d_%H%M"))
    os.makedirs(outdir, exist_ok=True)
    real_path = os.path.join(outdir, "realizations.json")
    cases_path = os.path.join(outdir, "cases.jsonl")

    if not a.resume and os.path.exists(cases_path):
        sys.exit(f"{cases_path} esiste gia': usare --resume o un altro --outdir")
    realizations = sample_realizations(seed=a.seed, count=a.count)
    if a.resume and os.path.exists(real_path):
        with open(real_path, encoding="utf-8") as f:
            saved = json.load(f)
        if saved["seed"] != a.seed or saved["count"] != a.count:
            sys.exit("Seed o count diversi dalla campagna da riprendere")
    with open(real_path, "w", encoding="utf-8") as f:
        json.dump({"seed": a.seed, "count": a.count, "realizations": realizations}, f, indent=2)

    tasks = [t for t in build_tasks(realizations, os.path.join(outdir, "series"))
             if t["test"] in a.tests and t["config"] in a.configs]
    for t in tasks:
        t["dry_run"] = a.dry_run
        t["tag"] = a.tag
    done = load_done(cases_path) if a.resume else []
    if a.rerun and not a.resume:
        sys.exit("--rerun richiede --resume (stessa cartella della campagna)")
    selected = {t["case_id"] for t in tasks}
    done_ids = {r["case_id"] for r in done if "metrics" in r}
    if a.rerun:
        # da rifare: i casi selezionati non ancora rieseguiti con questa etichetta
        # (cosi' una riesecuzione interrotta riprende senza ripetere i casi gia' rifatti)
        redone = {r["case_id"] for r in done
                  if "metrics" in r and a.tag and r.get("tag") == a.tag}
        done_ids = (done_ids - selected) | (redone & selected)
    todo = [t for t in tasks if t["case_id"] not in done_ids]
    print(f"Campagna in {outdir}\n  casi totali {len(tasks)}, gia' fatti {len(done_ids)}, "
          f"da eseguire {len(todo)}, processi {a.workers}")

    results = [r for r in done if "metrics" in r]  # in ordine: il piu' recente prevale
    t0 = time.time()
    if todo:
        ctx = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx,
                                 max_tasks_per_child=1) as pool, \
                open(cases_path, "a", encoding="utf-8") as log:
            futures = {pool.submit(run_task, t): t for t in todo}
            for i, fut in enumerate(as_completed(futures), 1):
                t = futures[fut]
                try:
                    r = fut.result()
                except Exception as exc:  # crash del processo figlio
                    r = {k: t[k] for k in ("case_id", "realization", "config", "test")}
                    r["error"] = f"processo terminato: {exc}"
                log.write(json.dumps(r) + "\n")
                log.flush()
                results.append(r)
                el = time.time() - t0
                eta = el / i * (len(todo) - i)
                status = "ERRORE " + r["error"] if "error" in r else "ok"
                print(f"[{i:3d}/{len(todo)}] {r['case_id']:28s} {status}   "
                      f"trascorso {el / 60:5.1f} min, stima residua {eta / 60:5.1f} min",
                      flush=True)

    summary = analyze(results, realizations)
    summary.update(seed=a.seed, count=a.count, tests=a.tests)
    with open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print_summary(summary)
    print(f"\nRisultati in {outdir}")


if __name__ == "__main__":
    main()