# modules/mc_worker.py
"""Esecuzione di un singolo caso Monte Carlo in un processo separato (spawn).

Ogni caso gira in un processo nuovo (max_tasks_per_child = 1): SOFA parte pulito e la
memoria viene liberata alla fine del caso.
"""

import csv
import os
import traceback

HISTORY_KEYS = ("t", "phase", "u", "u_target", "u1", "dLs", "T", "k", "z", "E",
                "tip_x", "tip_y", "tip_z", "dY", "F_contact")


def _save_series(path, history):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(HISTORY_KEYS)
        w.writerows(zip(*(history[k] for k in HISTORY_KEYS)))


def run_task(task):
    """task: dict con case_id, realization, config, test, overrides, series_dir."""
    from modules.protocols import CONFIGS, contact_test, free_cycle

    result = {k: task[k] for k in ("case_id", "realization", "config", "test")}
    if task.get("dry_run"):  # verifica della pipeline senza SOFA (metriche sintetiche)
        import math
        seed = sum(map(ord, task["case_id"]))
        if task["test"] == "free":
            yt = task["overrides"].get("yt_offset", 0.6) / 0.6 - 1
            result["metrics"] = {"H_pct": 15.0 * (1 + 0.8 * yt) + 0.1 * math.sin(seed)}
        else:
            result["metrics"] = {"stable": True, "in_range_10_50mN": task["test"] == "soft"}
        return result
    try:
        overrides = {**CONFIGS[task["config"]], **task["overrides"]}
        if task["test"] == "free":
            metrics, history = free_cycle(overrides)
        else:
            metrics, history = contact_test(task["test"], overrides)
        _save_series(os.path.join(task["series_dir"], f"{task['case_id']}.csv"), history)
        result["metrics"] = metrics
    except Exception as exc:  # il caso fallito viene registrato, la campagna prosegue
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["traceback"] = traceback.format_exc()
    return result
