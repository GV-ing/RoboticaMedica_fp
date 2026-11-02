# modules/monte_carlo.py
"""Campionamento, elenco dei casi e analisi statistica del Monte Carlo (Progetto 5).

Traccia: per le configurazioni (a) molla fissa e (b) a due livelli, 20 realizzazioni dei
parametri incerti (modulo uniforme 10-20%, segno casuale). Per ognuna: H (ciclo libero) e
stabilita' del contatto (rigido e soffice). Regressione lineare di H sulle perturbazioni
normalizzate -> diagramma tornado. Si riportano media, deviazione standard e caso peggiore di
H, indice di robustezza (peggiore / nominale) e numero di realizzazioni con contatto instabile.
Target: due parametri dominanti; per la configurazione piu' robusta peggiore <= 1.3 nominale
e contatto stabile in tutte le 20 realizzazioni.
"""

import numpy as np

from modules.config import MC_NUM_REALIZATIONS, MC_UNCERTAINTY_RANGE, NOMINAL_PARAMS, UNCERTAIN_GROUPS

SPRING_CONFIGS = ("fixed", "switched")
TESTS = ("free", "rigid", "soft")
NOMINAL = "nominal"
ROBUSTNESS_TARGET = 1.3


# ------------------------------------------------------------------ campionamento
def sample_realizations(seed=2026, count=MC_NUM_REALIZATIONS):
    """Realizzazioni riproducibili: un delta per gruppo, |delta| ~ U(10%, 20%), segno casuale."""
    if count < 1:
        raise ValueError("count deve essere almeno 1")
    low, high = MC_UNCERTAINTY_RANGE
    rng = np.random.default_rng(seed)
    out = []
    for i in range(count):
        deltas, overrides = {}, {}
        for group, keys in UNCERTAIN_GROUPS.items():
            d = float(rng.uniform(low, high) * rng.choice((-1.0, 1.0)))
            deltas[group] = d
            for key in keys:
                overrides[key] = NOMINAL_PARAMS[key] * (1.0 + d)
        out.append({"realization": f"r{i:02d}", "deltas": deltas, "overrides": overrides})
    return out


def build_tasks(realizations, series_dir, include_direct=True):
    """Elenco dei casi: (nominale + realizzazioni) x configurazione x prova, + diretta nominale."""
    items = [{"realization": NOMINAL, "deltas": {g: 0.0 for g in UNCERTAIN_GROUPS},
              "overrides": {}}] + list(realizations)
    tasks = []
    for r in items:
        for cfg in SPRING_CONFIGS:
            for test in TESTS:
                tasks.append({"case_id": f"{r['realization']}_{cfg}_{test}",
                              "realization": r["realization"], "config": cfg, "test": test,
                              "overrides": r["overrides"], "series_dir": series_dir})
    if include_direct:
        tasks.append({"case_id": f"{NOMINAL}_direct_free", "realization": NOMINAL,
                      "config": "direct", "test": "free", "overrides": {},
                      "series_dir": series_dir})
    cost = {"soft": 3, "rigid": 2, "free": 1}  # i piu' lunghi per primi: carico bilanciato
    return sorted(tasks, key=lambda t: -cost[t["test"]])


# ------------------------------------------------------------------ analisi
def regression_tornado(deltas, H):
    """Regressione lineare H = b0 + sum_g b_g * x_g, con x_g = delta_g / 20% in [-1, 1].

    b_g = variazione di H [punti %] per perturbazione massima (+20%) del parametro g.
    Restituisce i gruppi ordinati per |b_g| (tornado), con coefficiente standardizzato e R^2.
    """
    names = list(UNCERTAIN_GROUPS)
    X = np.array([[d[g] / MC_UNCERTAINTY_RANGE[1] for g in names] for d in deltas], float)
    y = np.asarray(H, float)
    if X.shape[0] <= X.shape[1] + 1:
        raise ValueError("Campioni insufficienti per la regressione")
    A = np.column_stack([np.ones(len(y)), X])
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ beta
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else float("nan")
    sx, sy = X.std(axis=0, ddof=1), y.std(ddof=1)
    rows = [{"parameter": g, "coef_pct_per_20pct": float(b),
             "standardized": float(b * s / sy) if sy > 0 else float("nan")}
            for g, b, s in zip(names, beta[1:], sx)]
    rows.sort(key=lambda r: abs(r["coef_pct_per_20pct"]), reverse=True)
    return {"intercept": float(beta[0]), "r2": r2, "coefficients": rows}


def analyze(results, realizations):
    """Statistiche per configurazione e verifica dei target a partire dai risultati dei casi."""
    ok = {(r["realization"], r["config"], r["test"]): r["metrics"]
          for r in results if "metrics" in r}
    failed = [r["case_id"] for r in results if "error" in r]
    summary = {"failed_cases": failed, "configs": {}}
    H_direct = ok.get((NOMINAL, "direct", "free"), {}).get("H_pct")
    summary["H_direct_nominal_pct"] = H_direct
    for cfg in SPRING_CONFIGS:
        nominal = ok.get((NOMINAL, cfg, "free"), {}).get("H_pct")
        H, D, unstable, out_of_range, used = [], [], {"rigid": 0, "soft": 0}, {"rigid": 0, "soft": 0}, []
        for r in realizations:
            key = (r["realization"], cfg)
            if (key + ("free",)) in ok:
                H.append(ok[key + ("free",)]["H_pct"])
                D.append(r["deltas"])
                used.append(r["realization"])
            for obj in ("rigid", "soft"):
                m = ok.get(key + (obj,))
                if m is None:
                    continue
                unstable[obj] += 0 if m["stable"] else 1
                out_of_range[obj] += 0 if m["in_range_10_50mN"] else 1
        c = {"n_realizations": len(H), "H_nominal_pct": nominal}
        if H:
            Ha = np.array(H)
            c.update(H_mean_pct=float(Ha.mean()),
                     H_std_pct=float(Ha.std(ddof=1)) if Ha.size > 1 else 0.0,
                     H_worst_pct=float(Ha.max()),
                     worst_realization=used[int(Ha.argmax())])
            if nominal:
                c["robustness_index"] = c["H_worst_pct"] / nominal
            if H_direct:
                c["H_nominal_over_direct"] = nominal / H_direct if nominal else None
            if len(H) > len(UNCERTAIN_GROUPS) + 1:
                c["tornado"] = regression_tornado(D, H)
                c["dominant_parameters"] = [row["parameter"]
                                            for row in c["tornado"]["coefficients"][:2]]
        c["unstable_contact"] = unstable
        c["unstable_total"] = unstable["rigid"] + unstable["soft"]
        c["contact_out_of_range_10_50mN"] = out_of_range
        n_contact = sum(1 for r in realizations for o in ("rigid", "soft")
                        if (r["realization"], cfg, o) in ok)
        c["n_contact_tests"] = n_contact
        c["meets_target"] = bool(
            c.get("robustness_index") is not None
            and c["robustness_index"] <= ROBUSTNESS_TARGET
            and c["unstable_total"] == 0
            and len(H) == len(realizations)
            and n_contact == 2 * len(realizations))
        summary["configs"][cfg] = c

    candidates = [cfg for cfg, c in summary["configs"].items() if "robustness_index" in c]
    if candidates:
        # prima chi rispetta il target, poi meno casi instabili, poi indice di robustezza minore
        best = min(candidates, key=lambda k: (not summary["configs"][k]["meets_target"],
                                              summary["configs"][k]["unstable_total"],
                                              summary["configs"][k]["robustness_index"]))
        summary["most_robust_config"] = best
    return summary