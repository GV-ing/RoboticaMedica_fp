# modules/metrics.py
"""Metriche del Progetto 5 (pure, testabili senza SOFA).

Isteresi H (traccia + glossario):
    H = max |dY_load(u) - dY_unload(u)| / dY_max * 100   [%]
    confronto a pari u su un ciclo; dY_max = deflessione massima a fine carico.
    Riferimento: attuazione diretta (k = 50 N/mm) -> si riporta anche H / H_diretta.

Stabilita' del contatto (traccia):
    forza che oscilla meno del 30% (picco-picco) rispetto alla media per 5 s,
    senza perdita di contatto.
"""

import numpy as np


def _branch(u, y):
    """Ordina per u e media i valori con u ripetuto (es. fermi del carrello)."""
    order = np.argsort(u, kind="stable")
    u_s, y_s = u[order], y[order]
    uniq, inv = np.unique(u_s, return_inverse=True)
    y_mean = np.bincount(inv, weights=y_s) / np.bincount(inv)
    return uniq, y_mean


def hysteresis_H(u, dY, phase, n_grid=400):
    """Indice di isteresi H [%] da un ciclo con fasi 'load' / 'hold' / 'unload'.

    Restituisce un dizionario con H_pct, dY_max, u del massimo scarto e i rami interpolati.
    """
    u = np.asarray(u, dtype=float)
    dY = np.asarray(dY, dtype=float)
    phase = np.asarray(phase)
    if not (u.shape == dY.shape == phase.shape) or u.ndim != 1:
        raise ValueError("u, dY e phase devono essere vettori della stessa lunghezza")
    if not (np.all(np.isfinite(u)) and np.all(np.isfinite(dY))):
        raise ValueError("u o dY contengono valori non finiti")
    load = phase == "load"
    unload = phase == "unload"
    if load.sum() < 3 or unload.sum() < 3:
        raise ValueError("Servono almeno 3 campioni sia in carico sia in scarico")

    uL, yL = _branch(u[load], dY[load])
    uU, yU = _branch(u[unload], dY[unload])
    lo, hi = max(uL[0], uU[0]), min(uL[-1], uU[-1])
    if not hi > lo:
        raise ValueError("I rami di carico e scarico non si sovrappongono in u")
    grid = np.linspace(lo, hi, n_grid)
    gap = np.abs(np.interp(grid, uL, yL) - np.interp(grid, uU, yU))

    dY_max = float(np.max(dY[load | (phase == "hold")]))
    if not dY_max > 0:
        raise ValueError("dY_max non positivo: il ciclo non ha prodotto deflessione")
    i = int(np.argmax(gap))
    return {
        "H_pct": float(100.0 * gap[i] / dY_max),
        "gap_max_mm": float(gap[i]),
        "u_at_gap_max_mm": float(grid[i]),
        "dY_max_mm": dY_max,
        "u_max_mm": float(np.max(u)),
    }


def contact_stability(t, F, phase, duration=5.0, ptp_max=0.30, loss_force=1e-4):
    """Stabilita' del contatto sugli ultimi `duration` secondi della fase 'hold'."""
    t = np.asarray(t, dtype=float)
    F = np.asarray(F, dtype=float)
    phase = np.asarray(phase)
    hold = phase == "hold"
    if hold.sum() < 2:
        raise ValueError("Nessuna fase 'hold' registrata")
    t_h, F_h = t[hold], F[hold]
    span = float(t_h[-1] - t_h[0])
    window = t_h >= t_h[-1] - duration + 1e-9
    Fw = F_h[window]
    mean = float(np.mean(Fw))
    ptp = float(np.ptp(Fw))
    lost = bool(np.any(Fw < loss_force))
    complete = span >= duration - 2 * np.min(np.diff(t_h)) if t_h.size > 1 else False
    ratio = ptp / mean if mean > 0 else float("inf")
    return {
        "F_mean_N": mean,
        "F_ptp_N": ptp,
        "ptp_over_mean": float(ratio),
        "contact_lost": lost,
        "hold_duration_s": span,
        "window_complete": bool(complete),
        "in_range_10_50mN": bool(0.010 <= mean <= 0.050),
        "stable": bool(complete and not lost and ratio < ptp_max),
    }
