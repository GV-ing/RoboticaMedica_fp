"""Monte Carlo orchestration and post-processing for the SOFA tendon robot."""

import json
import gc
import sys
from pathlib import Path
from typing import Callable

import numpy as np

from modules.config import (
    MC_NUM_REALIZATIONS,
    MC_UNCERTAINTY_RANGE,
    NOMINAL_PARAMS,
    UNCERTAIN_PARAMS,
)


def sample_realizations(seed=2026, count=MC_NUM_REALIZATIONS):
    """Return independent parameter sets with signed 10-20% perturbations."""
    if count < 1:
        raise ValueError("count deve essere almeno 1")
    low, high = MC_UNCERTAINTY_RANGE
    if not 0 <= low <= high:
        raise ValueError("Intervallo di incertezza non valido")

    rng = np.random.default_rng(seed)
    realizations = []
    for run_id in range(count):
        params = NOMINAL_PARAMS.copy()
        perturbations = {}

        def draw(key):
            magnitude = rng.uniform(low, high)
            delta = float(magnitude * rng.choice((-1.0, 1.0)))
            params[key] *= 1.0 + delta
            perturbations[key] = delta

        while True:
            params["E_load"] = NOMINAL_PARAMS["E_load"]
            params["E_unload"] = NOMINAL_PARAMS["E_unload"]
            draw("E_load")
            draw("E_unload")
            if params["E_load"] > params["E_unload"]:
                break
            params["E_load"] = NOMINAL_PARAMS["E_load"]
            params["E_unload"] = NOMINAL_PARAMS["E_unload"]

        for key in UNCERTAIN_PARAMS:
            if key not in {"E_load", "E_unload"}:
                draw(key)

        realizations.append({
            "run_id": run_id,
            "params": params,
            "relative_perturbations": perturbations,
        })
    return realizations


def hysteresis_area_Nmm(tip_displacement_mm, tendon_force_N):
    """Compute the closed-cycle force-displacement area in N mm."""
    displacement = np.asarray(tip_displacement_mm, dtype=float)
    force = np.asarray(tendon_force_N, dtype=float)
    if displacement.ndim != 1 or force.ndim != 1:
        raise ValueError("I campioni di isteresi devono essere vettori")
    if displacement.size < 3 or displacement.size != force.size:
        raise ValueError("Storico di isteresi insufficiente o disallineato")
    if not np.all(np.isfinite(displacement)) or not np.all(np.isfinite(force)):
        raise ValueError("Lo storico di isteresi contiene valori non finiti")
    return float(abs(np.sum(
        0.5 * (force[:-1] + force[1:]) * np.diff(displacement)
    )))


def contact_metrics(contact_force_N, threshold_N, target_force_N, tolerance_N):
    """Summarize threshold detection and force stability during contact."""
    force = np.asarray(contact_force_N, dtype=float)
    if force.ndim != 1 or force.size == 0 or not np.all(np.isfinite(force)):
        raise ValueError("Forze di contatto vuote, non vettoriali o non finite")
    if np.any(force < 0):
        raise ValueError("Le forze di contatto devono essere non negative")

    in_contact = force >= threshold_N
    detected_forces = force[in_contact]
    if detected_forces.size == 0:
        return {
            "contact_fraction": 0.0,
            "contact_force_mean_N": None,
            "contact_force_std_N": None,
            "target_force_error_mean_N": None,
            "stable": False,
        }

    mean_force = float(np.mean(detected_forces))
    force_std = float(np.std(detected_forces, ddof=1)) if detected_forces.size > 1 else 0.0
    target_error = float(np.mean(np.abs(detected_forces - target_force_N)))
    return {
        "contact_fraction": float(np.mean(in_contact)),
        "contact_force_mean_N": mean_force,
        "contact_force_std_N": force_std,
        "target_force_error_mean_N": target_error,
        "stable": bool(force_std <= tolerance_N and target_error <= tolerance_N),
    }


def summarize_hysteresis(H_values, H_nominal):
    values = np.asarray(H_values, dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Valori H assenti o non finiti")
    if not np.isfinite(H_nominal) or H_nominal == 0:
        raise ValueError("H nominale nullo o non finito")
    return {
        "mean_Nmm": float(np.mean(values)),
        "std_Nmm": float(np.std(values, ddof=1)) if values.size > 1 else 0.0,
        "worst_Nmm": float(np.max(values)),
        "robustness_index": float(np.max(values) / H_nominal),
    }


def tornado_coefficients(relative_perturbations, H_values):
    """Fit standardized least-squares sensitivities for the tornado plot."""
    if not relative_perturbations:
        raise ValueError("Nessuna perturbazione disponibile")
    names = list(relative_perturbations[0])
    X = np.asarray(
        [[row[name] for name in names] for row in relative_perturbations],
        dtype=float,
    )
    y = np.asarray(H_values, dtype=float)
    if X.ndim != 2 or X.shape[0] != y.size or X.shape[0] <= X.shape[1] + 1:
        raise ValueError("Campioni insufficienti rispetto ai parametri del modello")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(y)):
        raise ValueError("Dati della regressione non finiti")

    X_std = X.std(axis=0, ddof=1)
    y_std = y.std(ddof=1)
    if np.any(X_std == 0) or y_std == 0:
        raise ValueError("Varianza nulla: regressione di sensibilita' non definita")

    design = np.column_stack((np.ones(X.shape[0]), X))
    beta, *_ = np.linalg.lstsq(design, y, rcond=None)
    standardized = beta[1:] * X_std / y_std
    return [
        {"parameter": name, "standardized_coefficient": float(coefficient)}
        for name, coefficient in sorted(
            zip(names, standardized),
            key=lambda item: abs(item[1]),
            reverse=True,
        )
    ]


def _evaluate_case(simulation_result, params):
    required = {"tip_displacement_mm", "tendon_force_N", "contact_force_N"}
    missing = required.difference(simulation_result)
    if missing:
        raise ValueError(f"Risultato simulazione incompleto: {sorted(missing)}")
    hysteresis = hysteresis_area_Nmm(
        simulation_result["tip_displacement_mm"],
        simulation_result["tendon_force_N"],
    )
    contact = contact_metrics(
        simulation_result["contact_force_N"],
        params["F_contact_th"],
        params["contact_target_force"],
        params["contact_stability_tolerance"],
    )
    return {"H_Nmm": hysteresis, **contact}


def _show_progress(completed, total, label, *, stream=None, interactive=None):
    stream = sys.stderr if stream is None else stream
    interactive = stream.isatty() if interactive is None else interactive
    width = 30
    filled = int(width * completed / total)
    bar = f"[{'#' * filled}{'-' * (width - filled)}]"
    message = f"{bar} {completed}/{total} ({100 * completed / total:5.1f}%) {label}"
    if interactive:
        print(f"\r{message}", end="", file=stream, flush=True)
        if completed == total:
            print(file=stream, flush=True)
    else:
        print(message, file=stream, flush=True)


def make_sofa_simulator(contact_force_provider):
    """Build a SOFA case runner using a calibrated contact-force sensor callback.

    The callback receives the SOFA root node after each animation step and must
    return a scalar sensor reading; ``force_scale`` converts it to newtons.
    """
    if not callable(contact_force_provider):
        raise TypeError("contact_force_provider deve essere una funzione")

    def simulate_case(params, steps):
        import Sofa.Core
        import Sofa.Simulation

        from main_scene import createScene

        params = params.copy()
        progress_callback = params.pop("_progress_callback", None)
        params["contact_force_provider"] = contact_force_provider
        root = Sofa.Core.Node("root")
        createScene(root, params)
        Sofa.Simulation.init(root)
        controller = root.getObject("ActuationController")
        if controller is None:
            raise RuntimeError("ActuationController non trovato nella scena SOFA")
        controller.autonomous_mode = True

        dt = float(root.dt.value)
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError(f"Passo temporale SOFA non valido: {dt}")
        try:
            report_every = max(1, steps // 20)
            for step_index in range(steps):
                Sofa.Simulation.animate(root, dt)
                completed_steps = step_index + 1
                if progress_callback is not None and (
                    completed_steps % report_every == 0
                    or completed_steps == steps
                ):
                    progress_callback(completed_steps, steps)
            if len(controller.history_contact_force_N) != steps:
                raise RuntimeError(
                    "Campioni forza di contatto mancanti: verificare "
                    "onAnimateEndEvent e il callback del sensore."
                )
            return {
                "tip_displacement_mm": list(controller.history_tip_y),
                "tendon_force_N": list(controller.history_tendon_force_N),
                "contact_force_N": list(controller.history_contact_force_N),
            }
        finally:
            del controller
            del root
            gc.collect()

    return simulate_case


def run_monte_carlo(
    simulate_case: Callable,
    *,
    seed=None,
    count=MC_NUM_REALIZATIONS,
    output_path=None,
):
    """Run nominal and Monte Carlo cases for soft/rigid and both spring modes.

    ``simulate_case(params, steps)`` must run one complete load/unload cycle in
    SOFA and return tip displacement, tendon force, and measured contact force
    arrays. Contact-force samples must be calibrated to newtons.
    """
    seed = NOMINAL_PARAMS["random_seed"] if seed is None else seed
    realizations = sample_realizations(seed=seed, count=count)
    results = {}
    total_cases = 4 * (count + 1)
    completed_cases = 0
    _show_progress(completed_cases, total_cases, "Avvio simulazioni")

    for contact_object in ("soft", "rigid"):
        for switched_mode in (False, True):
            scenario = f"{contact_object}_{'switched' if switched_mode else 'fixed'}"
            base_params = NOMINAL_PARAMS.copy()
            base_params["contact_object"] = contact_object
            base_params["switched_mode"] = switched_mode
            dt = float(base_params["dt"])
            rate = float(base_params["auto_rate"])
            stroke = float(base_params["u_max"])
            if dt <= 0 or rate <= 0 or stroke <= 0:
                raise ValueError("dt, auto_rate e u_max devono essere positivi")
            steps = 2 * int(np.ceil(stroke / (rate * dt)))
            if steps < 3:
                raise ValueError("Il ciclo completo richiede almeno tre passi")

            def report_case_progress(completed_steps, total_steps, case_label):
                progress = completed_steps / total_steps
                _show_progress(
                    completed_cases + progress,
                    total_cases,
                    f"{scenario}: {case_label} ({completed_steps}/{total_steps} passi)",
                )

            nominal_params = base_params.copy()
            nominal_params["_progress_callback"] = lambda done, total: report_case_progress(
                done, total, "nominale"
            )
            nominal = _evaluate_case(simulate_case(nominal_params, steps), base_params)
            completed_cases += 1
            _show_progress(
                completed_cases,
                total_cases,
                f"{scenario}: nominale",
            )
            cases = []
            perturbation_rows = []
            for realization in realizations:
                params = realization["params"].copy()
                params["contact_object"] = contact_object
                params["switched_mode"] = switched_mode
                params["_progress_callback"] = lambda done, total, run_id=realization["run_id"]: report_case_progress(
                    done,
                    total,
                    f"realizzazione {run_id + 1}/{count}",
                )
                case = _evaluate_case(simulate_case(params, steps), params)
                cases.append({
                    "run_id": realization["run_id"],
                    **case,
                })
                perturbation_rows.append(realization["relative_perturbations"])
                completed_cases += 1
                _show_progress(
                    completed_cases,
                    total_cases,
                    f"{scenario}: realizzazione {realization['run_id'] + 1}/{count}",
                )

            H_values = [case["H_Nmm"] for case in cases]
            results[scenario] = {
                "nominal": nominal,
                "realizations": cases,
                "summary": summarize_hysteresis(H_values, nominal["H_Nmm"]),
                "tornado": tornado_coefficients(perturbation_rows, H_values),
            }

    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8") as output_file:
            json.dump(results, output_file, indent=2, ensure_ascii=False)
    return results
