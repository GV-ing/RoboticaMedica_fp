"""Run the N5 Monte Carlo campaign inside the SOFA Python environment."""

import argparse
import importlib
from pathlib import Path

from modules.config import NOMINAL_PARAMS
from modules.contact_force import read_contact_force
from modules.monte_carlo import make_sofa_simulator, run_monte_carlo


def load_provider(specification):
    module_name, separator, function_name = specification.partition(":")
    if not separator or not module_name or not function_name:
        raise ValueError("Il provider va specificato come package.module:function")
    provider = getattr(importlib.import_module(module_name), function_name)
    if not callable(provider):
        raise TypeError(f"{specification} non identifica una funzione")
    return provider


def main():
    parser = argparse.ArgumentParser(
        description="Esegue i casi Monte Carlo N5 in ambiente SOFA."
    )
    parser.add_argument(
        "--contact-force-provider",
        default="modules.contact_force:read_contact_force",
        help=(
            "Funzione package.module:function; riceve il root SOFA e restituisce "
            "la forza normale in unita' solver (convertita con force_scale). "
            "Default: modules.contact_force:read_contact_force."
        ),
    )
    parser.add_argument(
        "--output",
        default=NOMINAL_PARAMS["mc_output"],
        help="Percorso JSON dei risultati (default: mc_results.json).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=NOMINAL_PARAMS["random_seed"],
    )
    args = parser.parse_args()

    provider = load_provider(args.contact_force_provider)
    simulator = make_sofa_simulator(provider)
    results = run_monte_carlo(
        simulator,
        seed=args.seed,
        output_path=Path(args.output),
    )
    for scenario, data in results.items():
        summary = data["summary"]
        print(
            f"{scenario}: H media={summary['mean_Nmm']:.6g} Nmm, "
            f"dev.std={summary['std_Nmm']:.6g} Nmm, "
            f"worst={summary['worst_Nmm']:.6g} Nmm, "
            f"robustezza={summary['robustness_index']:.6g}"
        )


if __name__ == "__main__":
    main()
