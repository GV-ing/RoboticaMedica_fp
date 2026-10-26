# modules/config.py
"""
Parametri fisici, geometrici e di incertezza per il Progetto 5:
"Incertezza parametrica: Monte Carlo e analisi di sensibilita'"
Unita' coerenti in SOFA: mm, N, MPa (N/mm^2), s.
"""

NOMINAL_PARAMS = {
    # Geometria Tubo Continuo (Baseline)
    "od": 1.93,                 # Diametro esterno [mm]
    "thickness": 0.22,          # Spessore parete [mm] (ID = 1.49 mm)
    "length_total": 100.0,      # Lunghezza totale [mm]
    "length_shaft": 40.0,       # Albero prossimale vincolato [mm]
    "length_flex": 60.0,        # Parte flessibile distale [mm]

    # Materiale Nitinol
    "E_load": 77000.0,          # Modulo di Young in carico [MPa] (77 GPa)
    "E_unload": 65000.0,        # Modulo di Young in scarico [MPa] (65 GPa)
    "poisson": 0.33,            # Coefficiente di Poisson [-]
    "Fh_mat": 0.40,             # Isteresi del materiale [N] (+/- 0.40 N)

    # Tendine
    "yt_offset": 0.60,          # Offset eccentrico tendine lungo Y [mm]

    # Molla in Serie
    "L0_spring": 20.0,          # Lunghezza a riposo molla [mm]
    "k_spring": 1.2,            # Rigidezza fissa (Config. a) [N/mm]
    "k_low": 0.6,               # Rigidezza bassa (Config. b - moto libero) [N/mm]
    "k_high": 4.8,              # Rigidezza alta (Config. b - in contatto) [N/mm]
    "k_direct": 50.0,           # Attuazione diretta di riferimento [N/mm]
    "c_spring": 1.5,            # Smorzamento molla [N*s/mm]
    "Fh_spring": 0.45,          # Isteresi molla [N] (+/- 0.45 N)
    "switched_mode": False,     # False: config (a) k fissa; True: config (b) due livelli

    # Attuazione e Controllo
    "u_step": 0.05,             # Passo comando motore per pressione tasto [mm]
    "manual_rate": 1.0,         # Velocita' massima risposta ai comandi manuali [mm/s]
    "u_max": 8.0,               # Corsa massima in modo autonomo [mm]
    "auto_rate": 1.0,           # Velocita' di u in modo autonomo [mm/s]
    "F_contact_th": 0.005,      # Soglia rilevamento contatto: 5 mN = 0.005 N

    # Ambiente di Contatto
    "contact_object": "soft",   # "soft" | "rigid" | "none"
    "obj_size": 20.0,           # Dimensione caratteristica oggetto [mm]
    "obj_gap": 0.10,            # Distanza iniziale punta-oggetto lungo Y [mm]
    "obj_overlap_x": 5.0,       # Quanto l'oggetto supera la punta lungo X [mm]
    "obj_grid_n": 5,            # Nodi per lato della griglia FEM dell'oggetto soffice
    "E_soft": 0.005,            # Oggetto soffice: 5 kPa = 0.005 MPa
    "poisson_soft": 0.45,       # Poisson oggetto soffice [-]
    "mu_friction": 0.2,         # Coefficiente di attrito [-]
    "alarm_distance": 0.50,     # LocalMinDistance alarmDistance [mm]
    "contact_distance": 0.02,   # LocalMinDistance contactDistance [mm]
    "contact_target_force": 0.020,  # Forza di prova nella stabilita' di contatto: 20 mN (ipotesi)
    "force_scale": 1.0,         # Fattore di calibrazione lambda -> N (da verificare)
    "contact_stability_tolerance": 0.005,  # Tolleranza assoluta sul target forza [N]
    "random_seed": 2026,        # Seed per riproducibilita' delle realizzazioni
    "mc_output": "mc_results.json",

    # Simulazione e output
    "dt": 0.01,                 # Passo di integrazione temporale [s]
    "live_plot": True,          # Grafici in tempo reale (processo matplotlib separato)
    "plot_every": 5,            # Invia un punto ai grafici ogni N passi
    "print_every": 25,          # In modo autonomo stampa lo stato ogni N passi
}

# Range di perturbazione Monte Carlo (Progetto 5): modulo uniforme tra 10% e 20% con segno casuale
MC_UNCERTAINTY_RANGE = (0.10, 0.20)
MC_NUM_REALIZATIONS = 20

UNCERTAIN_PARAMS = (
    "E_load",
    "E_unload",
    "Fh_mat",
    "Fh_spring",
    "c_spring",
    "k_spring",
    "k_low",
    "k_high",
    "yt_offset",
)