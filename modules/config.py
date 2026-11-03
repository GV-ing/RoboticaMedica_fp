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

    # Discretizzazione (studio di convergenza: tools/mesh_convergence.py)
    "mesh_h_flex": 0.40,        # Passo assiale FEM sulla parte flessibile [mm]
    "mesh_n_shaft": 4,          # Elementi sull'albero vincolato (non deformabile)
    "mesh_n_theta": 12,         # Suddivisioni circonferenziali (ininfluenti sulla flessione)
    "visual_h": 2.0,            # Passo assiale della superficie visiva [mm]
    "collision_h": 1.0,         # Passo assiale della superficie di collisione [mm]
    "collision_margin": 5.0,    # La collisione parte da (L - obj_size - margin) [mm]

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
    "c_spring": 1.5,            # Smorzamento [N*s/mm]: agisce sul moto del tendine (-c*du1/dt)
    "Fh_spring": 0.45,          # Isteresi molla [N] (+/- 0.45 N)
    "switched_mode": False,     # False: config (a) k fissa; True: config (b) due livelli

    # Attuazione e Controllo
    "u_step": 0.35,             # Passo comando motore per pressione tasto [mm] (traccia: 0.35 mm)
    "manual_rate": 1.0,         # Velocita' massima risposta ai comandi manuali [mm/s]
    "u_max": 15.0,              # Corsa massima del carrello [mm] (glossario: 0-15 mm)
    "auto_rate": 1.0,           # Velocita' di u nel ciclo autonomo [mm/s] (glossario: 0.5-2 mm/s)
    "auto_hold_s": 0.0,         # Mantenimento a fine carico nel ciclo autonomo [s]
    "auto_repeat": False,       # Ciclo autonomo ripetuto (True) o singolo (False)
    "T_max": 10.0,              # Saturazione della tensione [N] (glossario)
    "max_force_rate": 5.0,      # |dT/dt| massimo [N/s] (glossario, rate-limiter di sicurezza)
    "k_switch_instant": True,   # Commutazione 0.6 <-> 4.8 N/mm istantanea (indicazione del docente)
    "dk_dt_max": 2.0,           # |dk/dt| [N/mm/s] (glossario), usato solo se k_switch_instant = False
    "dahl_sigma0": 10.0,        # Rigidezza di Dahl dell'isteresi [1/mm] (glossario: 5-20)
    "hold_duration": 5.0,       # Finestra di valutazione della stabilita' [s] (traccia)
    "auto_stop_force": None,    # Se impostata [N], il carico autonomo si ferma a questa forza
    "contact_stop_force": 0.030,  # Prova di contatto: arresto a 30 mN (centro del range) o a fine corsa
    "contact_settle_s": 2.5,    # Assestamento prima della finestra (k sale 0.6->4.8 in ~2.1 s)
    "contact_rate": 2.0,        # Velocita' di carico nella prova di contatto [mm/s] (glossario: max 2)
    "stability_ptp_max": 0.30,  # Oscillazione picco-picco massima / media (traccia: 30%)
    "contact_loss_force": 1e-4, # Sotto questa forza [N] il contatto si considera perso
    "F_contact_th": 0.005,      # Soglia rilevamento contatto: 5 mN = 0.005 N

    # Ambiente di Contatto
    "contact_object": "soft",   # "soft" | "rigid" | "none"
    "obj_size": 20.0,           # Dimensione caratteristica oggetto [mm]
    "obj_gap": 0.03,            # Gap punta-oggetto lungo Y [mm], UGUALE per rigido e soffice:
                                #   contatto a T ~ 0.3 N (c = 0.0304 mm/N), stessa geometria per il
                                #   confronto tra oggetti; H si misura comunque in moto libero
    "obj_overlap_x": 5.0,       # Quanto l'oggetto supera la punta lungo X [mm]
    "obj_grid_n": 5,            # Nodi per lato della griglia FEM dell'oggetto soffice
    "soft_mesh": "graded",      # "graded" (fine al contatto) | "uniform" (obj_grid_n per lato)
    "soft_h_min": 1.6,          # Passo minimo al contatto [mm]: compromesso costo/accuratezza (soft_mesh_convergence)
    "soft_h_max": 5.0,          # Passo massimo lontano dal contatto [mm]
    "soft_growth": 1.4,         # Fattore di crescita geometrica del passo
    "E_soft": 0.005,            # Oggetto soffice: 5 kPa = 0.005 MPa
    "poisson_soft": 0.45,       # Poisson oggetto soffice [-]
    "mu_friction": 0.2,         # Coefficiente di attrito [-]
    "alarm_distance": 0.30,     # LocalMinDistance alarmDistance [mm] (> gap, < dimensione elementi)
    "contact_distance": 0.02,   # LocalMinDistance contactDistance [mm]
    "contact_target_force": 0.020,  # Forza di prova nella stabilita' di contatto: 20 mN (ipotesi)
    "force_scale": 1.0,         # Fattore di calibrazione lambda -> N (da verificare)
    "contact_stability_tolerance": 0.005,  # Tolleranza assoluta sul target forza [N]
    "random_seed": 2026,        # Seed per riproducibilita' delle realizzazioni
    "mc_workers": 6,            # Processi paralleli del Monte Carlo (12 thread, 15 GB RAM)

    # Masse da densita' (unita' mm-N-s: massa in tonnellate, densita' in t/mm^3)
    "rho_tube": 6.45e-9,        # Nitinol 6450 kg/m^3
    "rho_soft": 1.0e-9,         # Oggetto soffice ~ tessuto/acqua 1000 kg/m^3

    # Modelli di collisione (Triangle | Line | Point). Verificato con tools/check_scene.py:
    # i modelli completi sono piu' veloci (432 vs 691 ms/passo) e stabili; i ridotti
    # (solo Point vs Triangle) oscillano a forza elevata (61 mN p-p a T = 10 N).
    "tube_collision_models": ("Triangle", "Line", "Point"),
    "object_collision_models": ("Triangle", "Line", "Point"),

    # Simulazione e output
    "dt": 0.01,                 # Passo di integrazione temporale [s]
    "live_plot": True,          # Grafici in tempo reale nella GUI (tools/live_plot_viewer.py)
    "record_csv": True,         # Registrazione CSV nella GUI (results/gui_<data_ora>.csv)
    "plot_every": 5,            # Invia un punto ai grafici ogni N passi
    "print_every": 50,          # In modo autonomo stampa lo stato ogni N passi
}

# Range di perturbazione Monte Carlo (Progetto 5): modulo uniforme tra 10% e 20% con segno casuale
MC_UNCERTAINTY_RANGE = (0.10, 0.20)
MC_NUM_REALIZATIONS = 20

# Parametri incerti della tabella del Progetto 5, raggruppati: ogni gruppo riceve UNA
# perturbazione (modulo uniforme 10-20%, segno casuale) applicata a tutte le sue chiavi.
#   E : E_load ed E_unload con lo stesso fattore (rapporto carico/scarico invariato)
#   k : k fissa e i due livelli 0.6 / 4.8 con lo stesso fattore ("k" e' un solo parametro)
UNCERTAIN_GROUPS = {
    "E": ("E_load", "E_unload"),
    "Fh_mat": ("Fh_mat",),
    "Fh_spring": ("Fh_spring",),
    "c": ("c_spring",),
    "k": ("k_spring", "k_low", "k_high"),
    "yt": ("yt_offset",),
}
