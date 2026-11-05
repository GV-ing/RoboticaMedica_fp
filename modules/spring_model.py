# modules/spring_model.py
"""Modello della molla in serie (SEA) tra carrello motore e tendine. Puro Python (testabile).

Riferimenti: traccia del Progetto 5 e glossario del corso (GLOSSARIO_VARIABILI_E_PARAMETRI).

    dL_s   = max(0, u - u1)                                    allungamento molla [mm]
    z      : stato di Dahl, dz/du = sigma0 * (1 - z*sgn(du))   (z in [-1, +1])
    F_hyst = z * (F_h-molla + F_h-mat)                         isteresi totale [N]
    T*     = k * dL_s + F_hyst - c * du1/dt                    (c agisce sul moto del tendine)
    T      = T* limitata in velocita' (|dT/dt| <= max_force_rate) e saturata in [0, T_max]

Configurazione (b) a due livelli: k = k_free in moto libero e k = k_contact quando la forza di
contatto supera la soglia F_contact_th; k_free_level ("high"/"low") sceglie quale dei due livelli
(k_high = 4.8, k_low = 0.6 N/mm) si usa in moto libero. La commutazione e' ISTANTANEA (indicazione del docente);
resta il limite |dT/dt| <= max_force_rate sulla tensione (glossario, sicurezza). La variante
con |dk/dt| <= dk_dt_max e' disponibile con k_switch_instant = False.

Scelte di progetto (concordate):
  * smorzamento c fuori dalla legge della molla, sul moto del tendine (u1);
  * sigma0 = 10 mm^-1 (centro del range 5-20 del glossario);
  * 0.6 N/mm in moto libero, 4.8 N/mm in contatto.
"""

import math


class SeriesSpring:
    def __init__(self, params):
        self.switched = bool(params.get("switched_mode", False))
        self.k_fixed = float(params["k_spring"])
        self.k_low = float(params["k_low"])
        self.k_high = float(params["k_high"])
        level = params.get("k_free_level", "low")
        if level not in ("low", "high"):
            raise ValueError("k_free_level deve essere 'low' oppure 'high'")
        # Livelli scelti dagli studenti: k in moto libero e k dopo il rilevamento del contatto
        self.k_free = self.k_high if level == "high" else self.k_low
        self.k_contact = self.k_low if level == "high" else self.k_high
        self.F_th = float(params["F_contact_th"])
        self.k_switch_instant = bool(params.get("k_switch_instant", True))
        self.dk_dt_max = float(params["dk_dt_max"])
        self.F_hyst_max = float(params["Fh_spring"]) + float(params["Fh_mat"])
        self.sigma0 = float(params["dahl_sigma0"])
        self.c = float(params["c_spring"])
        self.T_max = float(params["T_max"])
        self.max_force_rate = float(params["max_force_rate"])
        for name in ("k_fixed", "k_low", "k_high", "sigma0", "T_max", "max_force_rate"):
            if not getattr(self, name) > 0:
                raise ValueError(f"{name} deve essere positivo")
        if self.F_hyst_max < 0 or self.c < 0 or self.dk_dt_max <= 0:
            raise ValueError("Isteresi e smorzamento non negativi, dk_dt_max positivo")
        self.reset()

    # ------------------------------------------------------------------ stato
    def reset(self):
        self.z = 0.0
        self.k = self.k_free if self.switched else self.k_fixed
        self.T = 0.0
        self.T_star = 0.0
        self.dL_s = 0.0
        self.u_prev = 0.0
        self.u1_prev = 0.0
        self.direction = 1.0

    def k_target(self, contact_force_N):
        if not self.switched:
            return self.k_fixed
        return self.k_contact if contact_force_N >= self.F_th else self.k_free

    # ------------------------------------------------------------------ Dahl
    def update_dahl(self, du):
        """Aggiornamento esatto di Dahl per uno spostamento du a segno costante."""
        if du > 0.0:
            self.z = 1.0 - (1.0 - self.z) * math.exp(-self.sigma0 * du)
            self.direction = 1.0
        elif du < 0.0:
            self.z = -1.0 + (1.0 + self.z) * math.exp(self.sigma0 * du)
            self.direction = -1.0
        self.z = min(1.0, max(-1.0, self.z))

    # ------------------------------------------------------------------ passo
    def max_motor_rate(self, rate):
        """Velocita' del carrello compatibile con |dT/dt| <= max_force_rate (k * du/dt)."""
        return min(rate, self.max_force_rate / max(self.k, 1e-9))

    def step(self, u, u1, contact_force_N, dt):
        """Avanza di dt con comando motore u e spostamento del tendine u1. Restituisce T [N]."""
        if dt <= 0:
            raise ValueError("dt deve essere positivo")
        self.update_dahl(u - self.u_prev)

        k_new = self.k_target(contact_force_N)
        if self.k_switch_instant:
            self.k = k_new                      # commutazione istantanea (docente)
        else:
            lim = self.dk_dt_max * dt           # variante con slew-rate del glossario
            self.k += max(-lim, min(lim, k_new - self.k))

        self.dL_s = max(0.0, u - u1)
        u1_dot = (u1 - self.u1_prev) / dt
        if self.dL_s <= 0.0:
            self.T_star = 0.0  # molla scarica: il cavo non trasmette forza
        else:
            self.T_star = self.k * self.dL_s + self.z * self.F_hyst_max - self.c * u1_dot

        dT = self.T_star - self.T
        lim_T = self.max_force_rate * dt
        self.T = min(self.T_max, max(0.0, self.T + max(-lim_T, min(lim_T, dT))))

        self.u_prev = u
        self.u1_prev = u1
        return self.T
