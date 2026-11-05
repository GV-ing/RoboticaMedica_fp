# Progetto 5 — Incertezza parametrica: Monte Carlo e analisi di sensibilità

Robotica Medica 2026/27 (prof.ssa F. Ficuciello). Tubo in nitinol azionato da un tendine con
molla in serie, simulato in **SOFA v26.06 (SofaPython3, Python 3.12) su Ubuntu 22.04**.

Si confrontano due configurazioni della molla: **(a)** fissa con k = 1.2 N/mm e **(b)** a due
livelli, **4.8 N/mm in moto libero e 0.6 N/mm in contatto**, con commutazione istantanea a 5 mN.
Per ognuna si simulano il caso nominale e 20 realizzazioni con i parametri perturbati del
±(10–20)% con segno casuale. Su ogni realizzazione si calcolano l'isteresi H e la stabilità del
contatto; una regressione lineare dà poi il diagramma tornado.

## Risultati (campagna `results/mc_finale`, seed 2026)

> **Attenzione:** la colonna (b) di questa campagna usa l'ordine precedente
> (0.6 N/mm libero, 4.8 N/mm in contatto, `k_free_level = "low"`). Con l'ordine attuale
> vanno rieseguiti i 63 casi della configurazione (b).

| | (a) fissa 1.2 N/mm | (b) due livelli 0.6/4.8 N/mm |
|---|---:|---:|
| H nominale [% di ΔYmax] | 14.38 | 16.00 |
| H medio ± dev. std | 15.19 ± 1.59 | 18.01 ± 2.24 |
| H caso peggiore | 18.09 (r19) | 21.46 (r19) |
| Indice di robustezza (peggiore / nominale) | **1.258** | 1.341 |
| H nominale / H attuazione diretta (26.75 %) | 0.54 | 0.60 |
| Contatti instabili (rigido + soffice, su 40) | **0** | 0 |
| Due parametri dominanti (tornado) | Fh_mat, Fh_spring | k, Fh_spring |
| Target (RI ≤ 1.3 e 20/20 stabili) | **raggiunto** | non raggiunto |

La configurazione più robusta è la **(a) fissa**. Le figure sono in `results/figures/`:
`fig_tornado`, `fig_H_realizzazioni`, `fig_forze_contatto`, `fig_isteresi_nominale` e
`fig_contatto_nominale`.

## Requisiti della traccia e dove sono soddisfatti

| Requisito | Implementazione |
|---|---|
| Superficie OBJ per i modelli visivo e di collisione, mesh volumetrica per la meccanica | `modules/geometry_gen.py`: tetraedri del tubo graduati (0.4 mm nella parte flessibile), OBJ visivo e di collisione separati in `meshes/` |
| FEM con E diverso in carico e scarico | `TetrahedronFEMForceField`; il controller passa da E_load a E_unload secondo il verso del moto (`reinit()`) |
| `BarycentricMapping`, `CableConstraint` | `main_scene.py`: tendine con yt = 0.6 mm, ancorato alla punta |
| Molla in serie T = k·ΔLs ± isteresi | `modules/spring_model.py`: isteresi di Dahl (σ0 = 10/mm, ±(Fh_mat + Fh_spring)), smorzamento c, \|dT/dt\| ≤ 5 N/s |
| Controller da tastiera | `modules/controller.py` (`Sofa.Core.Controller`); stampa u, T, ΔLs e la punta a ogni comando |
| Grafici in tempo reale + CSV | `tools/live_plot_viewer.py` (struttura del Modulo 06 del corso) + `results/gui_<data>.csv` |
| Oggetti rigido e soffice (5 kPa), `FrictionContactConstraint` μ = 0.2 | `main_scene.py`; soffice esaedrico graduato, stessa distanza iniziale (0.03 mm) per i due oggetti |
| H = max ΔY carico−scarico a pari u, in % di ΔYmax, rispetto a k = 50 N/mm | `modules/metrics.py: hysteresis_H`; caso `direct` nella campagna |
| Contatto stabile: picco-picco < 30 % della media per 5 s, senza perdita | `modules/metrics.py: contact_stability` |
| 2 configurazioni × (1 nominale + 20 realizzazioni), commutazione non ritarata | `modules/monte_carlo.py`, `tools/run_monte_carlo.py` |
| Media, dev. std, caso peggiore, indice di robustezza, tornado, numero di contatti instabili | `results/mc_finale/summary.json`, `tools/make_figures.py` |

## Installazione e ambiente

Prerequisiti: il binario SOFA v26.06 per Python 3.12 in
`~/SOFA_v26.06.00_Linux_Python3.12/`, con questo progetto in `.../SOFA_Project`.
Per i grafici serve il Python di sistema con matplotlib:

```bash
sudo apt install python3-matplotlib python3-tk
cd ~/SOFA_v26.06.00_Linux_Python3.12/SOFA_Project
source env.sh                       # variabili di SOFA per gli script python3.12
python3.12 -m unittest discover -s tests -v      # atteso: OK
```

## 1. Scena interattiva (GUI)

```bash
~/SOFA_v26.06.00_Linux_Python3.12/SOFA_Project/run_gui.sh                  # soffice, molla fissa
./run_gui.sh --object rigid --spring switched     # rigido, molla a due livelli
./run_gui.sh --object none                        # moto libero
```

`run_gui.sh` imposta l'ambiente da solo e si può lanciare da qualunque cartella.
Si aprono la finestra di SOFA e quella dei grafici. Premere **Play**, poi:

| Tasto (finestra SOFA) | Azione |
|---|---|
| Freccia su / `U` / `+` | carico: u + 0.35 mm |
| Freccia giù / `D` / `-` | scarico: u − 0.35 mm |
| `R` | reset di comando e cavo |
| `M` / `A` | modo manuale ↔ ciclo autonomo |
| `R` (finestra dei grafici) | azzera i loop e l'energia dissipata |

La finestra dei grafici mostra:
- a sinistra, il loop della molla T–ΔLs (legge di Hooke, banda ±0.85 N, energia dissipata) e il loop della punta ΔY–u con il valore di H;
- a destra, nel tempo: X, Y, ΔY e ΔX della punta, u e ΔLs, T, la forza di contatto (soglia 5 mN e range 10–50 mN) e k.

Note sulla GUI (vanno dichiarate nel video):
- **Mesh più leggere** (`gui_fast`, `gui_overrides` in `config.py`): tubo con `mesh_h_flex` 0.8 mm (circa l'85 % della rigidezza convergente), soffice con griglia uniforme 5, componenti `Parallel…`. `dt` = 0.01 s e il solver dei vincoli (500 iterazioni, tolleranza 1e-7) restano quelli di consegna.
- **Deformazione mostrata ×50** (`visual_amplification`); la forma reale resta visibile in trasparenza.
- La simulazione va più lenta del tempo reale: il terminale stampa ogni 50 passi una riga `[TEMPO]` con il rapporto rispetto al tempo reale.

## 2. Prove nominali (senza GUI)

```bash
source env.sh
python3.12 -u tools/run_cycle.py --free         # cicli liberi: H per fixed, switched, direct
python3.12 -u tools/run_cycle.py --contact      # contatto rigido e soffice, per entrambe le molle
```

Ciclo libero: carico fino a u_max = 15 mm a 1 mm/s, poi scarico.
Prova di contatto: carico fino a 30 mN o a fine corsa, assestamento di 2.5 s, finestra di
stabilità di 5 s.

## 3. Campagna Monte Carlo (127 casi)

I 127 casi sono (1 nominale + 20 realizzazioni) × 2 molle × 3 prove (libera, rigido, soffice),
più 1 caso di attuazione diretta.

```bash
source env.sh
python3.12 -u tools/run_monte_carlo.py --dry-run                    # elenco dei casi
python3.12 -u tools/run_monte_carlo.py --workers 3 --outdir results/mc_finale
python3.12 -u tools/run_monte_carlo.py --workers 3 --outdir results/mc_finale --resume   # dopo un'interruzione
# rifare solo la (b) in una copia, per conservare la campagna con l'ordine precedente:
cp -r results/mc_finale results/mc_k48_libero
python3.12 -u tools/run_monte_carlo.py --workers 3 --outdir results/mc_k48_libero --resume --rerun --configs switched --tag k48_libero
env -u PYTHONPATH python3 tools/make_figures.py --campaign results/mc_finale --out results/figures
```

- Ogni caso viene aggiunto a `cases.jsonl` appena finisce, quindi `--resume` riparte dal primo caso mancante.
- La campagna completa richiede diverse ore. In una VM conviene usare 2–3 worker: con più worker si sono visti blocchi della VM.
- Le realizzazioni perturbano sei gruppi di parametri, ognuno con un solo fattore: E (E_load ed E_unload), Fh_mat, Fh_spring, c, k (k_spring, k_low, k_high) e yt.

## Verifiche del modello

```bash
source env.sh
python3.12 -u tools/mesh_convergence.py       # convergenza della mesh del tubo
python3.12 -u tools/check_cable.py            # T imposta = T trasmessa (rapporto 1.00)
python3.12 -u tools/check_force_readout.py    # lettura della forza di contatto da constraintForces
```

- Rigidezza del tubo: ΔY a 10 N vale 0.304 mm contro 0.320 mm della teoria della trave (−4.8 %).
- Convenzione SOFA: con `FreeMotionAnimationLoop`, valore del cavo e forze di vincolo sono impulsi, quindi si impone T·dt e si legge λ/dt.

## Struttura

```text
main_scene.py          scena SOFA (createScene); GUI se chiamata da runSofa
run_gui.sh, env.sh     launcher della GUI / ambiente per gli script
modules/               config, geometry_gen, spring_model, controller, contact_force,
                       metrics, protocols, monte_carlo, mc_worker, live_plot,
                       amplified_view, plot_style
tools/                 run_cycle, run_monte_carlo, make_figures, live_plot_viewer, check_*
tests/                 test di molla, metriche e mesh (senza SOFA)
results/               mc_finale/ (summary.json, cases.jsonl, realizations.json), figures/
```

## Ipotesi e limiti

- **Unità:** mm, N, s, MPa, massa in tonnellate.
- **Glossario contro traccia:** dove i valori differiscono (E_unload 66 contro 65 GPa, ν 0.30 contro 0.33, k 1.13 contro 1.2 N/mm, parete 0.215 contro 0.22 mm) si sono usati i valori della traccia.
- **Mesh del soffice non del tutto convergente:** dimezzando il passo la forza cambia di circa −15 %. I confronti tra configurazioni restano validi; le forze assolute sul soffice sono approssimate.
- **Smorzamento c quasi ininfluente su H:** il ciclo è lento, 1 mm/s.
- **Stabilità facile da soddisfare:** il modello è quasi statico, quindi la stabilità non discrimina molto tra le configurazioni.
- **Forza di contatto fuori dal range 10–50 mN:**
  - molla (b) sul rigido in 16 casi su 20 (nominale 53.5 mN), per il salto istantaneo a 4.8 N/mm;
  - soffice in 7 casi su 20 per configurazione (circa 10 mN).

  La traccia indica questo range come range di interesse: non è un target del progetto 5.
