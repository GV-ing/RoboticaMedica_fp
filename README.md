# Progetto N5 — Robot a tendine con molla in serie

Questo progetto simula in SOFA un tubo flessibile in nitinol azionato da un
tendine. Confronta una molla a rigidezza fissa con una molla a due livelli,
esegue simulazioni con parametri perturbati e raccoglie metriche di isteresi e
contatto.

> **Stato dei risultati:** il codice genera dati utili per confrontare il
> modello, ma non certifica da solo il comportamento fisico del dispositivo.
> In particolare, la conversione della reazione di contatto in newton (`force_scale`),
> la definizione accademica di H e il criterio di stabilità devono essere
> verificati con le specifiche del corso.

## Contenuto del progetto

```text
SOFA_Project/
├── main_scene.py             # Costruzione della scena SOFA
├── run_monte_carlo.py        # Avvio campagna Monte Carlo
├── mc_results.json           # Risultati dell'ultima campagna completata
├── meshes/
│   └── tube_surface.obj      # Superficie del tubo, generata all'avvio
├── modules/
│   ├── config.py             # Parametri nominali e parametri incerti
│   ├── controller.py         # Comando, molla virtuale e registrazione dati
│   ├── contact_force.py      # Lettura delle reazioni di contatto SOFA
│   ├── geometry_gen.py       # Generazione mesh tetraedrica e OBJ
│   └── monte_carlo.py        # Campionamento, simulazioni e metriche
└── tests/
    └── test_project_n5.py    # Test delle funzioni di analisi e delle mesh
```

Il launcher SOFA usato in questa installazione si trova in
`~/SOFA_v26.06.00_Linux_Python3.12/SOFA/run_sofa.sh`.

## Requisiti e ambiente

Il progetto è stato predisposto per SOFA 26.06, SofaPython3 e Python 3.12.
La campagna Monte Carlo deve usare il Python che può importare i moduli
`Sofa.Core` e `Sofa.Simulation`. Nell'installazione in uso è `/usr/bin/python3`;
le librerie SOFA sono rese disponibili tramite le variabili d'ambiente sotto.

Aprire un terminale e preparare l'ambiente:

```bash
cd ~/SOFA_v26.06.00_Linux_Python3.12/SOFA_Project

export SOFA_ROOT="$HOME/SOFA_v26.06.00_Linux_Python3.12/SOFA_v26.06.00_Linux"
export PYTHONPATH="$SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages:$SOFA_ROOT/plugins/STLIB/lib/python3/site-packages:$SOFA_ROOT/plugins/SoftRobots/lib/python3/site-packages:$PWD:$PYTHONPATH"
export LD_LIBRARY_PATH="$SOFA_ROOT/lib:$SOFA_ROOT/plugins/SofaPython3/lib:$LD_LIBRARY_PATH"
```

Le variabili vanno impostate in ogni nuovo terminale prima di eseguire script
SOFA. Se si usa un'altra installazione, sostituire i percorsi corrispondenti.

## Comandi

### Avvio interattivo della scena

Il launcher apre la scena grafica: non è necessario avviarlo prima della
campagna Monte Carlo.

```bash
cd ~/SOFA_v26.06.00_Linux_Python3.12/SOFA
./run_sofa.sh ../SOFA_Project/main_scene.py
```

Comandi tastiera del controller:

| Tasto | Azione |
|---|---|
| Freccia su, `U` o `+` | Aumenta il comando |
| Freccia giù, `D` o `-` | Riduce il comando |
| `R` | Riporta comando e cavo a zero |
| `M` o `A` | Attiva/disattiva il ciclo autonomo |

Il movimento manuale è limitato da `manual_rate` e il comando massimo da
`u_max`. Il controller stampa alcuni valori di stato nella console.

### Campagna Monte Carlo

Dopo aver preparato l'ambiente nella sezione precedente, eseguire:

```bash
cd ~/SOFA_v26.06.00_Linux_Python3.12/SOFA_Project
/usr/bin/python3 -u run_monte_carlo.py --output mc_results.json
```

Il runner crea e simula le scene direttamente: **non** avviare anche la scena
grafica. La barra di avanzamento riporta i passi del caso corrente e il
progresso totale. Con i parametri nominali ci sono 4 scenari, ognuno con una
prova nominale e 20 perturbazioni, per un totale di 84 casi. Ogni ciclo
nominale pianificato contiene 1.600 passi. La campagna può quindi richiedere
tempo.

Opzioni:

```bash
/usr/bin/python3 -u run_monte_carlo.py --help
/usr/bin/python3 -u run_monte_carlo.py --seed 123 --output risultati.json
```

Se si omette `--output`, il percorso predefinito è `mc_results.json`. Il file
JSON viene scritto al termine della campagna; una nuova esecuzione sullo
stesso percorso lo sostituisce.

### Test e mesh

Eseguire i test dall'interno di `SOFA_Project`:

```bash
python3 -m unittest discover -s tests -v
```

Rigenerare manualmente la mesh (la scena la rigenera comunque all'avvio):

```bash
python3 modules/geometry_gen.py
```

Questo comando aggiorna `meshes/tube_surface.obj`.

Controllare che Python 3.12 trovi i binding SOFA:

```bash
/usr/bin/python3 -c 'import Sofa.Core, Sofa.Simulation; print("Import SOFA OK")'
```

Se il comando non funziona, verificare `SOFA_ROOT`, `PYTHONPATH` e
`LD_LIBRARY_PATH`.

## Nodi della scena SOFA

La scena è costruita dalla funzione `createScene` in `main_scene.py`. Le
opzioni fisiche e geometriche vengono prese da `NOMINAL_PARAMS` in
`modules/config.py`; eventuali `custom_params` possono sovrascriverle.

```text
root
├── NitinolTube
│   ├── Tendon
│   ├── VisualModel
│   └── CollisionModel
├── SoftContactObject oppure RigidContactObject (opzionale)
│   ├── CollisionModel
│   └── VisualModel
└── ActuationController
```

### `root`

- Imposta il passo temporale `dt` e la gravità nulla.
- Carica i plugin SOFA richiesti.
- Contiene `FreeMotionAnimationLoop` e `BlockGaussSeidelConstraintSolver`.
- Il solver usa `computeConstraintForces=True` per rendere disponibili le
  reazioni vincolari.
- Definisce la pipeline di collisione e la risposta con attrito.

### `NitinolTube`

Contiene il modello meccanico del tubo:

- mesh tetraedrica (`MeshTopology`);
- gradi di libertà 3D (`MechanicalObject`);
- massa e integrazione implicita;
- `TetrahedronFEMForceField` con modulo iniziale `E_load` e coefficiente di
  Poisson configurato;
- `BoxROI` e vincolo che bloccano la porzione prossimale;
- correzione del solver lineare.

Il controller cambia il modulo di Young tra `E_load` e `E_unload` in base al
verso di movimento. Questo non è un modello costitutivo completo del nitinol:
nel codice l'isteresi del materiale entra anche come termine di forza nella
legge semplificata del tendine.

### `NitinolTube/Tendon`

Contiene i punti del tendine, `CableConstraint` e la mappatura sul tubo.
L'offset laterale iniziale è `yt_offset = 0.60 mm`. Il controller aggiorna il
valore di forza assegnato al cavo.

La molla in serie è **virtuale**, calcolata in `TendonSpringController`: non è
un elemento FEM o un componente molla indipendente nella scena. Il calcolo
usa l'allungamento, lo smorzamento, i termini di isteresi e la rigidezza
selezionata.

### `NitinolTube/VisualModel`

Carica l'OBJ e lo visualizza. La mappatura segue il tubo meccanico.

### `NitinolTube/CollisionModel`

Usa la superficie triangolata del tubo, oltre a modelli di linea e punto,
per individuare collisioni con l'oggetto.

### `SoftContactObject`

Creato se `contact_object = "soft"`:

- griglia regolare di elementi esaedrici;
- FEM con `E_soft = 0.005 MPa` (5 kPa);
- faccia posteriore bloccata;
- superficie di collisione triangolata mappata sulla griglia FEM.

### `RigidContactObject`

Creato se `contact_object = "rigid"`: usa una superficie di collisione
triangolata mantenuta fissa. Non ha un modello FEM deformabile.

### Nessun oggetto

Con `contact_object = "none"` la scena non crea l'oggetto di contatto; la
lettura della forza di contatto restituisce zero. Questo modo è utile per
prove senza contatto, ma non è incluso nei quattro scenari Monte Carlo
predefiniti.

### `ActuationController`

Risponde agli eventi di animazione e ai tasti. In modalità autonoma genera il
ciclo di comando, applica una limitazione di velocità, aggiorna la forza del
tendine e registra posizione, tensione e forza di contatto. Nella versione
attuale tenta anche di mantenere il target di contatto di 20 mN per 5 secondi
dopo il rilevamento del contatto.

## Parametri e modello della molla

Unità dichiarate dal progetto: mm, N, MPa (equivalenti a N/mm²) e s.

Valori nominali principali:

| Parametro | Valore | Significato |
|---|---:|---|
| `E_load` | 77.000 MPa | Modulo usato durante il carico |
| `E_unload` | 65.000 MPa | Modulo usato durante lo scarico |
| `yt_offset` | 0,60 mm | Offset laterale del tendine |
| `k_spring` | 1,2 N/mm | Rigidezza configurazione fissa |
| `k_low` | 0,6 N/mm | Livello scelto in moto libero |
| `k_high` | 4,8 N/mm | Livello scelto quando è rilevato contatto |
| `F_contact_th` | 0,005 N | Soglia di commutazione (5 mN) |
| `Fh_mat` | 0,40 N | Termine di isteresi attribuito al materiale |
| `Fh_spring` | 0,45 N | Termine di isteresi attribuito alla molla |
| `c_spring` | 1,5 N·s/mm | Smorzamento |
| `E_soft` | 0,005 MPa | Modulo dell'oggetto soffice |

In configurazione fissa, il controller applica sempre `k_spring`. In
configurazione switched, seleziona `k_low` sotto soglia e `k_high` sopra
soglia. Le rigidezze perturbate sono assegnate all'inizio di ogni
realizzazione; la soglia rimane il valore nominale di 5 mN.

## Campionamento e analisi Monte Carlo

`modules/monte_carlo.py` genera 20 realizzazioni riproducibili usando il seed
configurato (predefinito 2026). I parametri elencati in `UNCERTAIN_PARAMS`
sono perturbati individualmente: il modulo della variazione è uniforme tra
10% e 20%, con segno casuale.

La lista attuale comprende:

```text
E_load, E_unload, Fh_mat, Fh_spring, c_spring,
k_spring, k_low, k_high, yt_offset
```

Il campionamento ripete l'estrazione di `E_load` e `E_unload` finché
`E_load > E_unload`; quindi questi due campioni non sono completamente
indipendenti. La lista effettiva dei parametri incerti va confrontata con la
tabella ufficiale dell'esercizio.

Per ogni scenario si simula un caso nominale e le 20 realizzazioni. Gli
scenari sono:

| Nome JSON | Oggetto | Molla |
|---|---|---|
| `soft_fixed` | Soffice | Fissa |
| `soft_switched` | Soffice | Due livelli |
| `rigid_fixed` | Rigido | Fissa |
| `rigid_switched` | Rigido | Due livelli |

### Metrica H

Il codice calcola un integrale discreto della forza del tendine rispetto allo
spostamento laterale della punta:

```text
H = abs( somma_i 0.5 * (T_i + T_(i+1)) * (y_(i+1) - y_i) )
```

L'unità risultante è N·mm. Verificare che questa sia la definizione di H
richiesta nel materiale del corso. Inoltre, l'integrale numerico dipende dal
fatto che la traiettoria registrata rappresenti davvero un ciclo completo.

### Stabilità del contatto

La forza è letta da `modules/contact_force.py` usando la componente Y dei
`lambda` sulla superficie di contatto e viene moltiplicata per `force_scale`.
Il valore predefinito è `1.0` e **non costituisce una calibrazione
sperimentale**.

Per la metrica nel JSON:

1. un campione è in contatto se la forza è almeno `F_contact_th` (5 mN);
2. si calcolano media, deviazione standard e scarto medio assoluto dal target
   (`contact_target_force`, nominalmente 20 mN) sui campioni sopra soglia;
3. il caso risulta `stable` se deviazione standard e scarto dal target sono
   entrambi entro `contact_stability_tolerance` (nominalmente 5 mN).

Questa è la regola implementata, non necessariamente la definizione ufficiale
del corso. La routine del controller che legge il provider, nella versione
attuale, converte eventuali eccezioni in forza zero: un guasto di lettura può
quindi apparire come assenza di contatto invece che come errore evidente.

### Statistiche e tornado

Per ciascuno scenario vengono salvati:

- media, deviazione standard e massimo di H sulle 20 realizzazioni;
- indice di robustezza `H_worst / H_nominal`;
- coefficienti di regressione standardizzati, ordinati per valore assoluto.

Il tornado contiene coefficienti numerici nel JSON, non un grafico immagine.
Le due voci con valore assoluto maggiore sono le candidate dominanti. Con
solo 20 campioni e nove parametri, la regressione è esplorativa: non prova
causalità né garantisce una graduatoria statisticamente stabile.

Il JSON non contiene un campo aggregato `unstable_count`; si ricava contando
le realizzazioni con `stable == false`.

## Leggere il JSON e verificare i target

Ogni scenario contiene `nominal`, `realizations`, `summary` e `tornado`.
Dopo aver completato una campagna, questo comando stampa i target di H e il
numero di casi instabili:

```bash
python3 - <<'PY'
import json

with open("mc_results.json", encoding="utf-8") as file:
    results = json.load(file)

for scenario, data in results.items():
    summary = data["summary"]
    cases = data["realizations"]
    unstable = sum(not case["stable"] for case in cases)
    nominal_h = data["nominal"]["H_Nmm"]
    robust_h_ok = summary["worst_Nmm"] <= 1.3 * nominal_h
    contact_ok = unstable == 0
    dominant = data["tornado"][:2]

    print(
        f"{scenario}: H media={summary['mean_Nmm']:.5g}, "
        f"std={summary['std_Nmm']:.5g}, "
        f"worst={summary['worst_Nmm']:.5g}, "
        f"indice={summary['robustness_index']:.4f}, "
        f"instabili={unstable}/20, "
        f"H_target={'OK' if robust_h_ok else 'NO'}, "
        f"contatto={'OK' if contact_ok else 'NO'}"
    )
    print("  due coefficienti tornado principali:", dominant)
PY
```

Il target complessivo è raggiunto solo da una configurazione che soddisfa
**entrambi** i criteri:

```text
H_worst / H_nominal <= 1.3
numero di realizzazioni instabili = 0 su 20
```

## Limiti da tenere presenti

1. **Calibrazione delle forze:** verificare `force_scale` prima di interpretare
   la soglia di 5 mN o il target di 20 mN come forze fisiche in newton.
2. **Cavo:** se `CableConstraint.displacement` non è disponibile, il controller
   avvisa e usa uno spostamento geometrico approssimato.
3. **Molla:** la legge della molla è calcolata dal controller, non modellata
   come componente meccanico separato.
4. **Ciclo e mantenimento:** il controller può aggiungere fino a 5 secondi di
   mantenimento del contatto; la durata prevista dal runner è calcolata dal
   solo movimento di andata e ritorno. Verificare nei dati che il ciclo sia
   completo prima di usare H come area di isteresi.
5. **Eccezioni nel sensore:** un errore del callback di forza viene convertito
   in zero nella versione corrente del controller. Controllare i warning e
   la serie `contact_force_N` per distinguere mancato contatto da lettura
   fallita.
6. **Incertezza:** verificare che la lista dei nove parametri perturbati e il
   vincolo `E_load > E_unload` coincidano con la tabella e le regole del corso.
7. **Oggetto rigido/soffice:** sono modelli semplificati per confronto; la
   geometria, il fissaggio e le proprietà del contatto possono influenzare
   sensibilmente i risultati.

I test automatici controllano funzioni di campionamento e analisi, provider e
validità della mesh. Non certificano da soli la calibrazione delle forze, la
correttezza delle ipotesi fisiche o il raggiungimento dei target di progetto.
