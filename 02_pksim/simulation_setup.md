# PK-Sim Simulation Setup: Efgartigimod PBPK

## Overview

This guide walks through setting up an efgartigimod (VYVGART) simulation in PK-Sim
to generate predictions that can be compared against our custom R-based PBPK model.
PK-Sim is the PBPK platform from the Open Systems Pharmacology (OSP) suite
(https://www.open-systems-pharmacology.org/) and includes a dedicated large-molecule
(antibody/protein) module with built-in FcRn recycling.

## 1. Installation

PK-Sim runs on Windows only. If you are on macOS or Linux, you will need a
Windows virtual machine (Parallels, VMware, UTM, or VirtualBox).

1. Go to https://www.open-systems-pharmacology.org/
2. Download the latest OSP Suite installer (PK-Sim + MoBi bundle)
3. Run the installer and follow the default options
4. Launch PK-Sim from the Start menu

Minimum requirements: Windows 10/11, 8 GB RAM, .NET Framework 4.7.2+.

## 2. Create a New Simulation

### 2.1 Start a new project

- File > New Project
- Save with a descriptive name (e.g., "Efgartigimod_10mgkg_IV")

### 2.2 Select simulation type

- Building Blocks > Simulation > Create New Simulation
- Select **Large Molecule** as the simulation type
  - This activates PK-Sim's two-pore formalism and FcRn recycling module
  - Small molecule simulations use a completely different distribution model

## 3. Individual (Population) Setup

### 3.1 Create the individual

- Building Blocks > Individual > Create Individual
- Species: **Human**
- Population: **European (ICRP database)**
- Gender: **Male** (to match our 70 kg reference)
- Age: **30 years** (healthy volunteer, matching Ulrichts 2018)
- Weight: **70 kg**
- Height: leave at population default (~176 cm)

### 3.2 Key physiology to verify

After creation, inspect the individual parameters and confirm:

| Parameter | Expected value | Where to check |
|-----------|---------------|----------------|
| Body weight | 70 kg | Individual > Anatomy |
| Plasma volume | ~3.1 L | Individual > Anatomy > Organ volumes |
| Cardiac output | ~5.5-6.0 L/min | Individual > Physiology > Cardiovascular |
| Hematocrit | ~0.45 | Individual > Physiology |

Note: PK-Sim derives organ volumes, blood flows, and lymph flows from the
individual's demographics using its internal physiology database (PK-Sim
Standard Individual). Our custom model uses Shah & Betts 2012 values for a
70 kg reference - these should be very similar but not identical to PK-Sim's
ICRP-derived values.

## 4. Compound Setup

### 4.1 Create compound

- Building Blocks > Compound > Create Compound
- Name: "Efgartigimod"
- Compound type: **Protein** (not small molecule)

### 4.2 Molecular properties

| Parameter | Value | Notes |
|-----------|-------|-------|
| Molecular weight | **54,000 Da** (54 kDa) | Fc fragment, not full IgG |
| Is small molecule | No | |
| Compound type | Protein | |

### 4.3 FcRn binding parameters

This is the most critical section. In PK-Sim's protein module:

| Parameter | Value | Notes |
|-----------|-------|-------|
| FcRn binding on (kass) | Set to achieve Kd = 20 nM at pH 6.0 | ABDEG-engineered high affinity |
| FcRn binding off (kdiss) | Calculated from kass and Kd | |
| Kd at endosomal pH (pH 6.0) | **20 nM** | Key parameter - 40x tighter than wild-type IgG |
| Kd at plasma pH (pH 7.4) | **~1000 nM** | Releases at neutral pH |

How to set Kd in PK-Sim:
- PK-Sim typically parameterizes FcRn interaction via kass (association rate) and
  kdiss (dissociation rate) at endosomal pH
- Use kass = 8.9e5 M-1 s-1, kdiss = 0.0178 s-1 to achieve Kd = 20 nM
- For pH 7.4: use very high kdiss to represent weak binding at neutral pH

### 4.4 No target binding

- Efgartigimod does not bind a membrane-bound target (it is not a therapeutic antibody)
- Do NOT add any target-mediated disposition
- Leave target binding parameters empty/disabled

### 4.5 Solubility and lipophilicity

- For large molecules, PK-Sim uses the two-pore formalism rather than
  partition coefficients
- Solubility: leave at default (not rate-limiting for IV protein)
- Lipophilicity: not applicable for protein distribution in PK-Sim

## 5. Administration Protocol

### 5.1 Create protocol

- Building Blocks > Administration Protocol > Create
- Route: **Intravenous** (IV infusion)
- Dose: **10 mg/kg** = **700 mg** for 70 kg individual
- Infusion duration: **1 hour** (typical clinical infusion time)
- Dosing schedule: **Weekly x 4 doses**
  - Dose 1: Day 0
  - Dose 2: Day 7
  - Dose 3: Day 14
  - Dose 4: Day 21

### 5.2 Schedule setup

In PK-Sim's protocol builder:
- Schema: Multiple dosing
- Start time: 0 h
- Dose: 700 mg
- Infusion time: 60 min
- Number of repetitions: 4
- Time between repetitions: 168 h (= 7 days)

## 6. Simulation Settings

### 6.1 Build the simulation

- Combine the Individual, Compound, and Protocol into a new Simulation
- Name: "Efgartigimod_10mgkg_weekly_x4"

### 6.2 Simulation duration and resolution

| Setting | Value |
|---------|-------|
| Simulation end time | **84 days** (2016 hours) |
| Output resolution | 60 points per hour during dosing, 1 per hour post-dose |

This matches our R model's 84-day simulation window with hourly resolution.

### 6.3 Key outputs to select

Configure the following outputs for export:
- Venous blood plasma concentration (efgartigimod)
- Peripheral venous blood plasma concentration
- Organ concentrations: liver, muscle, skin, lung, kidney, spleen
- Lymph node concentrations (if available)

### 6.4 IgG dynamics (important limitation)

**PK-Sim does not natively model the PD effect of FcRn blockade on endogenous IgG.**

PK-Sim's large-molecule module handles FcRn recycling for the administered
compound, but it does not include an endogenous IgG turnover model that
responds to FcRn occupancy. This is a key difference from our custom model.

To capture IgG reduction dynamics, you would need to either:
1. Couple PK-Sim with MoBi (the OSP model builder) to add an IgG turnover
   compartment with FcRn competition
2. Export efgartigimod PK from PK-Sim and feed it into our custom IgG
   turnover ODE in R

For this comparison, we focus on **efgartigimod PK** (not IgG PD).

## 7. Running and Exporting Results

### 7.1 Run the simulation

- Click "Run Simulation" (or press F5)
- Wait for completion (~seconds for an individual simulation)

### 7.2 Export to CSV

1. Right-click on the simulation results in the chart
2. Select "Export to Excel" or "Export to CSV"
3. Save as `pksim_efgartigimod_results.csv`

Expected CSV format (column names may vary by PK-Sim version):

```
Time [h], Organism|VenousBlood|Plasma|Efgartigimod|Concentration [mg/l], ...
```

### 7.3 Alternative: use OSPSuite-R

PK-Sim results can also be accessed programmatically via the ospsuite R package:

```r
# install.packages("ospsuite", repos = "https://open-systems-pharmacology.r-universe.dev")
library(ospsuite)

sim <- loadSimulation("path/to/Efgartigimod_simulation.pkml")
results <- runSimulation(sim)
outputValues <- getOutputValues(results)
```

## 8. Parameter Mapping: PK-Sim vs Custom Model

This table shows how our custom R model parameters correspond to PK-Sim's
internal parameterization.

| Our model (R) | PK-Sim equivalent | Notes |
|----------------|-------------------|-------|
| V_plasma (3.13 L) | Organism > VenousBlood > Plasma > Volume | PK-Sim derives from individual |
| V_tight (lumped) | Sum of muscle + skin + adipose ISF | PK-Sim uses per-organ volumes |
| V_leaky (lumped) | Sum of liver + lung + spleen + gut + kidney + heart ISF | PK-Sim uses per-organ volumes |
| L_total (0.12 L/h) | Lymph flow (sum of all organ lymph flows) | PK-Sim calculates from two-pore model |
| sigma_v_tight (0.95) | Vascular reflection coefficient (tight tissues) | PK-Sim: set per organ |
| sigma_v_leaky (0.85) | Vascular reflection coefficient (leaky tissues) | PK-Sim: set per organ |
| sigma_l (0.2) | Lymphatic reflection coefficient | Similar concept in two-pore model |
| CLp_efg (0.015 L/h) | Endosomal uptake clearance | PK-Sim: pinocytosis rate |
| FR_efg (0.0) | Fraction recycled | In PK-Sim: determined by FcRn binding kinetics |
| Kd_efg (20 nM) | FcRn Kd at endosomal pH | Direct input in PK-Sim |
| FcRn_total (49,800 nM) | Endosomal FcRn concentration | PK-Sim has its own estimate |
| CL_renal (0.004 L/h) | GFR-based renal clearance | PK-Sim derives from kidney physiology |

### Key structural differences

1. **Compartments**: Our model uses 3 lumped compartments (plasma, tight ISF,
   leaky ISF). PK-Sim uses a full organ-resolved model (15+ organs).

2. **FcRn recycling**: Our model uses competitive inhibition (Cheng-Prusoff)
   with a recycling fraction. PK-Sim uses explicit endosomal binding kinetics
   with the two-pore formalism.

3. **Lymph flow**: Our model uses a single total lymph flow split by ISF
   volume. PK-Sim calculates organ-specific lymph flows.

4. **Distribution**: Our model uses reflection coefficients. PK-Sim uses the
   two-pore model (small pore + large pore radii).

## 9. What to Expect

For efgartigimod 10 mg/kg IV in a 70 kg healthy volunteer:

| Metric | Our model | PK-Sim (expected) | Literature (Ulrichts 2018) |
|--------|-----------|-------------------|---------------------------|
| Cmax (after first dose) | ~240 ug/mL | ~200-250 ug/mL | ~250 ug/mL |
| Terminal t1/2 | ~4.8 days | ~4-6 days | ~4.8 days |
| AUC0-inf (first dose) | Computed | Computed | Not directly reported |
| Steady-state trough | Computed | Computed | Detectable at 7 days |

Differences of 20-50% between our lumped model and PK-Sim's organ-resolved
model are expected and acceptable, since both use different structural
assumptions for the same underlying physiology.

## 10. Troubleshooting

**"Large molecule" option not available**: Make sure you have PK-Sim v11+ and
selected "Protein" as compound type.

**Very different PK profiles**: Check that FcRn Kd is set at endosomal pH
(pH 6.0), not plasma pH. A common mistake is using the pH 7.4 Kd, which would
give very different recycling behavior.

**Efgartigimod half-life too long**: Make sure FR is effectively zero for
efgartigimod. In PK-Sim, this should happen automatically if the Kd at pH 7.4
is set to show weak/no binding (efgartigimod should NOT be recycled by FcRn -
it binds and gets degraded).

**No IgG reduction output**: As noted above, PK-Sim does not natively predict
the PD effect on endogenous IgG. Use MoBi coupling or our custom R model for
this part.
