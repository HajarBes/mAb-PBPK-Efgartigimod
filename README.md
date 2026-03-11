# mAb PBPK-PD — Efgartigimod FcRn Blockade Model

## Overview

Mechanistic PBPK-PD model for monoclonal antibody disposition, applied to efgartigimod (VYVGART) - an engineered Fc fragment that blocks FcRn-mediated IgG recycling. Built on the Shah & Betts 2012 minimal PBPK framework in R. The model predicts efgartigimod PK, endogenous IgG reduction via competitive FcRn blockade (Cheng-Prusoff), and explores Fc engineering scenarios for next-generation antibody therapeutics.

## Project Structure

```
mAb-PBPK-Efgartigimod/
├── src/
│   ├── 01_minimal_pbpk.R           # Shah & Betts 2012 minimal PBPK (3 compartments)
│   ├── 02_full_pbpk.R              # 12-tissue PBPK with per-tissue distribution
│   ├── 03_fcrn_recycling.R         # Detailed FcRn binding kinetics module
│   ├── 04_tmdd_extension.R         # Target-mediated drug disposition
│   ├── 05_efgartigimod_app.R       # Efgartigimod PK-PD: FcRn blockade + IgG reduction
│   ├── 06_sensitivity_analysis.R   # OAT sensitivity + Fc engineering scenarios
│   └── utils/
│       ├── physiology.R            # Human physiological parameter database
│       └── plotting.R              # Plotting functions
├── 04_applications/
│   ├── special_populations.R       # Pediatric, hepatic, renal, elderly, obese
│   └── fc_engineering_scenarios.R  # Kd sweep + half-life extension predictions
├── 02_pksim/
│   ├── simulation_setup.md         # PK-Sim setup guide for cross-platform comparison
│   ├── pksim_vs_scratch.R          # Overlay plots: custom model vs PK-Sim
│   └── platform_comparison.md      # PK-Sim vs SimCyp vs custom code
├── tests/
│   └── test_mass_balance.R         # Non-negativity, dose proportionality, validation
├── report/
│   ├── model_passport.md           # Context of use, assumptions, V&V, limitations
│   └── technical_report.md         # Full equations, parameter tables, references
├── results/                        # Generated plots
└── references/
    └── key_papers.md
```

## How to Run

Each script has a self-contained main block. From the project root:

```bash
# Minimal PBPK for typical IgG1
Rscript src/01_minimal_pbpk.R

# Efgartigimod PK-PD simulation
Rscript src/05_efgartigimod_app.R

# Sensitivity analysis + Fc engineering
Rscript src/06_sensitivity_analysis.R

# Special populations
Rscript 04_applications/special_populations.R

# Run tests (10/10 expected)
Rscript tests/test_mass_balance.R
```

## Requirements

- R >= 4.0
- Packages: `deSolve`, `ggplot2`, `tidyr`, `scales`, `patchwork`

```r
install.packages(c("deSolve", "ggplot2", "tidyr", "scales", "patchwork"))
```

## Key References

- Shah & Betts (2012). *J Pharmacokinet Pharmacodyn* 39:67-86 (minimal PBPK framework)
- Ulrichts et al. (2018). *J Clin Invest* 128(10):4372-4386 (efgartigimod clinical data)
- Garg & Balthasar (2007). *J Pharmacokinet Pharmacodyn* 34(5):687-709 (FcRn PBPK)
- Mager & Jusko (2001). *J Pharmacokinet Pharmacodyn* 28(6):507-532 (TMDD)

## Data

No proprietary or patient-level data is used. All parameters are from published literature. See `report/technical_report.md` for full parameter tables and sources.
