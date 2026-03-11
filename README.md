# Mechanistic PBPK Model for Monoclonal Antibody Disposition and FcRn-Mediated IgG Reduction

## Drug Development Question

**Can a physiologically-based pharmacokinetic model predict the impact of FcRn blockade on endogenous IgG levels, and how does FcRn binding affinity engineering guide the design of next-generation antibody therapeutics?**

This project builds a mechanistic PBPK framework for monoclonal antibodies, applied to efgartigimod (VYVGART) - an engineered Fc fragment that blocks FcRn-mediated IgG recycling. The model predicts:

- Efgartigimod plasma PK across dose levels
- Endogenous IgG reduction following FcRn blockade
- Optimal FcRn binding affinity for therapeutic efficacy
- PK/PD in special populations (pediatric, hepatic impairment, obesity, elderly)

## Why This Matters

The neonatal Fc receptor (FcRn) is the central mechanism governing IgG homeostasis. Understanding FcRn biology quantitatively is essential for:

1. **IgG-reducing therapies** (efgartigimod, rozanolixizumab) - predicting dose-response for autoimmune diseases
2. **Half-life extension** - engineering Fc domains for longer-acting mAbs
3. **Biosimilar development** - ensuring equivalent FcRn binding kinetics
4. **Special population predictions** - PBPK-based dose adjustments when clinical data is limited

## Project Structure

```
mAb-PBPK-Efgartigimod/
|-- src/
|   |-- 01_minimal_pbpk.R           Shah & Betts 2012 minimal PBPK (3 compartments)
|   |-- 02_full_pbpk.R              12-tissue PBPK with per-tissue FcRn recycling
|   |-- 03_fcrn_recycling.R         Detailed FcRn binding kinetics module
|   |-- 04_tmdd_extension.R         Target-mediated drug disposition
|   |-- 05_efgartigimod_app.R       Efgartigimod PK-PD: FcRn blockade + IgG reduction
|   |-- 06_sensitivity_analysis.R   Global SA + Fc engineering scenario analysis
|   |-- utils/
|       |-- physiology.R            Human physiological parameter database
|       |-- plotting.R              Publication-quality plotting functions
|-- 04_applications/
|   |-- special_populations.R       Pediatric, hepatic, renal, elderly, obese
|   |-- fc_engineering_scenarios.R  Kd sweep + half-life extension predictions
|-- tests/
|   |-- test_mass_balance.R         Mass balance, non-negativity, dose proportionality
|-- report/
|   |-- model_passport.md           Context of use, assumptions, verification, limitations
|   |-- technical_report.md         Full equations, parameter tables, references
|-- data/
|-- results/
|-- references/
```

## Model Complexity Ladder

| Level | Script | Compartments | FcRn | Use case |
|-------|--------|-------------|------|----------|
| Minimal PBPK | 01 | 3 (plasma + tight + leaky) | Constant FR | Rapid screening, dose ranging |
| Full PBPK | 02 | 12 tissues | Per-tissue FR | Tissue distribution, biodistribution |
| FcRn module | 03 | Endosomal sub-model | Explicit binding kinetics | FcRn saturation, competition |
| TMDD | 04 | 3 + target | Constant FR + TMDD | Target-binding mAbs |
| Efgartigimod | 05 | 3 + IgG turnover | Competitive inhibition | FcRn-blocking therapeutics |

## Key Results

- **IgG1 half-life**: Model predicts ~21 days (literature: 18-23 days)
- **Efgartigimod 10 mg/kg weekly x4**: ~75% IgG reduction at nadir (Ulrichts 2018: ~75%)
- **Fc engineering**: Lower FcRn Kd predicts deeper IgG suppression, with diminishing returns below ~10 nM
- **Half-life extension**: Enhanced FcRn binding (Kd 800 -> 50 nM) predicts >2x half-life increase
- **Special populations**: Hepatic impairment shows largest deviation from reference adult

## Requirements

```r
install.packages(c("deSolve", "ggplot2", "tidyr", "scales", "patchwork"))
```

## How to Run

Each script has a self-contained main block. From the project root:

```r
# Run minimal PBPK for typical IgG1
source("src/01_minimal_pbpk.R")

# Run efgartigimod PK-PD simulation
source("src/05_efgartigimod_app.R")

# Run sensitivity analysis
source("src/06_sensitivity_analysis.R")

# Run special populations
source("04_applications/special_populations.R")

# Run tests
source("tests/test_mass_balance.R")
```

## Scientific Basis

1. Shah DK, Betts AM (2012). Towards a platform PBPK model to characterize the plasma and tissue disposition of monoclonal antibodies. *J Pharmacokinet Pharmacodyn* 39:67-86.
2. Ulrichts P et al. (2018). Neonatal Fc receptor antagonist efgartigimod safely and sustainably reduces IgGs in humans. *J Clin Invest* 128(10):4372-4386.
3. Garg A, Balthasar JP (2007). Physiologically-based pharmacokinetic (PBPK) model to predict IgG tissue kinetics in wild-type and FcRn-knockout mice. *J Pharmacokinet Pharmacodyn* 34(5):687-709.
4. Mager DE, Jusko WJ (2001). General pharmacokinetic model for drugs exhibiting target-mediated drug disposition. *J Pharmacokinet Pharmacodyn* 28(6):507-532.
5. Pyzik M et al. (2019). The Neonatal Fc Receptor (FcRn): A Misnomer? *Front Immunol* 10:1540.

## Author

Hajar Besbassi, PhD - Senior Scientist, University of Antwerp (VAXINFECTIO)

## License

This project is for educational and portfolio purposes. Physiological parameters are from published literature. No proprietary data or software is used.
