# Model Passport: mAb PBPK-PD Model for Efgartigimod

## 1. Context of Use

**Primary question**: Can a mechanistic PBPK model predict the impact of FcRn blockade by efgartigimod on endogenous IgG reduction, and how does FcRn binding affinity engineering affect therapeutic efficacy?

**Secondary questions**:
- What dose of efgartigimod achieves >70% IgG reduction?
- How does FcRn binding affinity (Kd) relate to IgG suppression depth and duration?
- Are dose adjustments needed in special populations (pediatric, hepatic impairment, obesity)?
- What FcRn binding affinity optimizes half-life extension for non-FcRn-blocking mAbs?

**Decision supported**: Lead candidate selection for FcRn-targeting therapeutics, dose rationale for first-in-human studies, special population predictions.

**Regulatory context**: PBPK-based predictions for special populations and drug-drug interactions are accepted by FDA (2018 PBPK Guidance) and EMA (2018 PBPK Guideline). This model demonstrates the approach; it is not a regulatory submission.

## 2. Model Structure

### 2.1 Efgartigimod PK (Minimal PBPK)
- **Type**: 3-compartment PBPK (plasma + tight tissues + leaky tissues)
- **Framework**: Shah & Betts 2012, J Pharmacokinet Pharmacodyn 39:67-86
- **Compartments**:
  - Plasma (3.13 L)
  - Tight tissues: muscle + skin + adipose (lumped interstitial volume)
  - Leaky tissues: lung + liver + spleen + gut + kidney + heart (lumped interstitial volume)
- **Transport**: Convective extravasation (2-pore model) + lymphatic return
- **Elimination**: Pinocytosis from plasma -> endosomal degradation (no FcRn recycling for efgartigimod) + renal clearance
- **Lymph flow**: 0.12 L/h total, distributed by ISF volume fraction

### 2.2 Full Tissue PBPK
- **Type**: 12-tissue PBPK with individual compartments
- **Tissues**: lung, liver, spleen, gut, kidney, heart, muscle, skin, adipose, bone, brain, rest
- **Each tissue**: vascular + interstitial sub-compartments
- **FcRn recycling**: explicit per-tissue endosomal processing

### 2.3 FcRn Recycling Module
- **Mechanism**: Pinocytosis -> endosomal acidification -> FcRn binding (pH 6.0) -> recycling or degradation
- **Key innovation**: Concentration-dependent recycling fraction (Langmuir binding isotherm)
- **Parameters**: FcRn_total, Kd (pH 6.0), k_recycle, k_degradation

### 2.4 TMDD Extension
- **For**: Therapeutic mAbs with soluble/membrane targets
- **Equations**: Mager & Jusko 2001 full TMDD
- **States**: Free drug, free target, drug-target complex

### 2.5 IgG Turnover PD Model
- **Type**: Indirect response model
- **Mechanism**: Efgartigimod blocks FcRn -> reduced IgG recycling fraction -> increased IgG catabolism
- **Competition**: Cheng-Prusoff competitive inhibition (efgartigimod vs endogenous IgG for FcRn)
- **Key equation**: FR_eff = FR_base * (1 + IgG/Kd_igg) / (1 + IgG/Kd_igg + E/Kd_efg)
- **Baseline**: IgG ~10 mg/mL (66,667 nM), synthesis rate balanced with FcRn-dependent clearance

## 3. Parameter Sources

| Parameter | Value | Source |
|-----------|-------|--------|
| Tissue volumes | ICRP 2002 reference values | ICRP Publication 89 |
| Blood flows | Brown et al. 1997 | Toxicol Ind Health 13:407 |
| Vascular/interstitial fractions | Shah & Betts 2012 | Table 1 |
| FcRn total concentration | 49,800 nM | Garg & Balthasar 2007 |
| IgG-FcRn Kd (pH 6.0) | 800 nM | Literature consensus |
| Efgartigimod-FcRn Kd (pH 6.0) | 20 nM | Ulrichts et al. 2018 |
| IgG baseline | 66,667 nM (~10 mg/mL) | Clinical reference range |
| IgG half-life | 21 days | Standard for IgG1 |
| Pinocytosis CL (full model) | 0.0366 L/h | Shah & Betts 2012 |
| Pinocytosis CL (minimal model) | 0.017 L/h | Calibrated for 3-compartment |
| Recycling fraction (baseline) | 0.715 | Shah & Betts 2012 |

## 4. Verification & Validation

### 4.1 Internal Verification
- Non-negativity: all concentrations >= 0
- Functional tests: half-life range, dose proportionality, IgG reduction/recovery
- Dose proportionality: linear PK confirmed for IgG1 (no TMDD)
- Numerical stability: tested across 0.1-50 mg/kg dose range
- Baseline IgG stability: <1% drift over 60 days without drug

### 4.2 External Validation
- **IgG1 half-life**: Model predicts ~21 days (literature: 18-23 days)
- **Efgartigimod IgG reduction**: Model predicts ~75% reduction at nadir with 10 mg/kg weekly x4 (Ulrichts 2018 observed: ~75% reduction)
- **IgG recovery**: Model predicts return to >80% baseline by day 84

### 4.3 Sensitivity Analysis
- One-at-a-time (OAT) sensitivity with 50% parameter variation
- Most sensitive parameters: Kd_efg, FcRn_total, FR_igg_baseline
- Least sensitive: CL_renal, sigma_v

## 5. Assumptions & Limitations

### Assumptions
1. Rapid equilibrium between vascular and plasma compartments
2. Well-mixed compartments (no spatial gradients within tissues)
3. FcRn binding at quasi-steady state within endosome
4. Competitive inhibition model for efgartigimod-IgG competition
5. IgG subclass-independent behavior (model represents total IgG)
6. Allometric scaling for pediatric predictions
7. Linear pinocytosis (no saturation of fluid-phase endocytosis)

### Limitations
1. **No clinical calibration**: Parameters from literature only, not fitted to individual patient data
2. **Simplified endosomal kinetics**: Steady-state binding assumed, no explicit endosomal trafficking dynamics
3. **No IgG subclass differentiation**: IgG1-4 have different FcRn affinities
4. **No anti-drug antibody (ADA) modeling**: ADA can affect efgartigimod PK
5. **Special populations**: Scaling based on physiological principles, not validated against clinical data
6. **No disease state modeling**: gMG, CIDP patients may have altered IgG turnover
7. **Lumped tissue model**: Minimal PBPK lumps 12 tissues into 2 groups

### What this model is NOT
- Not a regulatory submission model
- Not fitted to individual patient data
- Not a replacement for population PK analysis of clinical trial data
- Not validated for dose selection without additional clinical verification

## 6. Software & Reproducibility

- **Language**: R (version >= 4.0)
- **Dependencies**: deSolve, ggplot2, tidyr, scales, patchwork
- **All code**: Open source, available in this repository
- **Run**: Each script is self-contained with `if (sys.nframe() == 0)` main blocks

## 7. References

1. Shah DK, Betts AM. Towards a platform PBPK model to characterize the plasma and tissue disposition of monoclonal antibodies in preclinical species and human. J Pharmacokinet Pharmacodyn. 2012;39(1):67-86.
2. Ng CM, Stefanich EG, Anand BS, et al. Pharmacokinetics/pharmacodynamics of nondepleting anti-CD4 monoclonal antibody (TRX1) in healthy human volunteers. Pharm Res. 2006;23(1):95-103.
3. Jones HM, et al. A Physiologically-Based Pharmacokinetic Model for the Prediction of Monoclonal Antibody Pharmacokinetics From In Vitro Data. CPT Pharmacometrics Syst Pharmacol. 2019;8(10):738-747.
4. Ulrichts P, et al. Neonatal Fc receptor antagonist efgartigimod safely and sustainably reduces IgGs in humans. J Clin Invest. 2018;128(10):4372-4386.
5. Garg A, Balthasar JP. Physiologically-based pharmacokinetic (PBPK) model to predict IgG tissue kinetics in wild-type and FcRn-knockout mice. J Pharmacokinet Pharmacodyn. 2007;34(5):687-709.
6. Mager DE, Jusko WJ. General pharmacokinetic model for drugs exhibiting target-mediated drug disposition. J Pharmacokinet Pharmacodyn. 2001;28(6):507-532.
7. Pyzik M, et al. The Neonatal Fc Receptor (FcRn): A Misnomer? Front Immunol. 2019;10:1540.
8. Datta-Mannan A. Mechanisms Influencing the Pharmacokinetics and Disposition of Monoclonal Antibodies and Peptides. Drug Metab Dispos. 2019;47(12):1100-1110.

## 8. Author

Hajar Besbassi, PhD
Senior Scientist, University of Antwerp (VAXINFECTIO)
hajar.besbassi@gmail.com
