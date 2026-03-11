# Technical Report: mAb PBPK-PD Model for Efgartigimod

## 1. Mathematical Framework

### 1.1 Minimal PBPK (Shah & Betts 2012)

The minimal model lumps 12 tissues into two groups based on capillary permeability:

**Tight tissues** (muscle, skin, adipose): High vascular reflection coefficient (sigma_v = 0.95), slow extravasation.

**Leaky tissues** (lung, liver, spleen, gut, kidney, heart): Lower reflection coefficient (sigma_v = 0.85), faster extravasation.

#### Governing equations

**Plasma** (pinocytosis from plasma, Shah & Betts 2012):

```
CL_net = CLp * (1 - FR)    # net clearance after FcRn recycling

dC_p/dt = (1/V_p) * [ -CL_net * C_p
                       + L_tight*(1-sigma_l)*C_tight
                       + L_leaky*(1-sigma_l)*C_leaky
                       - L_tight*(1-sigma_v_tight)*C_p
                       - L_leaky*(1-sigma_v_leaky)*C_p ]
```

**Tight tissues** (muscle, skin, adipose, bone, brain, rest - interstitial only):

```
dC_tight/dt = (1/V_tight) * [ L_tight*(1-sigma_v_tight)*C_p
                               - L_tight*(1-sigma_l)*C_tight ]
```

**Leaky tissues** (lung, liver, spleen, gut, kidney, heart - interstitial only):

```
dC_leaky/dt = (1/V_leaky) * [ L_leaky*(1-sigma_v_leaky)*C_p
                               - L_leaky*(1-sigma_l)*C_leaky ]
```

Where:
- `V_p`, `V_tight`, `V_leaky`: compartment volumes (L)
- `L`: lymph flow (L/h), total 0.12 L/h distributed by ISF volume
- `sigma_v`: vascular reflection coefficient (tight=0.95, leaky=0.85)
- `sigma_l`: lymphatic reflection coefficient (0.20)
- `CLp`: pinocytosis clearance from plasma (L/h)
- `FR`: fraction recycled via FcRn (0.715 for IgG1)
- Tissues only have convection in and lymph out (no pinocytosis from ISF)

### 1.2 FcRn Recycling Mechanism

The key biological process governing mAb PK:

1. **Pinocytosis**: Fluid-phase uptake of IgG from plasma by vascular endothelium into endosomes
2. **Acidification**: Endosomal pH drops to ~6.0
3. **FcRn binding**: IgG-Fc binds FcRn with high affinity at pH 6.0
4. **Sorting**: FcRn-bound IgG recycled to cell surface; unbound IgG -> lysosomes
5. **Release**: At pH 7.4 (extracellular), IgG dissociates from FcRn

#### Concentration-dependent recycling

At high IgG concentrations, FcRn becomes saturated:

```
FR_effective = (k_recycle * Bound) / (k_recycle * Bound + k_deg * Free)

where Bound = FcRn_total * C_IgG / (Kd + C_IgG)  [Langmuir isotherm]
```

This creates **nonlinear PK at very high IgG doses** (e.g., IVIG therapy) but **linear PK at therapeutic mAb concentrations** (because FcRn is not saturated by a single mAb at therapeutic doses).

### 1.3 Efgartigimod Mechanism

Efgartigimod is unique: it uses FcRn biology as its therapeutic mechanism.

```
Efgartigimod + FcRn <-> Efgartigimod:FcRn  (Kd = 20 nM at pH 6.0)
IgG + FcRn <-> IgG:FcRn                     (Kd = 800 nM at pH 6.0)
```

**Competitive inhibition** (Cheng-Prusoff): Efgartigimod outcompetes IgG for FcRn binding (40x higher affinity), reducing IgG recycling. Since endogenous IgG at ~66,667 nM >> Kd_igg (800 nM), the effective IC50 is much higher than Kd_efg alone:

```
FR_IgG_effective = FR_baseline * (1 + IgG/Kd_igg) / (1 + IgG/Kd_igg + E/Kd_efg)

where:
  IgG/Kd_igg ~= 83 at baseline (66,667/800)
  E = efgartigimod plasma concentration (nM)
  Kd_efg = 20 nM (engineered ABDEG)
  Kd_igg = 800 nM (wild-type IgG1)
```

**IgG clearance** increases when FR drops:

```
k_elim_IgG = k_deg_base * (1 - FR_effective) / (1 - FR_baseline)
```

**IgG turnover:**

```
dIgG/dt = k_syn - k_elim_IgG * IgG
```

### 1.4 TMDD (for target-binding mAbs)

For mAbs with soluble or membrane-bound targets:

```
dAb/dt = [PBPK terms] - kon*Ab*R + koff*RC
dR/dt = ksyn - kdeg*R - kon*Ab*R + koff*RC
dRC/dt = kon*Ab*R - koff*RC - kint*RC
```

Creates nonlinear PK: target-mediated clearance at low concentrations, linear CL when target saturated.

## 2. Physiological Parameters

### 2.1 Tissue volumes (70 kg human, ICRP 2002)

| Tissue | Volume (L) | Fv | Fi |
|--------|-----------|-----|-----|
| Plasma | 3.13 | - | - |
| Lung | 0.76 | 0.55 | 0.30 |
| Liver | 1.69 | 0.21 | 0.16 |
| Spleen | 0.15 | 0.17 | 0.15 |
| Gut | 1.20 | 0.05 | 0.14 |
| Kidney | 0.28 | 0.07 | 0.13 |
| Heart | 0.31 | 0.26 | 0.14 |
| Muscle | 25.00 | 0.04 | 0.12 |
| Skin | 3.30 | 0.05 | 0.38 |
| Adipose | 10.00 | 0.05 | 0.14 |
| Bone | 3.80 | 0.04 | 0.10 |
| Brain | 1.40 | 0.03 | 0.09 |

Fv = vascular fraction, Fi = interstitial fraction

### 2.2 Transport parameters

| Parameter | Value | Unit | Source |
|-----------|-------|------|--------|
| Cardiac output | 362.4 | L/h | Shah & Betts 2012 |
| sigma_v (tight) | 0.95 | - | Shah & Betts 2012 |
| sigma_v (leaky) | 0.85 | - | Shah & Betts 2012 |
| sigma_l | 0.20 | - | Shah & Betts 2012 |
| Total lymph flow | 0.12 | L/h | Shah & Betts 2012 |
| CLp (pinocytosis, full) | 0.0366 | L/h | Shah & Betts 2012 |
| CLp (pinocytosis, minimal) | 0.017 | L/h | Calibrated to match t1/2 ~21d |
| FR (IgG1 baseline) | 0.715 | - | Shah & Betts 2012 |

### 2.3 FcRn parameters

| Parameter | Value | Unit | Source |
|-----------|-------|------|--------|
| FcRn_total | 49,800 | nM | Garg & Balthasar 2007 |
| Kd IgG-FcRn (pH 6.0) | 800 | nM | Consensus |
| Kd IgG-FcRn (pH 7.4) | ~1,000,000 | nM | Vaccaro 2005 |
| Kd Efgartigimod-FcRn (pH 6.0) | 20 | nM | Ulrichts 2018 |
| k_recycle | 0.693 | 1/h | Estimated (t1/2 ~1h) |
| k_degradation | 0.347 | 1/h | Estimated (t1/2 ~2h) |

## 3. Implementation Notes

- ODE solver: `deSolve::lsoda` (adaptive step-size, stiff-capable)
- Multiple dosing: event-based (deSolve events framework)
- Time resolution: 1 hour
- Simulation duration: typically 84 days (12 weeks)
- All scripts self-contained with source() dependencies
- Mass balance verified via dedicated test suite
