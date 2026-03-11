# PBPK Platform Comparison: PK-Sim vs SimCyp vs Custom Code

## Purpose

This document compares three approaches to PBPK modeling for monoclonal antibodies
and Fc-engineered proteins, with specific attention to the efgartigimod use case.
The goal is to help a modeler (or reviewer) understand when each platform is most
appropriate and what trade-offs are involved.

## The Three Approaches

### 1. PK-Sim (Open Systems Pharmacology)

**What it is**: Free, open-source PBPK platform developed by Bayer Technology
Services and the Open Systems Pharmacology (OSP) community. Includes a
dedicated large-molecule module with FcRn recycling.

**Strengths**:
- Open-source (GPLv2) - fully transparent equations and parameter databases
- Built-in human physiology database (ICRP, PK-Sim Standard Individual)
- Dedicated antibody/protein module with two-pore formalism and endosomal FcRn kinetics
- Validated against clinical data for multiple mAbs (trastuzumab, bevacizumab, etc.)
- MoBi integration for custom model extensions (add compartments, reactions, PD)
- OSPSuite-R package for scripted/batch simulations
- FDA and EMA have accepted PK-Sim-based submissions
- Active community: regular updates, qualification library, published tutorials
- Population simulation capability
- Free - no licensing cost

**Limitations**:
- Windows-only (macOS/Linux users need a VM)
- GUI can be slow for large parameter sweeps
- No native PD module for FcRn blockade effects on endogenous IgG
- Limited built-in statistical estimation (no population NLME fitting)
- Smaller user base than SimCyp in regulatory submissions
- Learning curve for MoBi (needed for non-standard mechanisms)

**Best for**:
- Academic and small-biotech PBPK work (no license cost)
- Transparent, reproducible modeling (open source)
- Standard mAb PK prediction (distribution, half-life)
- First-in-human dose prediction with built-in physiology
- Situations requiring full model transparency for regulatory review

### 2. SimCyp (Certara)

**What it is**: Commercial PBPK platform with the largest validation database.
Industry standard for regulatory PBPK submissions, especially for DDI and
special populations.

**Strengths**:
- Largest built-in validation database (hundreds of compounds, DDI pairs)
- Most extensive regulatory track record - used in the majority of PBPK-based
  FDA/EMA submissions
- Dedicated Biologics module (SimCyp Biologics) for mAbs and proteins
- Built-in FcRn recycling, TMDD, and target-mediated clearance
- Population simulator with covariate models (age, weight, genotype, disease)
- Integrated with Phoenix WinNonlin for PK analysis workflows
- Strong DDI prediction capability (not directly relevant for mAbs but valuable
  for ADCs and bispecifics metabolized by CYPs)
- Dedicated support team and training programs
- Virtual bioequivalence (VBE) module

**Limitations**:
- Expensive commercial license (typically $30,000-80,000+/year per seat)
- Proprietary - equations and database not fully transparent
- Windows-only
- Biologics module is newer and less extensively validated than small-molecule module
- Less flexible for novel mechanisms (efgartigimod FcRn blockade PD requires
  custom scripting or workarounds)
- No open-source option for reproducibility audits

**Best for**:
- Regulatory submissions (largest acceptance track record)
- DDI prediction (CYP-mediated, transporter-mediated)
- Virtual bioequivalence studies
- Large pharma with existing Certara infrastructure
- Pediatric extrapolation with validated population models
- Projects where regulatory precedent matters most

### 3. Custom Code (R, Python, MATLAB)

**What it is**: Building the PBPK model from scratch using ODEs in a
general-purpose programming language. This is the approach used in our
efgartigimod project (R + deSolve).

**Strengths**:
- Full control over every equation, parameter, and assumption
- Complete transparency - every line of code is visible and auditable
- Maximum flexibility for novel mechanisms (e.g., efgartigimod's competitive
  FcRn inhibition with Cheng-Prusoff equation, IgG turnover PD)
- Easy integration with other analyses (sensitivity analysis, Bayesian inference,
  optimization, machine learning)
- No license cost (using open-source languages)
- Cross-platform (R/Python run on any OS)
- Version control with git - full change history
- Can be structured as reproducible research (scripts, tests, documentation)
- Easy to scale on HPC clusters for population simulations
- Ideal for method development and publication

**Limitations**:
- No built-in physiology database (must curate parameters from literature)
- No pre-validated compound library
- Higher risk of coding errors (must write and validate own ODE system)
- Less regulatory precedent than platform-based submissions
- Requires strong programming skills
- No GUI - harder for non-programmers to use or review
- Must build own verification/validation framework
- More effort to document for regulatory review

**Best for**:
- Novel mechanisms not supported by existing platforms
- Method development and academic research
- Full-control scenarios where transparency is paramount
- Integration with bespoke statistical or ML workflows
- Situations where the model must be shared as open-source code
- Teaching and learning PBPK principles (forces understanding of every equation)

## Feature Comparison Table

| Feature | PK-Sim | SimCyp | Custom Code |
|---------|--------|--------|-------------|
| **Cost** | Free (open-source) | Commercial ($30K-80K+/yr) | Free |
| **OS support** | Windows only | Windows only | Any |
| **Source code available** | Yes (GitHub) | No (proprietary) | Yes (by definition) |
| **Built-in physiology** | Yes (ICRP database) | Yes (larger database) | No (manual curation) |
| **mAb/protein module** | Yes (two-pore + FcRn) | Yes (Biologics module) | Build from scratch |
| **FcRn recycling** | Built-in (endosomal kinetics) | Built-in | Must implement |
| **TMDD** | Via MoBi extension | Built-in | Must implement |
| **DDI prediction** | Limited | Extensive (industry standard) | Must implement |
| **Population simulation** | Yes | Yes (validated covariates) | Must implement |
| **Regulatory acceptance** | Yes (growing) | Yes (most extensive) | Case-by-case |
| **Customizability** | Moderate (MoBi for extensions) | Limited (scripting API) | Unlimited |
| **Novel mechanisms** | Requires MoBi | Requires workarounds | Native |
| **Reproducibility** | High (open-source) | Lower (proprietary) | Highest (code = model) |
| **Learning curve** | Moderate | Moderate | High (programming) |
| **Batch/scripted runs** | Yes (OSPSuite-R) | Yes (Lua scripting) | Native |
| **HPC scalability** | Limited | Limited | Full |
| **Community/support** | OSP community, GitHub | Certara support team | Self/collaborators |

## Efgartigimod-Specific Considerations

### Why custom code was necessary for this project

Efgartigimod presents a modeling challenge that is not straightforward in either
PK-Sim or SimCyp:

1. **Competitive FcRn inhibition**: Efgartigimod and endogenous IgG compete for
   the same FcRn binding sites. The PD effect (IgG reduction) depends on the
   ratio of efgartigimod to IgG at the endosomal level. Neither PK-Sim nor
   SimCyp natively models this two-species FcRn competition with a feedback loop
   on endogenous IgG turnover.

2. **IgG turnover dynamics**: The clinically relevant endpoint is IgG reduction
   over time, not just efgartigimod PK. This requires coupling efgartigimod PK
   to an endogenous IgG synthesis-degradation model where degradation rate
   depends on FcRn occupancy. This is a PD model that platform tools do not
   include out of the box.

3. **Cheng-Prusoff competition**: Our model uses the Cheng-Prusoff competitive
   inhibition framework to calculate effective IgG recycling fraction as a
   function of both efgartigimod and IgG concentrations. This bidirectional
   feedback (IgG level affects efgartigimod's effective potency, and
   efgartigimod affects IgG level) requires custom ODEs.

4. **Fc engineering scenarios**: Sweeping FcRn Kd from 1-10,000 nM to explore
   how binding affinity engineering affects both half-life and IgG reduction
   is trivial in custom code but requires significant workarounds in GUI-based
   platforms.

### What PK-Sim adds to this project

Despite the need for custom code for the PD component, PK-Sim provides value:

- **Benchmarking**: Comparing efgartigimod PK predictions (not PD) between our
  lumped 3-compartment model and PK-Sim's organ-resolved model validates our
  structural assumptions.
- **Tissue distribution**: PK-Sim predicts organ-level concentrations that our
  lumped model cannot resolve, which is relevant for safety assessments.
- **Credibility**: Showing agreement with an established platform strengthens
  confidence in the custom model.
- **Physiology validation**: PK-Sim's ICRP-derived physiological parameters
  serve as an independent check on our Shah & Betts 2012 parameter values.

### Recommended hybrid workflow

For a project like efgartigimod, the ideal approach combines platforms:

1. **PK-Sim**: Generate efgartigimod PK predictions using the validated
   large-molecule module. Use organ-resolved output for tissue distribution.
2. **Custom R code**: Implement the FcRn competition and IgG turnover PD model.
   Feed PK-Sim's efgartigimod concentration-time profile as a driver function,
   or compare both PK predictions head-to-head.
3. **SimCyp** (if available): Cross-validate with the Biologics module for
   additional regulatory credibility.

## Regulatory Acceptance

### PK-Sim

- FDA: Accepted in multiple INDs and NDAs. The OSP community maintains a
  qualification library with published verification/validation reports.
- EMA: Accepted. OSP tools are used by several EU pharma companies.
- PMDA (Japan): Limited but growing acceptance.
- Key advantage: Full model transparency allows regulators to inspect every
  equation and parameter.

### SimCyp

- FDA: Most widely used PBPK platform in regulatory submissions. FDA's Office
  of Clinical Pharmacology has published guidance citing SimCyp-based analyses.
- EMA: Extensively used. EMA's PBPK guideline (2018) references SimCyp-type
  analyses.
- PMDA: Accepted in multiple submissions.
- Key advantage: Largest regulatory track record and precedent.

### Custom Code

- FDA: Accepted on a case-by-case basis, provided adequate documentation
  (model equations, parameter justification, verification, sensitivity analysis).
- EMA: Same - requires thorough model qualification documentation.
- Key requirement: Must demonstrate model verification, validation against
  clinical data, and sensitivity analysis. A model passport or technical report
  (as included in this project) is essential.
- Key advantage: Full transparency. Key disadvantage: No pre-established
  regulatory precedent for the specific codebase.

## Practical Decision Framework

Use this flowchart to decide which approach fits your project:

**Is the mechanism standard (typical mAb PK, no novel PD)?**
- Yes -> PK-Sim (free) or SimCyp (if license available)
- No -> Continue below

**Does the project require DDI prediction?**
- Yes -> SimCyp (strongest DDI capability)
- No -> Continue below

**Is the mechanism novel (e.g., FcRn blockade PD, bispecific, ADC payload)?**
- Yes -> Custom code (full flexibility), with PK-Sim/SimCyp for PK benchmarking
- No -> PK-Sim or SimCyp based on availability

**Is this for a regulatory submission?**
- Yes, and regulatory precedent matters -> SimCyp
- Yes, and transparency matters -> PK-Sim (open-source)
- No (academic/internal) -> Custom code or PK-Sim

**Budget constraints?**
- No license budget -> PK-Sim or custom code
- License available -> SimCyp for regulatory, PK-Sim for flexibility

## References

1. Open Systems Pharmacology. PK-Sim and MoBi. https://www.open-systems-pharmacology.org/
2. Certara. SimCyp PBPK Simulator. https://www.certara.com/software/simcyp-pbpk/
3. Niederalt C et al. (2018). A generic whole body physiologically based pharmacokinetic
   model for therapeutic proteins in PK-Sim. J Pharmacokinet Pharmacodyn 45:235-257.
4. Jones H et al. (2015). Physiologically based pharmacokinetic modeling in drug
   discovery and development: a pharmaceutical industry perspective. Clin Pharmacol
   Ther 97:247-262.
5. FDA (2018). Physiologically Based Pharmacokinetic Analyses - Format and Content.
   Guidance for Industry.
6. EMA (2018). Guideline on the reporting of physiologically based pharmacokinetic
   (PBPK) modelling and simulation. EMA/CHMP/458101/2016.
7. Shah DK, Betts AM (2012). Towards a platform PBPK model to characterize the plasma
   and tissue disposition of monoclonal antibodies in preclinical species and human.
   J Pharmacokinet Pharmacodyn 39:67-86.
8. Abduljalil K et al. (2022). Prediction accuracy of PBPK models: assessment and
   reporting. Clin Pharmacol Ther 111:1279-1286.
