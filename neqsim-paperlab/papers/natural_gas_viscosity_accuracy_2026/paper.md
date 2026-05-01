# Accuracy of Viscosity Models for Natural Gas: Literature Synthesis and NeqSim Benchmarking Study Design

## Abstract
Natural-gas viscosity is a key transport property for pressure-drop prediction, well deliverability, compressor sizing, and process simulation. This paper consolidates published evidence on gas-viscosity model accuracy and translates it into a practical NeqSim benchmarking workflow suitable for journal publication. Literature indicates that classical correlations (e.g., Lee–Gonzalez–Eakin family) remain practical for lean-gas screening, while broader corresponding-states / dense-gas methods are usually needed for rich, acid-gas, and near-critical conditions. Reported errors vary strongly by domain; for example, modern empirical regressions trained on large databases report overall AARD near 2–3% in their fitted windows, while broader non-hydrocarbon extensions are often reported below approximately 6% AARD. Building on this literature, we define a reproducible NeqSim protocol to compare SRK-based, PR-based, and tuned viscosity paths against literature datasets using segmented diagnostics (AARD, RMSE, bias, P10/P50/P90 error, and dense-gas residual trends).

**Keywords:** natural gas viscosity, NeqSim, Lee-Gonzalez-Eakin, corresponding states, dense gas, benchmark

## 1. Introduction
Gas viscosity enters directly in Reynolds number, friction factor, and pressure-gradient calculations. In reservoir and process engineering, even modest viscosity bias can propagate into material-balance interpretation, tubing performance, and compressor duty calculations.

The central issue is not whether one model is “good,” but **where** it is good. Published comparative studies repeatedly show model ranking changes with pressure, temperature, and composition (especially heavy-end and acid-gas content). Therefore, model assessment must be regime-specific rather than based on a single global error value.

## 2. Literature Review of Similar Work
### 2.1 Foundational and widely used models
Classic natural-gas viscosity practice in petroleum engineering starts from Lee–Gonzalez–Eakin (LGE) style formulations and related reduced-property correlations. These methods are simple and robust for conventional lean-gas applications but are less reliable when extrapolated far outside their calibration envelope.

### 2.2 Later empirical and modified correlations
Later studies introduced expanded datasets and alternative functional forms for wider pressure/composition conditions. A frequently cited example reports a new empirical correlation with **AARD = 2.173%** on its validation dataset of 4089 points (Sanjari & Lay, 2012, Journal of Natural Gas Science and Engineering).

### 2.3 Processing-focused model evaluation
Londono et al. (2012, Fluid Phase Equilibria) compared viscosity methods for natural-gas processing and reported that both versions of their correlation achieved **overall AARD below 6%** for non-hydrocarbon components in the assessed scope.

### 2.4 HPHT and sour-gas context in recent reviews
Recent review literature on HPHT and sour-gas applications emphasizes that no universal model dominates all domains; model suitability depends on dense-gas behavior, compositional detail, and calibration coverage.

## 3. Evidence-Based Accuracy Expectations
Based on published comparative literature:
- **Lean methane-rich gas (moderate pressure):** classical correlations are often acceptable for screening-level studies.
- **Rich gas / heavy-end influence:** correlation bias increases unless characterization quality is high.
- **CO2/H2S-containing gas:** acid-gas-aware extensions or denser-physics transport paths are preferred.
- **Near-critical/high-density conditions:** highest error sensitivity; uncertainty reporting is mandatory.

## 4. NeqSim Benchmarking Methodology
### 4.1 Objective
Evaluate which NeqSim model path best reproduces literature reference data per regime, not just globally.

### 4.2 Model paths to compare
1. **Path A:** SRK EOS + default gas transport workflow.
2. **Path B:** PR EOS + default gas transport workflow.
3. **Path C:** tuned transport path (CSP/LBC/Pedersen style parameters where available from PVT data).

### 4.3 Dataset requirements
Each literature point should include:
- T [K], P [bar]
- full composition (mole fractions; include CO2/H2S/N2 where present)
- measured/reference viscosity
- source metadata and uncertainty

### 4.4 Reproducible run sequence in NeqSim
For each model path and data point:
1. Create fluid with identical composition.
2. Apply consistent mixing rule.
3. Run TP flash.
4. Initialize properties.
5. Record gas viscosity in a common unit (cP).

### 4.5 Metrics and statistical tests
Report by regime and globally:
- AARD (%), RMSE (cP), mean signed bias (%), max abs error (%), P10/P50/P90 abs error (%).
- Paired model-path comparison with signed-rank or paired t-test (as appropriate).

## 5. Journal-Ready Results Structure
### 5.1 Mandatory tables
- Regime-wise metric table (lean/rich/acid/near-critical)
- Global summary table
- Worst-case residual table (top 10 states)

### 5.2 Mandatory figures
- Parity plot (predicted vs literature viscosity)
- Abs-error vs pressure and vs reduced density
- Boxplot of abs error by regime

### 5.3 Required discussion points
- Where each model path wins/fails.
- Physical reason for bias trends.
- Engineering implications (pressure drop, flow, compressor work).
- Recommendation matrix for model selection.

## 6. Practical Recommendation Matrix
- **Screening and lean gas:** start with simpler correlation-oriented path; validate with sampled points.
- **Design and broad P–T–x envelope:** prefer EOS/transport workflow with regime validation.
- **Sour and near-critical service:** require segmented validation and explicit uncertainty bands.

## 7. Conclusions
Literature and engineering practice agree on one core point: there is no universal best natural-gas viscosity model across all pressure, temperature, and composition domains. Reported best-in-class results near 2–3% AARD are generally tied to specific fitted windows, while broader applicability studies often report larger but acceptable errors (commonly below about 6% in targeted scopes). A publishable NeqSim study should therefore prioritize regime-wise validation, transparent uncertainty, and reproducible data/model artifacts rather than a single headline metric.

## References
1. Lee, A.L., Gonzalez, M.H., Eakin, B.E. *The Viscosity of Natural Gases*. Journal of Petroleum Technology (classic foundational correlation work).
2. Sanjari, E., Lay, E.N. (2012). *An accurate empirical correlation for predicting natural gas viscosity*. Journal of Natural Gas Science and Engineering, 9, 283–289. (Reported AARD 2.173% on 4089-point dataset).
3. Londono, R.A., Archer, R.A., Blasingame, T.A. (2012). *Viscosity prediction for natural gas processing applications*. Fluid Phase Equilibria, 335, 138–149. (Processing-focused comparisons; reported non-hydrocarbon fits below 6% AARD).
4. Recent HPHT/sour-gas review articles on natural-gas viscosity correlation development and limitations.
5. NeqSim documentation on gas-property and viscosity workflows (SRK/PR transport usage and tuning options).
