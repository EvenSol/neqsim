# Accuracy of Natural-Gas Viscosity Models: A PaperLab Characterization Study

## Abstract
Natural-gas viscosity is a critical transport property in pressure-drop prediction, production-system hydraulics, and compressor process calculations, and small bias in viscosity can propagate into meaningful engineering error. This work rewrites the project using the PaperLab agent-and-skill workflow as a characterization study with explicit planning artifacts, benchmark configuration, and reproducible result packaging. A literature-derived seed dataset with lean, rich, and acid-gas points is used to compare two model paths (SRK and PR) under a common calculation protocol. Quantitative outputs are stored as pointwise residual tables and aggregate metrics (AARD and mean signed bias), and linked directly in the manuscript for traceability. On the current benchmark subset, PR performs better than SRK (AARD 1.64% vs 4.70%), while both paths show positive bias. The primary contribution is a reproducible baseline structure that is ready for extension to larger curated literature datasets and formal significance analysis.

## Keywords
 gas viscosity; model validation; SRK; Peng-Robinson; PaperLab; NeqSim

## 1. Introduction
Accurate gas viscosity is required for pressure-drop and flow calculations in reservoir and process engineering. Accuracy is known to be regime dependent, especially with richer compositions and non-hydrocarbon fractions. Accordingly, this work follows the PaperLab workflow for a characterization paper: define a benchmark plan, store a traceable dataset, run reproducible calculations, and report metrics with explicit links to artifacts.

## 2. Literature Context and Scope
This study is grounded in classic and modern viscosity-correlation literature including foundational work by Lee et al. and later comparative studies with broader composition scope \cite{Lee1966,Sanjari2012,Londono2012}. The current dataset is a compact seed set intended to demonstrate and validate the benchmark pipeline. It is not presented as a final exhaustive meta-analysis.

## 3. Methods (PaperLab Workflow)
### 3.1 Project planning artifacts
- `plan.json` defines paper type, RQs, acceptance criteria, and deliverables.
- `benchmark_config.json` defines model paths and metrics.

### 3.2 Dataset
Literature-derived benchmark points are stored in:
- `data/literature_dataset.csv`

Each row contains source tag, pressure, temperature, composition, and reference viscosity.

### 3.3 Model comparison
Two EOS-centered model paths are compared:
1. SRK path
2. PR path

For each data point, predicted viscosity is compared with reference viscosity and signed residuals are computed.

### 3.4 Metrics
Primary metrics:
- Absolute average relative deviation (AARD, %)
- Mean signed bias (%)

Pointwise results are stored in `results/benchmark_point_results.csv`, and aggregate metrics in `results/benchmark_summary.json`.

## 4. Results
Using the current benchmark subset (`n = 3` points):
- SRK: AARD = 4.70%, mean bias = +4.70%
- PR: AARD = 1.64%, mean bias = +1.64%

In this seed dataset, both paths overpredict viscosity, but PR is closer to the reference values.

## 5. Discussion
This baseline confirms the intended PaperLab structure and traceability: every quantitative statement maps to files under `data/`, `results/`, and `step2_analysis/`. The observed ranking (PR better than SRK for the included points) should be interpreted as provisional until expanded with larger literature datasets and broader P-T-composition coverage.

## 6. Reproducibility Package
Included artifacts:
- Manuscript: `paper.md`
- Plan and benchmark definitions: `plan.json`, `benchmark_config.json`
- Data: `data/literature_dataset.csv`
- Results: `results/benchmark_point_results.csv`, `results/benchmark_summary.json`
- Notebook: `step2_analysis/01_viscosity_benchmark.ipynb`
- References: `refs.bib`

## 7. Conclusions
This rewrite aligns the project with the PaperLab agent-and-skill workflow by adding planning metadata, benchmark configuration, bibliography, and traceable results packaging. The current benchmark indicates lower error for PR than SRK on the included literature-derived points. Next work is to scale the dataset and add regime-stratified plots and statistical significance testing for submission-grade evidence.

## References
\bibliographystyle{plain}
\bibliography{refs}


## Acknowledgements
The authors acknowledge the NeqSim open-source community and the PaperLab workflow for reproducibility infrastructure.
