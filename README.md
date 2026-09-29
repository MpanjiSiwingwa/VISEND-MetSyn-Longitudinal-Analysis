
**Reproducible R Analysis Pipeline for Longitudinal Assessment of Metabolic Syndrome Among People Living with HIV in the VISEND Trial, Zambia**

![R Version](https://img.shields.io/badge/R-4.5.1-blue?logo=r)
![Clinical Research](https://img.shields.io/badge/Clinical-HIV%20Research-red)
![Project Status](https://img.shields.io/badge/Status-Research-orange?logo=github)
![License: MIT](https://img.shields.io/badge/License-MIT-green?logo=open-source-initiative)

---

## 📑 Table of Contents

* [Project Overview](#-project-overview)
* [Study and Data Source](#-study-and-data-source)
* [Analysis Objectives](#-analysis-objectives)
* [Metabolic Syndrome Definition](#-metabolic-syndrome-definition)
* [Tools and Software](#-tools-and-software)
* [Data Cleaning and Preparation](#-data-cleaning-and-preparation)
* [Statistical Analysis](#-statistical-analysis)
* [Repository Structure](#-repository-structure)
* [Getting Started](#-getting-started)
* [Reproducibility](#-reproducibility)
* [System Requirements](#-system-requirements)
* [Figures and Tables](#-figures-and-tables)
* [Data Availability](#-data-availability)
* [Ethical Approval](#-ethical-approval)
* [Strengths](#-strengths)
* [Limitations](#-limitations)
* [Future Work](#-future-work)
* [Citation](#-citation)
* [Repository Citation](#-repository-citation)
* [License](#-license)
* [Contact](#-contact)
* [Status](#-status)

---

## 🔬 Project Overview

This repository contains the R analysis pipeline used to investigate the longitudinal epidemiology of metabolic syndrome (MetS) among people living with HIV participating in the Virological Impact of Switching from Efavirenz/Nevirapine-based first-line ART to Dolutegravir (VISEND) trial in Zambia.

The analysis evaluates changes in metabolic syndrome over 144 weeks of follow-up and examines associations with antiretroviral therapy regimen, follow-up time, and selected demographic and clinical characteristics.

The primary analysis script is:

```text
MetSyn_Longitudinal_Analysis.R
```

---

## 📌 Study and Data Source

The analysis uses longitudinal data from the VISEND trial, a randomized clinical trial conducted in Zambia evaluating antiretroviral therapy regimens containing dolutegravir.

The analysis includes repeated measurements collected during follow-up through Week 144.

Key variables include:

* Demographic characteristics
* Anthropometric measurements
* Blood pressure
* Fasting glucose
* Lipid measurements
* HIV viral load
* CD4 cell count
* Antiretroviral therapy regimen

The ART regimen groups evaluated in the analysis are:

* **PI-based** – protease inhibitor-based control regimen
* **TAFED** – tenofovir alafenamide/emtricitabine/dolutegravir
* **TLD** – tenofovir/lamivudine/dolutegravir

---

## 🎯 Analysis Objectives

The analysis pipeline was developed to:

1. Describe the longitudinal prevalence of metabolic syndrome over 144 weeks.
2. Evaluate the association between ART regimen and metabolic syndrome.
3. Assess changes in metabolic syndrome risk over follow-up.
4. Investigate non-linear temporal trends in metabolic syndrome.
5. Evaluate associations between demographic and clinical characteristics and metabolic syndrome.
6. Examine changes in individual metabolic syndrome components over time.
7. Conduct sensitivity and supplementary analyses to assess the robustness of the findings.

---

## 🧬 Metabolic Syndrome Definition

Metabolic syndrome was defined using the study-specific criteria applied in the VISEND analysis.

The metabolic syndrome components included:

* **Abdominal obesity**
* **Elevated fasting plasma glucose**
* **Reduced high-density lipoprotein cholesterol (HDL-C)**
* **Elevated triglycerides**
* **Elevated blood pressure**

Participants meeting the required number of abnormal components were classified as having metabolic syndrome.

Participants with metabolic syndrome at baseline were excluded from the longitudinal analysis of incident metabolic syndrome.

---

## 🛠️ Tools and Software

The analysis was conducted using **R**.

### R Packages

The pipeline uses the following R packages:

#### Data manipulation and cleaning
- `dplyr`
- `tidyr`
- `forcats`
- `janitor`
- `purrr`
- `stringr`

#### Statistical modelling
- `geepack`
- `emmeans`
- `ggeffects`
- `margins`
- `epitools`
- `splines`

#### Model tidying and results
- `broom`
- `knitr`
- `gtsummary`
- `gt`

#### Data visualization
- `ggplot2`
- `scales`
- `ggpubr`
- `gridExtra`
- `patchwork`

#### Publication and document outputs
- `flextable`
- `officer`
- `gtExtras`
- `sysfonts`
- `webshot`

The analysis script loads these packages at the beginning of the pipeline and uses them for data preparation, statistical modelling, visualisation, and generation of publication-ready outputs.

---

## 🧹 Data Cleaning and Preparation

The pipeline includes several data preparation steps:

* Importing and inspecting the VISEND longitudinal dataset
* Standardising study visit names
* Identifying and excluding participants with baseline metabolic syndrome
* Recoding ART regimen categories
* Standardising demographic and clinical variables
* Deriving metabolic syndrome and its individual components
* Preparing repeated longitudinal observations
* Creating variables required for regression modelling
* Handling missing observations for specific analyses
* Preparing datasets for descriptive and longitudinal statistical analyses

---

## 📊 Statistical Analysis

The pipeline includes descriptive and longitudinal statistical analyses.

### Descriptive Analysis

The analysis describes:

* Participant characteristics
* Metabolic syndrome prevalence over follow-up
* ART regimen-specific prevalence
* Changes in metabolic syndrome components
* Clinical and demographic characteristics associated with metabolic syndrome

### Longitudinal Analysis

Generalised estimating equations (**GEE**) are used to account for repeated observations from participants over time.

The primary longitudinal models evaluate:

* ART regimen
* Follow-up time
* Sex
* Age
* CD4 cell count

An exchangeable correlation structure is used to account for within-participant correlation.

### Non-linear Time Trends

The pipeline also evaluates whether the relationship between follow-up time and metabolic syndrome is non-linear.

Additional analyses include:

* Spline-based GEE models
* Piecewise GEE models
* Regimen-specific temporal trends
* Component-specific analyses

### Sensitivity Analyses

Additional analyses examine the robustness of the primary findings, including analyses incorporating relevant baseline and treatment-related characteristics.

---

## 📁 Repository Structure

The repository is organised as follows:

```text
VISEND-MetSyn-Longitudinal-Analysis/
│
├── MetSyn_Longitudinal_Analysis.R
├── README.md
├── R/
├── data/
├── figures/
├── tables/
└── results/
```

The participant-level dataset is **not included** in the public repository.

Generated figures, tables, and other analysis outputs may be organised in their respective directories.

---

## 🚀 Getting Started

### Clone the Repository

```bash
git clone https://github.com/MpanjiSiwingwa/VISEND-MetSyn-Longitudinal-Analysis.git
cd VISEND-MetSyn-Longitudinal-Analysis
```

### Install Required R Packages

Open R or RStudio and install the required packages:

```r
install.packages(c(
  "dplyr",
  "tidyr",
  "forcats",
  "ggplot2",
  "geepack",
  "emmeans",
  "ggeffects",
  "broom",
  "gtsummary",
  "gt",
  "flextable",
  "officer",
  "epitools"
))
```

### Run the Analysis

The primary analysis script can be run using:

```r
source("MetSyn_Longitudinal_Analysis.R")
```

Alternatively, from the terminal:

```bash
Rscript MetSyn_Longitudinal_Analysis.R
```

Before running the pipeline, the required dataset must be available in the appropriate local data location.

---

## 🔁 Reproducibility

The analysis uses a fixed random seed where applicable:

```r
set.seed(123)
```

The computational environment can be documented using:

```r
sessionInfo()
```

For projects using `renv`, the environment can be restored with:

```r
renv::restore()
```

Researchers reproducing the analysis should record their R version and package versions because changes in software environments may affect computational results.

---

## 💻 System Requirements

Recommended environment:

* **R ≥ 4.5.1**
* macOS, Linux, or Windows
* At least **8 GB RAM**
* Sufficient storage for the study dataset and generated outputs

Actual computational requirements will depend on the size of the dataset and the analyses being performed.

---

## 📈 Figures and Tables

The analysis pipeline generates figures and tables describing:

* Participant follow-up
* Metabolic syndrome prevalence over time
* ART regimen-specific trends
* Non-linear temporal patterns
* Adjusted associations between ART regimen and metabolic syndrome
* Predicted probabilities of metabolic syndrome
* Individual metabolic syndrome components
* Sensitivity and supplementary analyses

Final manuscript figures and tables should be generated from the version of the analysis script corresponding to the manuscript submission.

---

## 🔐 Data Availability

The VISEND clinical dataset contains sensitive participant information and is **not publicly available** in this repository because of ethical, privacy, and institutional data-sharing restrictions.

The repository therefore contains the **analysis code but not participant-level data**.

Researchers interested in accessing de-identified data for scientific collaboration should contact the relevant VISEND study investigators and obtain the required institutional and ethical approvals.

---

## ⚖️ Ethical Approval

The VISEND study received ethical approval from the **University of Zambia Biomedical Research Ethics Committee (UNZABREC)**.

All participants provided informed consent prior to enrolment in the study.

---

## 💪 Strengths

Key strengths of the analysis include:

* Longitudinal follow-up extending to **144 weeks**
* Repeated metabolic measurements during follow-up
* Evaluation of multiple ART regimen groups
* Assessment of both metabolic syndrome and its individual components
* Use of longitudinal statistical methods accounting for repeated observations
* Evaluation of non-linear temporal trends
* Focus on an African population that remains underrepresented in longitudinal cardiometabolic research

---

## ⚠️ Limitations

Important limitations include:

* Loss to follow-up over the 144-week study period may introduce attrition-related bias.
* Residual confounding from factors not captured in the analysis cannot be excluded.
* The study population may not be representative of all people living with HIV in Zambia or other African populations.
* Findings from this study should be interpreted within the context of the VISEND trial population and study design.

---

## 🔮 Future Work

Future research may build on this analysis by:

* Investigating biological mechanisms underlying metabolic changes during ART
* Evaluating detailed body composition and adiposity measures
* Examining genetic and pharmacogenomic determinants of metabolic response
* Investigating longer-term cardiovascular outcomes
* Evaluating integrated HIV and non-communicable disease care models
* Conducting collaborative analyses across African HIV cohorts
* Developing predictive models for individual cardiometabolic risk

---

## 📚 References

1. Venter WDF, Moorhouse M, Sokhela S, Fairlie L, Mashabane N, Masenya M, et al. Dolutegravir plus Two Different Prodrugs of Tenofovir to Treat HIV. *New England Journal of Medicine*. 2019;381(9):803–815.

2. Hurbans N, Naidoo P. Efficacy, safety, and tolerability of dolutegravir-based ART regimen in Durban, South Africa: a cohort study. *BMC Infectious Diseases*. 2024;24(1).

3. Zambia Consolidated Guidelines for Treatment and Prevention of HIV Infection. 2020.

4. Gebremedhin T, Ayenalem M, Adem M, Geremew D, Aleka Y, Kiflie A. Dolutegravir based therapy showed CD4+ T cell count recovery and viral load suppression among ART-naive HIV-positive individuals: a longitudinal evaluation. *Systematic Reviews in Pharmacy*. 2023;14(4):264–271.

---

## 📌 Citation

If you use this code, please cite:

> Siwingwa M. *VISEND Metabolic Syndrome Longitudinal Analysis*. Version 1.0.0. Zenodo. 2026. https://doi.org/10.5281/zenodo.23023581

### Associated manuscript

Siwingwa M, et al. *Longitudinal Changes in the Odds of Metabolic Syndrome Among People Living with HIV Receiving Dolutegravir-Based Antiretroviral Therapy: A 144-Week Analysis from the VISEND Trial in Zambia.* Manuscript in preparation, 2026.

---

## 📖 Repository Citation

Siwingwa M. *VISEND-MetSyn-Longitudinal-Analysis: Reproducible R Analysis Pipeline for Longitudinal Assessment of Metabolic Syndrome Among People Living with HIV in the VISEND Trial, Zambia.* GitHub repository. 2026.

**Repository:**
https://github.com/MpanjiSiwingwa/VISEND-MetSyn-Longitudinal-Analysis

---

## 📜 License

This project is licensed under the **MIT License**.

---

## 📬 Contact

**Mpanji Siwingwa**
PhD Researcher | Bioinformatics | HIV Research

* **GitHub:** https://github.com/MpanjiSiwingwa
* **Email:** [mpanjisiwingwa@gmail.com](mailto:mpanjisiwingwa@gmail.com)
* **LinkedIn:** https://linkedin.com/in/mpanji-siwingwa-b0272a74
* **ORCID:** https://orcid.org/0000-0002-3623-2108

---

## 📌 Status

This repository contains the analysis pipeline associated with the VISEND longitudinal metabolic syndrome analysis.

The repository may be updated as the manuscript undergoes further review and as additional validation and supplementary analyses are completed.
