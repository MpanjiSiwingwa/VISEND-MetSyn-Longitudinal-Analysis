# VISEND-MetSyn-Longitudinal-Analysis

**Reproducible Analysis Pipeline for Longitudinal Assessment of Metabolic Syndrome Risk in People Living with HIV on Dolutegravir-Based Therapy**

![R Version](https://img.shields.io/badge/R-4.5.1-blue?logo=r)
![Machine Learning](https://img.shields.io/badge/Machine_Learning-Pipeline-success?logo=r)
![Clinical Research](https://img.shields.io/badge/Clinical-HIV%20Research-red?logo=medrxiv)
![Project Status](https://img.shields.io/badge/Status-Active_Research-orange?logo=github)
![License: MIT](https://img.shields.io/badge/License-MIT-green?logo=open-source-initiative)
[![GitHub Profile](https://img.shields.io/badge/GitHub-MpanjiSiwingwa-black?logo=github)](https://github.com/MpanjiSiwingwa)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.20145268-blue?logo=zenodo)](https://doi.org/10.5281/zenodo.20145268)

---

## 📑 Table of Contents
- [Project Overview](#-project-overview)
- [Data Source](#-data-source)
- [Tools](#-tools)
- [Data Cleaning and Preparation](#-data-cleaning-and-preparation)
- [Exploratory Data Analysis](#-exploratory-data-analysis)
- [Data Analysis](#-data-analysis)
- [Getting Started](#-getting-started)
- [Reproducible Environment](#-reproducible-environment)
- [System Requirements](#-system-requirements)
- [Approximate Runtime](#-approximate-runtime)
- [Figures](#-figures)
- [Results](#-results)
- [Recommendations](#-recommendations)
- [Strengths](#-strengths)
- [Limitations](#-limitations)
- [Future Work](#-future-work)
- [References](#-references)
- [Citation](#-citation)
- [Repository Citation](#-repository-citation)
- [License](#-license)
- [Contact](#-contact)
- [Status](#-status)

---

## 📌 Project Overview
This project investigates the correlation between **Metabolic Syndrome (MetS)** and the use of **dolutegravir (DTG)-based antiretroviral regimens** in people living with HIV (PLHIV). Using longitudinal data from the VISEND clinical study, we assessed the prevalence and risk of MetS over 144 weeks, focusing on regimens such as **TDF/3TC/DTG (TLD)** and **TAF/XTC/DTG (TAFED)**.

---

## 📂 Data Source
The primary dataset was derived from the **VISEND trial**, specifically the file *VISEND_week_48.csv*, which contains 48‑week follow‑up data including:
- Anthropometric measurements  
- Lipid profiles  
- HIV viral load  
- CD4 count  

---

## 🛠️ Tools
- **[Excel](https://www.microsoft.com)** – Data capture and preliminary cleaning  
- **[REDCap](https://redcap.moh.gov.zm)** – Secure data storage and management  
- **[Python](https://www.python.org/downloads)** – Advanced data cleaning, statistical analysis, and visualization  

---

## 🧹 Data Cleaning and Preparation
1. Importing and inspecting raw data.  
2. Handling missing values and inconsistencies.  
3. Standardizing variable formats for analysis.  

---

## 🔍 Exploratory Data Analysis
Exploratory analysis focused on:
- Estimating overall MetS prevalence.  
- Comparing prevalence across DTG-based regimens.  
- Assessing whether all DTG regimens exhibit similar risk patterns.  

---

## 📈 Data Analysis
Analytical methods included regression modeling and stratified comparisons to evaluate associations between MetS and ART regimens.

```python
# Example placeholder for regression analysis
import statsmodels.api as sm
model = sm.GEE.from_formula("MetS ~ Regimen_Type + Age + Sex + CD4_count",
                            groups="Patient_ID", data=df)
result = model.fit()
print(result.summary())

---

## 🛠️ Getting Started

### 🚀 Clone the Repository

```bash
git clone https://github.com/MpanjiSiwingwa/Metabolic-syndrome-prediction-using-machine-learning.git
cd Metabolic-syndrome-prediction-using-machine-learning
Rscript pipeline_metabolic_syndrome_MLA.R
```

---

### ⚙️ Install Dependencies

```r
install.packages(c(
  # 📦 Core Data Manipulation & Cleaning
  "dplyr", "tidyr", "forcats", "janitor",
  
  # 📊 Visualization & Plot Formatting
  "ggplot2", "scales", "gtExtras", "sysfonts", "gridExtra", "ggpubr",
  
  # 📈 Modeling & Marginal Effects
  "geepack", "emmeans", "ggeffects", "margins",
  
  # 🧹 Model Tidying & Reporting
  "broom", "knitr", "gt", "gtsummary", "flextable", "officer",
  
  # 🧪 Epidemiological Tools
  "epitools"
))

))
```

---

### ▶️ Run the Pipeline

```r
source("pipeline_metabolic_syndrome_MLA.R")
```

Alternatively:

```bash
Rscript pipeline_metabolic_syndrome_MLA.R
```

---

## 📦 Reproducible Environment

Package versions were managed using:
- R 4.5.1
- `renv`
- `sessionInfo()`

Restore the computational environment using:

```r
renv::restore()
```

Export session information:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

---

## 💻 System Requirements

- R ≥ 4.5.1
- macOS, Linux, or Windows
- Recommended RAM: ≥8 GB
- Multi-core CPU recommended for model training

---

## ⏱️ Approximate Runtime

| Step | Estimated Runtime |
|---|---|
| Data cleaning | 2–5 min |
| Feature selection | 5–15 min |
| Model training | 20–60 min |
| SHAP analysis | 10–30 min |

---

## 📊 Figures

Key manuscript figures are available in the `figures/` directory.

### Main Figures
- Figure 1. Participant flow and incident metabolic syndrome cases over 144 weeks
- Figure 2. Cumulative incidence of metabolic syndrome over 144 weeks by ART regimen
- Figure 3. Adjusted odds ratios for metabolic syndrome stratified by ART regimen
- Figure 4. Predicted marginal probabilities of metabolic syndrome by ART regimen
- Figure 5. Sex-stratified predicted probabilities of metabolic syndrome by ART regimen
- Figure 6. Adjusted odds ratios for metabolic syndrome from primary and sensitivity analyses

### Supplementary Figures
- xxxxxx
---

### Results
---

- Median age was 44 years (IQR 38–51); 58.4% of participants were female.  
- Metabolic syndrome prevalence increased over follow‑up, peaking at week 120 (TLD: 28.9%; TAFED: 27.4%; PI‑control: 20.0%).  
- DTG‑based regimens were associated with higher MetS risk compared with PI‑based controls (TAFED RR 1.35; TLD RR 1.38).  
- Associations were stronger among women (TAFED RR 1.59; TLD RR 1.66).  
- Older age and higher CD4 count independently predicted increased risk.  


### Recommendations
---

- Integrate **routine metabolic screening** into HIV care, particularly for patients receiving TLD regimens.  
- Conduct **baseline and periodic assessments** of weight, waist circumference, lipid profiles, and glucose metabolism to enable early detection and management of metabolic complications.  
- Develop **sex- and age-specific interventions**, given the heightened risk among women and older adults.  
- Carefully balance **ART regimen selection** between virological efficacy and metabolic safety, especially in patients with pre-existing cardiometabolic risk factors.  
- Implement **lifestyle modification programs**, including culturally adapted dietary counselling, physical activity promotion, and smoking cessation, tailored to African contexts.  
- Strengthen **laboratory capacity** for metabolic monitoring and train healthcare providers in integrated HIV–NCD management.  
- Establish **referral systems** for specialized care to manage complex metabolic complications.  
- Include **metabolic endpoints** in ART program evaluations to monitor long-term cardiovascular outcomes.  
- Support **policy-level initiatives** to ensure sustainability of integrated HIV–NCD care models in resource-limited settings.  
- Advance **pharmacogenomics research** and collaborative efforts to build African genomic databases, enabling precision medicine approaches to metabolic risk in HIV care.  

---

## 🔐 Data Availability

The VISEND clinical dataset used in this study contains sensitive participant information and is not publicly available due to ethical and institutional restrictions.

Researchers interested in accessing de-identified data for scientific collaboration may contact the corresponding author subject to institutional approvals and data-sharing agreements.

All scripts required to reproduce the analyses are fully available in this repository.

---

## 🧾 Ethical Approval

The VISEND study received ethical approval from the University of Zambia Biomedical Research Ethics Committee (UNZABREC). All participants provided informed consent prior to enrollment.

---

## 📑 Reporting Standards

This repository and accompanying manuscript were developed in alignment with:
- TRIPOD reporting recommendations
- TRIPOD-AI guidance principles
- Transparent and reproducible machine learning practices in clinical research

---

## 💪 Strengths
- Large sample size, enhancing statistical power and reliability of findings.  
- Longitudinal design with extended 144‑week follow‑up, allowing assessment of long‑term metabolic outcomes.  
- Comprehensive clinical and laboratory assessments conducted at multiple time points.  
- Use of standardized metabolic syndrome criteria, improving comparability with other studies.  
- Among the few investigations to report long‑term metabolic outcomes of DTG‑based ART in a large African cohort.  
- Provides valuable context‑specific evidence relevant to HIV care in sub‑Saharan Africa.  


---

## ⚠️ Limitations

- The observational design precludes causal inference, and residual confounding from unmeasured factors such as diet, physical activity, or genetic predisposition cannot be excluded.  
- Attrition over 144 weeks may have introduced selection bias if participants lost to follow-up differed in their metabolic risk.  
- The study was conducted among urban Zambian adults, which may limit generalizability to rural settings or other African populations with different demographic and epidemiological profiles.  

---

## 🔮 Future Work

- Explore mechanistic pathways linking **DTG** to metabolic disturbances, including inflammatory markers, adipokines, and detailed body composition analyses.  
- Conduct **randomized controlled trials** comparing DTG with alternative ART regimens over extended follow-up in African populations to establish causal evidence.  
- Implement and evaluate **integrated HIV–NCD care models** that are feasible and scalable in resource-limited settings.  
- Investigate the role of **pharmacogenomics**, as genetic predisposition may explain interindividual variability in metabolic responses to DTG.  
- Support collaborative efforts to establish **African genomic databases** to advance precision medicine and metabolic research in HIV care.  


---

## 📚 References

1. 	Venter WDF, Moorhouse M, Sokhela S, Fairlie L, Mashabane N, Masenya M, et al. Dolutegravir plus Two Different Prodrugs of Tenofovir to Treat HIV. New England Journal of Medicine. 2019 Aug 29;381(9):803–15. 
2. 	Hurbans N, Naidoo P. Efficacy, safety, and tolerability of dolutegravir-based ART regimen in Durban, South Africa: a cohort study. BMC Infect Dis. 2024 Dec 1;24(1). 
3. 	Zambia Consolidated Guidelines for Treatment and Prevention of HIV Infection. 2020. 
4. 	Gebremedhin T, Ayenalem M, Adem M, Geremew D, Aleka Y, Kiflie A. Dolutegravir Based Therapy Showed CD4 + T Cell Count Recovery and Viral Load Suppression among Art Naive HIV Positive Individuals: A Longitudinal Evaluation. Systematic Review Pharmacy. 2023;14(4):264–71.

 ---

## 📌 Citation

Siwingwa M, et al. *Longitudinal Assessment of Metabolic Syndrome Risk in People Living with HIV on Dolutegravir-Based Antiretroviral Therapy: A 144-Week Analysis from the VISEND Trial in Zambia.* Manuscript in preparation, 2026.

---

## 📖 Repository Citation

Siwingwa M. *Metabolic Syndrome Prediction Using Machine Learning in HIV Cohorts Receiving Dolutegravir-Based ART* [GitHub repository]. 2026.

Available at:  
https://github.com/MpanjiSiwingwa/Metabolic-syndrome-prediction-using-machine-learning

---

## 📜 License

This project is licensed under the MIT License.

---

## 📬 Contact

**Mpanji Siwingwa**  
PhD Researcher | Machine Learning | HIV Research | Bioinformatics

- GitHub: https://github.com/MpanjiSiwingwa
- Email: mpanjisiwingwa@gmail.com
- LinkedIn: https://linkedin.com/in/mpanji-siwingwa-b0272a74
- ORCID: https://orcid.org/0000-0002-3623-2108

---

## 📌 Status

This repository accompanies an ongoing research project and will continue evolving as additional analyses and validation studies are completed.
