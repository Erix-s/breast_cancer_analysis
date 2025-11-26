Proposal
================

## Team Member:

Eric Xu <sx2402@cumc.columbia.edu>  
Alexandria Taylor <at3583@cumc.columbia.edu>  
Kevin Guzman <kg3006@cumc.columbia.edu>  
Edmund Wong <ecw2178@cumc.columbia.edu>  
Jiseung Lee <jl6458@cumc.columbia.edu>

## Proposal(draft):

*This project aims to investigate survival outcomes in a pressing
women’s health issue, emphasizing the demographic and biological factors
that shape disease prognosis and health equity.*

**Primary Topic: Breast Cancer Survival and Risk Patterns**  
Breast cancer remains one of the most prevalent disease worldwide,
accounting for a substantial proportion of cancer incidence and health
burden. Despite advances in early detection and treatment, disparities
in outcomes is often observed across demographic and biological groups.
Clinical outcomes are influenced by multiple factors, including tumor
size, stage, lymph node involvement, and hormone receptor status.

Understanding how patient characteristics, tumor biology, and receptor
profiles influence survival and risk is essential for improving risk
stratification and guiding treatment decisions. This project aims to
explore the relationship between the following factors predicting
survival outcomes among patients with breast cancer.

*This project starts with a public, de-identified dataset derived from
the 2017 SEER (Surveillance, Epidemiology, and End Results) Program,
containing information from 4,024 breast cancer cases with 16
demographic, clinical, and molecular variables. The data include
survival time, tumor stage and grade, hormone receptor status, lymph
node involvement, and other key prognostic factors.*

In further investigation, usage of this dataset or a full SEER dataset
with more detailed covariates would be considered. The main goal is to
explore survival differences by demographic group, tumor
characteristics, and receptor status, as well as broader examination of
disease burden and disparities in breast cancer outcomes.

The project may also incorporate publicly available SEER incidence data
(1975–2022) to contextualize survival findings within trends in breast
cancer occurrence and population-level risk patterns over time.

*If SEER data set is not applicable or not informative, an alternative
plan would be applied.*

**Alternative Topic: Infant and Maternal Mortality in the United
States**  
As an alternative or complementary direction, this project may
investigate patterns and disparities in infant and maternal mortality
using national data from the CDC/NCHS Birth–Infant Death and Maternal
Mortality datasets. These datasets provide detailed demographic and
health information on mothers and infants across the United States,
offering a rich foundation for understanding the social, economic, and
medical determinants of survival in early life and maternal health.

This topic is of major public health importance, as both infant and
maternal mortality reflect the overall effectiveness of healthcare
systems, access to care, and equity in health outcomes. Analyses could
focus on identifying demographic patterns, geographic disparities, and
population subgroups at elevated risk, contributing to ongoing national
efforts to reduce preventable deaths and improve maternal–child health
outcomes.

**Project Vision**

Both directions align closely with the goals of applied epidemiologic
and data science research. The overarching aim is to use population data
to describe health disparities across female health burden and
systematic determinants. With this project, we aim to identify key
determinants of disease outcomes, and highlight areas for intervention.
By integrating data management, exploratory insight, and public health
relevance, this project will contribute to understanding critical
patterns in population health and inform future work in epidemiology and
health equity.

## Final Product Plans

Our team will produce a fully reproducible, Github-hosted project that
identify predictors of survival in breast cancer. Deliverables will
include.

- **Final Written Report**: A comprehensive report detailing workflow
  following the required structure (Motivation, Related Work, Data,
  Exploratory Analysis, Additional Analysis, Discussion).

- **Interactive Project Website**: Built using R Markdown, hosted via
  GitHub Pages, and design to communicate results to a broad audience.

- **Two-minute Video Summary**: A concise video presentation summarizing
  key findings, methodology, and public health implications of the
  project that is embedded into the project webpage.

- **GitHub Repository**: A well-organized repository containing all
  code, data processing scripts, and documentation to ensure full
  reproducibility of the analysis.

## Planned Analysis/Visualization/Coding Challanges

### Planned Analysis:

- Descriptive statistic and distribution of demographics and clincial
  variables

- Correlation and association test between tumor characteristics/lymph
  node characteristics/hormone receptor status and survival outcomes

- Kaplan-Meier survival analysis

- Cox proportional hazards regression

### Visualizations:

- Kaplan-Meier survival curves stratified by key variables (e.g., tumor
  stage, hormone receptor status)

- Boxplots of tumor size by stage

- Heatmap of variables correlation

### Coding Challenges:

- Data cleaning and preprocessing of SEER dataset (if able to be)

- Cleaning categorial variables in publicly available data set

- Reformating data for survival analysis

- Managing missing or inconsistent data enteries.

## Project Timeline

| Week / Date | Milestone | Key Tasks |
|:---|:---|:---|
| Week 1 (Nov 11–17) | Data Cleaning & Setup | Import data, recode variables, handle missing values, set up<br>GitHub repo & R/Quarto structure. |
| Week 2 (Nov 18–24) | Exploratory Analysis & Visualization | Descriptives & plots; correlations; refine questions; draft<br>report sections. |
| Week 3 (Nov 25–Dec 1) | Modeling & Interpretation | Kaplan–Meier & Cox models; test associations; interpret<br>results; finalize figures & tables. |
| Week 4 (Dec 2–Dec 6) | Final Deliverables & Submission | Finish report & website; record 2-min screencast; polish<br>repo; submit by Dec 6 (11:59 PM). |
