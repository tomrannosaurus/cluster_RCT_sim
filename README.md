# Resource Optimization in a Cluster Randomized Control Trial

## Background
Cluster randomized control trials (cRCTs) require careful balance between the number of clusters and observations per cluster to optimize limited research budgets. Using the framework of a genetic sequencing study, where clusters represent individuals and within-cluster observations represent technical replicates, we explored how to optimize this balance through simulation. We investigated how optimal designs vary with the relative costs of adding new clusters versus new observations within clusters, while accounting for between-cluster variation ($\gamma$) and within-cluster variation ($\sigma$).

## Methods
Using the ADEMP framework, we simulated cRCTs under both normal and Poisson distributions to identify configurations that minimize the standard error (SE) of the treatment effect ($\beta$). Our approach included:
- Univariate analyses examining individual parameter effects
- Factorial designs exploring parameter interactions
- Cost ratios ($c_1$/$c_2$) ranging from 20:19 to 100:1
- Various ICC scenarios through different $\gamma$ and $\sigma$ combinations
- Budget constraints enforcing trade-offs between cluster number (G) and replicates (R)

## Key Findings

### Univariate Analysis
- Demonstrated predictable relationships between individual parameters and optimal designs
- Revealed initial evidence of boundary conditions where optimal configurations change abruptly

### Factorial Analysis
- Uncovered complex interactions between cost structures and variance components
- Showed that cost ratio effects dominate variance effects within tested parameter ranges
- Identified similar optimization patterns between normal and Poisson cases, with Poisson showing more defined optimal regions
- Demonstrated that optimal designs cluster more distinctly in high-variance scenarios

### Main Insights
1. Higher cost ratios ($c_1$/$c_2$) lead to:
   - Greater sensitivity to variance parameters
   - More volatile optimal configurations
   - Non-linear relationships between variance and optimal designs

2. ICC effects:
   - More pronounced in high total variance scenarios
   - Create consistent monotonic relationships with SE in Poisson cases
   - Generate variable relationships in normal cases

## Conclusions
While individual parameter effects follow expected patterns, their interactions create complex optimization landscapes. Cost ratios primarily drive optimal design choices, but variance parameters introduce important nuances. The study emphasizes the critical importance of accurate variance estimation before trial implementation, as suboptimal design choices may significantly impact precision.

## Limitations
- Boundary effects in optimization
- Model convergence challenges in some scenarios
- Limited exploration of variance parameter ranges
- Focus on specific cost ratio ranges


[The full report can be found here](https://github.com/tomrannosaurus/cluster_RCT_sim/blob/main/cluster_RCT_sim.pdf)


## Setup

- R Version: 4.3.1
- Package Versions:
   - tidyverse: 2.0.0
   - ggplot2: 3.4.3
   - lme4 1.1-34
   - ggrepel 0.9.6
   - latex2exp 0.9.6
   - parallel 4.3.1
   - doParallel 1.0.17

## Files

- Report Quarto: cluster_RCT_sim.qmd
- Tex files: column-commands.tex, float-setup.tex, geometry-settings.tex, table-packages.tex, title-settings.tex
- References file: references.bib