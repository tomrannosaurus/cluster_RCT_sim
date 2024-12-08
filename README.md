README

# Resource Optimization in a Cluster Randomized Control Trial

## Background: 
Researchers conducting cluster randomized control trials (cRCTs), including genetic sequencing studies where a sample is tested multiple times, must navigate a trade-off between the number of clusters and their size to optimize their budgets while producing the most meaningful results possible. In this simulation study, we explored cost optimization in a cRCT of the treatment effect on outcome Y. Borrowing the framework of a genetic sequencing study, we conceptualized clusters as individuals and samples within clusters as technical replicates. We examined the effect of changing the cost of adding individuals relative to the cost of adding technical replicates, as well as the effect of changing variance between individuals ($\gamma$) and variance of a given technical replicate from an individual’s true value ($\sigma$). 

## Methods
We used the ADEMP framework to examine which optimal combinations of cluster size and cluster number minimize standard error (SE) of the treatment effect ($\beta$) of outcome variable Y. We performed both normal and Poisson distributions. Moreover, in addition to simulations varying one parameter at a time, we performed factorial designs varying multiple parameters. 

## Results 

For a given problem space, different scenarios varying cost and variance parameters at various levels were produced and analyzed. The directionality of the effects of these parameters on the optimal number of individuals and technical replicates were as expected. Varying the underlying parameters for the data generating process also had predictable effects on our ability to estimate treatment effect.  

When cost ratio was high, the specific value of $\sigma$ and $\gamma$ changed the optimal number of individuals and technical replicates more than when cost ratio was low. Moreover, in high cost ratio scenarios, the optimal r and G to minimize SE of the estimate changes in a non-linear way when varying ($\sigma$ and $\gamma$).

## Conclusion 

*Note - i am not sure that everything captures the difference between univariate and factorial results/conclusions*

*Note - is the term 'variance' being used correctly? should it be changed to ICC*

Overall, the univariate cases yield generally expected results, but factorial design revealed unexpectedly complicated ones. Cost ratio appeared to be more consequential than $\sigma$ and $\gamma$ with respect to the optimal number of individuals and number of technical replicates, at least within the values of $\sigma$ and $\gamma$ that we tested. 
Notably, the optimal number of individuals and number of technical replicates depended on variance in a much more complex way than expected, with specific values of variance changing optimal cluster size and number in a non-uniform, non-linear way. The unexpected complexity of this relationship suggests that estimating variance between and within clusters prior to beginning a cRCT is critical to minimizing standard error. 

Limitations included boundary effects and the volatility of generalized linear models. [*note - optional - as well as lack of justification for range*] Future directions for research could be including cluster size and number combinations that come very close to minimizing SE or using a Bayesian modeling approach.


[The full report can be found here](**add URL**)


## Setup

- R Version: 4.3.1
- Package Versions:
   - tidyverse: 2.0.0
   - knitr: 1.45
   - kableExtra: 1.3.4
   - ggplot2: 3.4.3
   - naniar 1.0.0
   - visdat 0.6.0
   - car 3.1-2
   - lme4 1.1-34
   - ggpubr 0.6.0
   - outliers: 0.15
   - reshape2: 1.4.4
   - moments: 0.14.1
   - vcd: 1.4-12
   - glmnet 4.1-8
   - mice 3.16.0
   - caret 6.0-94
   - pROC 1.18.4

## Files

- Report Quarto: [*add quarto*].qmd
- Helper functions code: _helpers.R
- Tex files: column-commands.tex, float-setup.tex, geometry-settings.tex, table-packages.tex, title-settings.tex
- References file: references.bib