## ----setup, include = FALSE, message = FALSE--------------------------------------------
# --- Preamble ---
# Date of last update: Dec. 8, 2024
# R Version: 4.3.1
# Package Versions:
#   tidyverse: 2.0.0
#   ggplot2: 3.4.3
#   lme4 1.1-34
#   ggrepel 0.9.6
#   latex2exp 0.9.6
#   parallel 4.3.1
#   doParallel 1.0.17


setwd("~/GitHub/cluster_RCT_sim")

# Knitr Engine Setup
knitr::opts_chunk$set(message=F, 
                      warning=F, 
                      error=F, 
                      echo=F, 
                      fig.pos = "H" ,
                      fig.align = 'center')

# Packages
options(kableExtra.latex.load_packages = FALSE) # Required to avoid floatrow error
library(knitr)
library(kableExtra)
library(ggplot2)
library(ggrepel)
library(latex2exp)
library(lme4)
library(tidyverse)
library(parallel)   
library(doParallel) 

# seed
set.seed(42)




## ----core code, eval = F----------------------------------------------------------------
## 
## 
## ######### Section: Simulation Code #########
## 
## ####### Subsection: DGP #######
## 
## generate_clustered_data <- function(
## #' Generate Clustered Gaussian Data
## #'
## #' Simulates clustered data with a Gaussian outcome, including random effects
## #' at the cluster level and measurement error at the observation level.
## #'
## #' @param n_clusters Num of clusters (min 2). Ensures treatment variation
## #' across clusters
## #' @param n_replicates Num of replicates per cluster
## #' @param alpha Intercept of the model
## #' @param beta Treatment effect on the outcome
## #' @param gamma Standard deviation of the cluster-level random effects
## #' @param sigma Standard deviation of the measurement error
## #' @param seed Optional. Random seed for reproducibility
## #'
## #' @return A data frame containing simulated data, including observed outcome
## #' and associated covariates
##     n_clusters = 100, # num clusters (G)
##     n_replicates = 3, # num of replicates per cluster (R)
##     alpha = 0, # intercept
##     beta = 0.5, # treatment effect
##     gamma = 1, # SD of cluster random effects
##     sigma = 0.5, # SD of measurement error
##     seed = NULL
## ) {
##   if (!is.null(seed)) set.seed(42)
## 
##   if (n_clusters < 2) {  #BUGFIX
##     stop("Number of clusters must be at least 2 to ensure treatment variation.")
##   }
## 
##   clusters <- 1:n_clusters
## 
##   # Make sure at least one treated and one control cluster BUGFIX
##   treatment <- rep(c(0, 1), length.out = n_clusters)
##   treatment <- sample(treatment)  # shuffled to randomize treatment assignment
## 
##   epsilon <- rnorm(n_clusters, 0, gamma) # cluster-level random effects
## 
##   data <- expand.grid(
##     cluster = clusters,
##     replicate = 1:n_replicates
##   ) %>%
##     mutate(
##       treatment = treatment[cluster],
##       epsilon = epsilon[cluster],
##       mu = alpha + beta * treatment + epsilon, # true mean
##       error = rnorm(n(), 0, sigma), # measurment error
##       Y = mu + error #obs outcome
##     )
## 
##   return(data)
## }
## 
## 
## 
## generate_clustered_poisson_data <- function(
## #' Generate Clustered Poisson Data
## #'
## #' Simulates clustered data with a Poisson-distributed outcome, including random
## #' effects at the cluster level affecting the log mean of the distribution.
## #'
## #' @param n_clusters Number of clusters (minimum 2). Ensures treatment variation
## #' across clusters
## #' @param n_replicates Number of replicates per cluster
## #' @param alpha Intercept of the log mean model
## #' @param beta Treatment effect on the log mean of the outcome
## #' @param gamma SD of the cluster-level random effects on the log scale
## #' @param seed Optional. Random seed for reproducibility
## #'
## #' @return A data frame containing simulated data, including observed outcome
## #' and associated covariates
##     n_clusters = 100,
##     n_replicates = 3,
##     alpha = 0,
##     beta = 0.5,
##     gamma = 0.5,
##     seed = NULL
## ) {
##   if (!is.null(seed)) set.seed(42)
## 
##   if (n_clusters < 2) {
##     stop("Number of clusters must be at least 2 to ensure treatment variation.")
##   }
## 
##   clusters <- 1:n_clusters
## 
##   # Make sure at least one treated and one control cluster BUGFIX
##   treatment <- rep(c(0, 1), length.out = n_clusters)
##   treatment <- sample(treatment)  # shuffled to randomize treatment assignment
## 
##   epsilon <- rnorm(n_clusters, 0, gamma) # cluster-level random eff on log scale
## 
##   data <- expand.grid(
##     cluster = clusters,
##     replicate = 1:n_replicates
##   ) %>%
##     mutate(
##       treatment = treatment[cluster],
##       epsilon = epsilon[cluster],
##       log_mu = alpha + beta * treatment + epsilon,
##       mu = exp(log_mu),
##       Y = rpois(n(), lambda = mu)
##     )
## 
##   return(data)
## }
## 
## 
## 
## 
## 
## fit_mixed_model <- function(data, family = "gaussian") {
## #' Fit Mixed Model
## #'
## #' Fits a mixed-effects model or standard regression model based on the specified
## #' family and data structure. Supports Gaussian and Poisson families.
## #'
## #' @param data Input dataset containing outcome and covariates
## #' @param family Model family, either "gaussian" or "poisson"
## #'
## #' @return A list containing estimated treatment effect, standard error, confidence
## #' interval, the fitted model object, and a flag indicating convergence
##   default_return <- list(
##     beta_hat = NA_real_,
##     beta_se = NA_real_,
##     conf_int = c(NA_real_, NA_real_),
##     model = NULL,
##     converged = FALSE
##   )
## 
##   if (family == "gaussian") {
##     if (all(table(data$cluster) == 1)) { # use lm when R is 1
##       model <- lm(Y ~ treatment, data = data)
##       converged <- TRUE
##     } else {
##       model <- lmer(Y ~ treatment + (1 | cluster),
##                     data = data,
##                     control = lmerControl(optimizer = "bobyqa",
##                                           calc.derivs = FALSE,
##                                           optCtrl = list(maxfun = 1e5,
##                                                          rel.tol = 1e-5)))
##       converged <- TRUE
##     }
## 
##   } else if (family == "poisson") {
##     if (all(table(data$cluster) == 1)) {
##       model <- glm(Y ~ treatment, family = poisson, data = data)
##       converged <- TRUE
##     } else {
##       model <- glmer(Y ~ treatment + (1 | cluster),
##                      family = poisson,
##                      data = data,
##                      nAGQ = 0,
##                      control = glmerControl(optimizer = "bobyqa",
##                                             calc.derivs = FALSE,
##                                             optCtrl = list(maxfun = 1e5,
##                                                            rel.tol = 1e-5)))
##       converged <- TRUE
##     }
##   }
## 
##   if (!converged || is.null(model)) {
##     return(default_return)
##   }
## 
##   # grab results
##   if (inherits(model, "lm") || inherits(model, "glm")) {
##     beta_hat <- coef(model)["treatment"]
##     beta_se <- summary(model)$coefficients["treatment", "Std. Error"]
##     conf_int <- tryCatch(confint(model)["treatment", ],
##                          error = function(e) c(NA_real_, NA_real_))
##   } else {
##     beta_hat <- fixef(model)["treatment"]
##     beta_se <- sqrt(vcov(model)["treatment", "treatment"])
##     conf_int <- tryCatch(confint(model, method = "Wald")["treatment", ],
##                          error = function(e) c(NA_real_, NA_real_))
##   }
## 
##   list(
##     beta_hat = beta_hat,
##     beta_se = beta_se,
##     conf_int = conf_int,
##     model = model,
##     converged = TRUE
##   )
## }
## 
## 
## 
## 
## 
## 
## 
## 
## 
## get_feasible_designs <- function(
## #' Get Feasible Designs
## #'
## #' Identifies feasible combinations of clusters and replicates within a given
## #' budget, considering fixed and additional per-replicate costs.
## #'
## #' @param budget Total budget available
## #' @param cost_first Cost per cluster for the first replicate
## #' @param cost_additional Cost per additional replicate per cluster
## #' @param min_clusters Min num of clusters (default 7, to ensure model fit)
## #' @param max_replicates Max num of replicates per cluster
## #'
## #' @return A data frame containing feasible designs with columns: num of clusters,
## #' num of replicates, and total cost
##     budget, # B
##     cost_first, # c1
##     cost_additional, # c2
##     min_clusters = 7, #(note: models tend to fail to fit below 7)
##     max_replicates = 200
## ) {
## 
##   # init storage
##   G_values <- c()
##   R_values <- c()
##   total_costs <- c()
## 
##   # step 1: maximum possible G given B and c1
##   max_possible_G <- floor(budget / cost_first)
## 
##   # step 2: loop to test all G for possible R
##   for(g in min_clusters:max_possible_G) {
##     # remaining budget after G
##     budget_remaining <- budget - (g * cost_first)
## 
##     # find R that fits into budget given G
##     R_max <- floor(budget_remaining / (g * cost_additional)) + 1
##     R_max <- min(R_max, max_replicates)
## 
##     if(R_max >= 1) { # BUGFIX ensuring R > 0
##       # calc total cost
##       total_cost <- g * cost_first + g * (R_max - 1) * cost_additional
## 
##       # double check total cost < budget
##       if(total_cost <= budget) {
##         G_values <- c(G_values, g)
##         R_values <- c(R_values, R_max)
##         total_costs <- c(total_costs, total_cost)
##       }
##     }
##   }
##   designs <- data.frame(G = G_values,
##                         R = R_values,
##                         cost = total_costs)
## 
##   return(designs)
## }
## 
## 
## 
## 
## ####### Subsection: Run Simulations #######
## 
## ## ---## ---## ---## ---
## # setup parallel
## cl <- makeCluster(detectCores())
## clusterEvalQ(cl, {
##   library(lme4)
##   library(tidyverse)
## })
## 
## clusterExport(cl, varlist = c("generate_clustered_data",
##                               "generate_clustered_poisson_data",
##                               "fit_mixed_model"),
##               envir = environment())
## ## ---## ---## ---## ---
## 
## evaluate_design <- function(
## #' Evaluate Design
## #'
## #' Simulates data for a given design (clusters and replicates) and evaluates
## #' model performance using specified parameters and family.
## #'
## #' @param G Num of clusters
## #' @param R Num of replicates per cluster
## #' @param n_sims Num of simulations to perform
## #' @param alpha Intercept of the model
## #' @param beta Treatment effect
## #' @param gamma SD of cluster-level random effects
## #' @param sigma SD of measurement error (Gaussian family only)
## #' @param family Model family, either "gaussian" or "poisson"
## #' @param cl Optional parallel cluster object for parallel simulations
## #'
## #' @return A list containing mean beta estimate, empirical SE, mean model SE,
##     G,
##     R,
##     n_sims = 100,
##     alpha = 1,
##     beta = 0.5,
##     gamma = 1,
##     sigma = 0.5,
##     family = "gaussian",
##     cl = NULL) {
## 
##   # Calculate expected mu for Poisson case
##   # Using moment generating function of normal for epsilon
##   if(family == "poisson") {
##     control_mu = exp(alpha + gamma^2/2)  # X=0 case
##     treat_mu = exp(alpha + beta + gamma^2/2)  # X=1 case
##     avg_mu = (control_mu + treat_mu)/2
##   }
##   # parallel handling (under 10 parallel not worth overhead)
##   run_parallel <- !is.null(cl) && n_sims > 10
## 
##   # BUGFIX: local variables to avoid conflicts
##   local_G <- G
##   local_R <- R
##   local_alpha <- alpha
##   local_beta <- beta
##   local_gamma <- gamma
##   local_sigma <- sigma
##   local_family <- family
## 
## 
##   sim_func <- function(i) {
##     if (local_family == "gaussian") {
##       sim_data <- generate_clustered_data(
##         n_clusters = local_G,
##         n_replicates = local_R,
##         alpha = local_alpha,
##         beta = local_beta,
##         gamma = local_gamma,
##         sigma = local_sigma
##       )
##     } else if (local_family == "poisson") {
##       sim_data <- generate_clustered_poisson_data(
##         n_clusters = local_G,
##         n_replicates = local_R,
##         alpha = local_alpha,
##         beta = local_beta,
##         gamma = local_gamma
##       )
##     }
## 
##     results <- fit_mixed_model(sim_data, family = local_family)
##     list(beta_hat = results$beta_hat, beta_se = results$beta_se)
##   }
## 
## 
##   if (run_parallel) {
##     results_list <- parLapply(cl, 1:n_sims, sim_func)
##   } else {
##     results_list <- lapply(1:n_sims, sim_func)
##   }
## 
##   beta_hats <- sapply(results_list, `[[`, "beta_hat")
##   ses <- sapply(results_list, `[[`, "beta_se")
## 
##   # return different values based on family
##   base_results <- list(
##     mean_beta_hat = mean(beta_hats, na.rm = TRUE),
##     empirical_se = sd(beta_hats, na.rm = TRUE),
##     mean_model_se = mean(ses, na.rm = TRUE),
##     coverage = mean(abs(beta_hats - beta) <= 1.96 * ses, na.rm = TRUE)
##   )
## 
##   if(family == "gaussian") {
##     base_results$measurement_var <- sigma
##   } else {
##     base_results$measurement_var <- sqrt(avg_mu)
##   }
## 
##   return(base_results)
## }
## 
## 
## 
## 
## 
## 
## find_optimal_design <- function(
## #' Find Optimal Design
## #'
## #' Identifies the optimal design (clusters and replicates) under a specified
## #' budget by evaluating feasible designs through simulation.
## #'
## #' @param budget Total budget available
## #' @param cost_first Cost per cluster for the first replicate
## #' @param cost_additional Cost per additional replicate per cluster
## #' @param n_sims Num of simulations per design
## #' @param parameters List of model parameters
## #' @param family Model family, either "gaussian" or "poisson"
## #' @param cl Optional parallel cluster object for parallel simulations
## #'
## #' @return A list with two components: all_results (data frame of
## #' evaluated designs) and optimal (design with lowest empirical SE)
##     budget,
##     cost_first,
##     cost_additional,
##     n_sims = 100,
##     parameters = list(alpha = 1, beta = 0.5, gamma = 1, sigma = 0.5),
##     family = "gaussian",
##     cl = NULL) {
## 
##   designs <- get_feasible_designs(budget, cost_first, cost_additional)
## 
##   results <- list()
##   for(i in 1:nrow(designs)) {
##     cat(sprintf("Evaluating design %d of %d: G=%d, R=%d\n",
##                 i, nrow(designs), designs$G[i], designs$R[i]))
## 
##     eval <- evaluate_design(
##       G = designs$G[i],
##       R = designs$R[i],
##       n_sims = n_sims,
##       alpha = parameters$alpha,
##       beta = parameters$beta,
##       gamma = parameters$gamma,
##       sigma = parameters$sigma,
##       family = family,
##       cl = cl
##     )
## 
##     results[[i]] <- c(designs[i,],
##                       empirical_se = eval$empirical_se,
##                       model_se = eval$mean_model_se,
##                       coverage = eval$coverage,
##                       measurement_var = eval$measurement_var)
##   }
## 
##   results_df <- do.call(rbind.data.frame, results)
##   optimal <- results_df[which.min(results_df$empirical_se),]
## 
##   return(list(
##     all_results = results_df,
##     optimal = optimal
##   ))
## }
## 
## 
## 
## 
## ####### Subsection: Factorial Simulation #######
## sim_factorial <- function(
## #' Simulate Factorial Combinations
## #'
## #' Simulates and evaluates optimal designs across factorial combinations of
## #' cost, measurement variance, and cluster-level variance for a given budget.
## #'
## #' @param budget Total budget available
## #' @param sigma_range Range of measurement variances (for Gaussian family)
## #' @param gamma_range Range of SDs for cluster-level random effects
## #' @param c1_range Range of costs per cluster for the first replicate
## #' @param c2_range Range of costs per additional replicate per cluster
## #' @param family Model family, either "gaussian" or "poisson"
## #' @param alpha Intercept of the model
## #' @param beta Treatment effect
## #' @param n_sims Num of simulations per design
## #' @param cl Optional parallel cluster object for parallel simulations
## #'
## #' @return A list with two components: optimal_results (data frame of optimal
## #' designs for each combination) and all_scenarios_results (data frame of results
## #' for all designs across all combinations)
##     budget = 1000,
##     sigma_range = seq(0.2, 1, by = 0.2),
##     gamma_range = seq(0.5, 2, by = 0.3),
##     c1_range = c(20, 50, 100),
##     c2_range = c(1, 5, 10),
##     family = "gaussian",
##     alpha = 1,
##     beta = 0.5,
##     n_sims = 1000,
##     cl = NULL
## ) {
##   if(family == "gaussian") {
##     param_grid <- expand.grid(
##       measurement_var = sigma_range,
##       gamma = gamma_range,
##       c1 = c1_range,
##       c2 = c2_range
##     )
##   } else {
##     # for Poisson case, compute mu range based on gamma
##     compute_mu <- function(gamma) {
##       control_mu <- exp(alpha + gamma^2/2)  # X=0 case
##       treat_mu <- exp(alpha + beta + gamma^2/2)  # X=1 case
##       (control_mu + treat_mu)/2  # Average mu
##     }
## 
##     # mu values for each gamma
##     mu_values <- sapply(gamma_range, compute_mu)
## 
##     param_grid <- expand.grid(
##       measurement_var = sqrt(mu_values),  # sqrt to match scale of sigma
##       gamma = gamma_range,
##       c1 = c1_range,
##       c2 = c2_range
##     )
##   }
##   # add computed vars
##   param_grid$cost_ratio <- param_grid$c1 / param_grid$c2
## 
##   # ICC differently for Poisson vs normal
##   if(family == "gaussian") {
##     param_grid$ICC <- (param_grid$gamma^2) /
##                      (param_grid$gamma^2 + param_grid$measurement_var^2)
##   } else {
##     param_grid$ICC <- param_grid$gamma^2 /
##                      (param_grid$gamma^2 + log(param_grid$measurement_var^2))
##   }
## 
## 
##   # init overall storage obj
##   results_optimal_list <- vector("list", nrow(param_grid))
##   results_all_list <- vector("list", nrow(param_grid))
## 
##   for(i in 1:nrow(param_grid)) {
##     cat(sprintf("\nTesting combination %d of %d:\n", i, nrow(param_grid)))
## 
##     # Adjust printing based on family
##     if(family == "gaussian") {
##       cat(sprintf("sigma=%.2f, gamma=%.2f, c1=%d, c2=%d \n",
##                   param_grid$measurement_var[i], param_grid$gamma[i],
##                   param_grid$c1[i], param_grid$c2[i]))
##     } else {
##       cat(sprintf("sqrt(mu)=%.2f, gamma=%.2f, c1=%d, c2=%d \n",
##                   param_grid$measurement_var[i], param_grid$gamma[i],
##                   param_grid$c1[i], param_grid$c2[i]))
##     }
## 
##     parameters <- list(
##       alpha = alpha,
##       beta = beta,
##       gamma = param_grid$gamma[i],
##       sigma = if(family == "gaussian") param_grid$measurement_var[i] else NULL
##     )
## 
##     optimal <- find_optimal_design(
##       budget = budget,
##       cost_first = param_grid$c1[i],
##       cost_additional = param_grid$c2[i],
##       n_sims = n_sims,
##       parameters = parameters,
##       family = family,
##       cl = cl
##     )
## 
##     # save optimal results with appropriate naming
##     results_optimal_list[[i]] <- data.frame(
##       measurement_var = param_grid$measurement_var[i],
##       gamma = param_grid$gamma[i],
##       c1 = param_grid$c1[i],
##       c2 = param_grid$c2[i],
##       cost_ratio = param_grid$cost_ratio[i],
##       ICC = param_grid$ICC[i],
##       optimal_G = if(!is.null(optimal$optimal$G)) optimal$optimal$G else NA_real_,
##       optimal_R = if(!is.null(optimal$optimal$R)) optimal$optimal$R else NA_real_,
##       empirical_se = if(!is.null(optimal$optimal$empirical_se)) optimal$optimal$empirical_se else NA_real_,
##       model_se = if(!is.null(optimal$optimal$model_se)) optimal$optimal$model_se else NA_real_,
##       coverage = if(!is.null(optimal$optimal$coverage)) optimal$optimal$coverage else NA_real_,
##       total_cost = if(!is.null(optimal$optimal$cost)) optimal$optimal$cost else NA_real_
##     )
## 
##     # save all results
##     all_res_scenario <- optimal$all_results
##     all_res_scenario$measurement_var <- param_grid$measurement_var[i]
##     all_res_scenario$gamma <- param_grid$gamma[i]
##     all_res_scenario$c1 <- param_grid$c1[i]
##     all_res_scenario$c2 <- param_grid$c2[i]
##     all_res_scenario$cost_ratio <- param_grid$cost_ratio[i]
##     all_res_scenario$ICC <- param_grid$ICC[i]
##     results_all_list[[i]] <- all_res_scenario
##   }
## 
##   optimal_res <- do.call(rbind, results_optimal_list)
##   all_res <- do.call(rbind, results_all_list)
## 
##   return(list(
##     optimal_results = optimal_res,
##     all_scenarios_results = all_res
##   ))
## }
## 
## 
## 
## 
## 
## ####### Subsection: Univariate Simulation #######
## 
## compute_design_cost <- function(G, R, c1, c2) {
##   total_cost <- G * c1 + G * (R - 1) * c2
##   return(total_cost)
## }
## 
## # NOTE: Normal case only, would need to alter fixed_params for Poisson
## explore_parameter_effect <- function(
## #' Explore Parameter Effect
## #'
## #' Conducts univariate simulations to assess the effect of varying the number of
## #' clusters G or replicates R while keeping other parameters fixed.
## #' Optionally considers budget constraints and adjusts designs accordingly.
## #'
## #' @param param_name Name of the parameter to vary ("G" or "R")
## #' @param param_range Range of values for the parameter
## #' @param fixed_params List of fixed parameter values for all other variables
## #' @param budget Optional total budget to constrain designs
## #' @param cost_first Cost per cluster for the first replicate
## #' @param cost_additional Cost per additional replicate per cluster
## #' @param n_sims Num of simulations to perform for each parameter value
## #' @param family Model family, either "gaussian" or "poisson" NOTE: BROKEN
## #' @param cl Optional parallel cluster object for parallel simulations
## #'
## #' @return A data frame containing simulation results for each parameter value,
## #' including variance estimate, empirical SE, model SE, coverage, and design
## #' characteristics (G, R, total_cost, and budget status)
##   param_name,
##   param_range,
##   fixed_params = list(
##     G = 20,
##     R = 5,
##     alpha = 1,
##     beta = 0.5,
##     gamma = 1,
##     sigma = 0.5
##   ),
##   budget = NULL,
##   cost_first = NULL,
##   cost_additional = NULL,
##   n_sims = 100,
##   family = "gaussian",
##   cl = NULL
## ) {
##   results <- data.frame(
##     param_value = param_range,
##     var_estimate = NA,
##     empirical_se = NA,
##     model_se = NA,
##     coverage = NA,
##     G_used = NA,
##     R_used = NA,
##     total_cost = NA,
##     budget = budget,
##     cost_first = cost_first,
##     cost_additional = cost_additional,
##     is_maximal = NA
##   )
## 
##   for(i in seq_along(param_range)) {
##     current_params <- fixed_params
##     current_params[[param_name]] <- param_range[i]
## 
##     if(!is.null(budget)) {
##       if(param_name == "G") {
##         # when G varies, calc max R
##         R_val <- 1 + floor((budget - param_range[i] * cost_first) /
##                              (param_range[i] * cost_additional))
##         R_val <- max(1, R_val)
##         current_params$R <- R_val
##         results$R_used[i] <- R_val
##         results$G_used[i] <- param_range[i]
##       } else if(param_name == "R") {
##         # when R varies, calc max G
##         G_val <- floor(budget / (cost_first + (param_range[i] - 1) * cost_additional))
##         G_val <- max(2, G_val)  # least 2 clusters BUGFIX
##         current_params$G <- G_val
##         results$G_used[i] <- G_val
##         results$R_used[i] <- param_range[i]
##       }
## 
##       #current total cost / budget check
##       results$total_cost[i] <- compute_design_cost(results$G_used[i],
##                                                    results$R_used[i],
##                                                    cost_first,
##                                                    cost_additional)
##       if(results$total_cost[i] > budget) {
##         results$is_maximal[i] <- "over_budget"
##       } else {
##         # BUGFIX: budget boundary check
##         cost_plus_G <- compute_design_cost(results$G_used[i] + 1,
##                                            results$R_used[i],
##                                            cost_first,
##                                            cost_additional)
##         cost_plus_R <- compute_design_cost(results$G_used[i],
##                                            results$R_used[i] + 1,
##                                            cost_first,
##                                            cost_additional)
## 
##         results$is_maximal[i] <- if((cost_plus_G > budget) || (cost_plus_R > budget)) {
##           "maximal" } else { "non_maximal" }
##       }
##     }
## 
##     sim_results <- evaluate_design(
##       G = current_params$G,
##       R = current_params$R,
##       n_sims = n_sims,
##       alpha = current_params$alpha,
##       beta = current_params$beta,
##       gamma = current_params$gamma,
##       sigma = current_params$sigma,
##       family = family,
##       cl = cl
##     )
## 
##     results$var_estimate[i] <- sim_results$empirical_se^2
##     results$empirical_se[i] <- sim_results$empirical_se
##     results$model_se[i] <- sim_results$mean_model_se
##     results$coverage[i] <- sim_results$coverage
## 
##   }
## 
##   return(results)
## }
## 
## 
## 
## # NOTE: Normal case only, would need to alter fixed_params for Poisson
## explore_model_parameter <- function(
## #' Explore Model Parameter
## #'
## #' Simulates and evaluates the effect of varying a model parameter
## #' (alpha, beta, sigma, or gamma) while keeping others fixed.
## #'
## #' @param param_name Name of the parameter to vary
## #' @param param_range Range of values for the parameter
## #' @param fixed_params List of fixed parameter values
## #' @param budget Total budget for design
## #' @param cost_first Cost per cluster for the first replicate
## #' @param cost_additional Cost per additional replicate per cluster
## #' @param n_sims Num of simulations to perform
## #' @param family Model family, either "gaussian" or "poisson" NOTE: BROKEN
## #' @param cl Optional parallel cluster object
## #'
## #' @return Data frame of results for each parameter value, including
## #' variance estimate, empirical SE, model SE, coverage, and budget status
##     param_name,
##     param_range,
##     fixed_params = list(
##       G = 20,
##       R = 5,
##       alpha = 1,
##       beta = 0.5,
##       gamma = 1,
##       sigma = 0.5
##     ),
##     budget = 1000,
##     cost_first = 50,
##     cost_additional = 2,
##     n_sims = 100,
##     family = "gaussian",
##     cl = NULL
## ) {
##   results <- data.frame(
##     param_value = param_range,
##     var_estimate = NA,
##     empirical_se = NA,
##     model_se = NA,
##     coverage = NA,
##     G_used = fixed_params$G,
##     R_used = fixed_params$R,
##     total_cost = NA,
##     budget = budget,
##     cost_first = cost_first,
##     cost_additional = cost_additional,
##     is_maximal = NA
##   )
## 
##   for(i in seq_along(param_range)) {
##     current_params <- fixed_params
##     current_params[[param_name]] <- param_range[i]
## 
##     #current total cost / budget check
##     results$total_cost[i] <- compute_design_cost(fixed_params$G,
##                                                  fixed_params$R,
##                                                  cost_first,
##                                                  cost_additional)
## 
##     if(results$total_cost[i] > budget) {
##       results$is_maximal[i] <- "over_budget"
##     } else {
##       # BUGFIX: budget boundary check
##       cost_plus_G <- compute_design_cost(fixed_params$G + 1,
##                                          fixed_params$R,
##                                          cost_first,
##                                          cost_additional)
##       cost_plus_R <- compute_design_cost(fixed_params$G,
##                                          fixed_params$R + 1,
##                                          cost_first,
##                                          cost_additional)
## 
##       results$is_maximal[i] <- if((cost_plus_G > budget) || (cost_plus_R > budget)) {
##         "maximal" } else { "non_maximal" }
##     }
## 
##     sim_results <- evaluate_design(
##       G = fixed_params$G,
##       R = fixed_params$R,
##       n_sims = n_sims,
##       alpha = current_params$alpha,
##       beta = current_params$beta,
##       gamma = current_params$gamma,
##       sigma = current_params$sigma,
##       family = family,
##       cl = cl
##     )
## 
##     results$var_estimate[i] <- sim_results$empirical_se^2
##     results$empirical_se[i] <- sim_results$empirical_se
##     results$model_se[i] <- sim_results$mean_model_se
##     results$coverage[i] <- sim_results$coverage
##   }
## 
##   return(results)
## }
## 
## # NOTE: Normal case only, would need to alter fixed_params for Poisson
## explore_cost_parameter <- function(
## #' Explore Cost Parameter
## #'
## #' Simulates and evaluates the effect of varying cost parameters (c1 or c2)
## #' while keeping design and model parameters fixed.
## #'
## #' @param param_name Name of the cost parameter to vary (c1 or c2)
## #' @param param_range Range of values for the parameter
## #' @param fixed_params List of fixed parameter values
## #' @param budget Total budget for design
## #' @param other_cost Value of the other cost parameter
## #' @param n_sims Num of simulations to perform
## #' @param family Model family, either "gaussian" or "poisson" NOTE: BROKEN
## #' @param cl Optional parallel cluster object
## #'
## #' @return Data frame of results for each parameter value, including
## #' variance estimate, empirical SE, model SE, coverage, and budget status
##     param_name,
##     param_range,
##     fixed_params = list(
##       G = 20,
##       R = 5,
##       alpha = 1,
##       beta = 0.5,
##       gamma = 1,
##       sigma = 0.5
##     ),
##     budget = 1000,
##     other_cost = NULL,
##     n_sims = 100,
##     family = "gaussian",
##     cl = NULL
## ) {
##   results <- data.frame(
##     param_value = param_range,
##     var_estimate = NA,
##     empirical_se = NA,
##     model_se = NA,
##     coverage = NA,
##     G_used = fixed_params$G,
##     R_used = fixed_params$R,
##     total_cost = NA,
##     budget = budget,
##     is_maximal = NA
##   )
## 
##   for(i in seq_along(param_range)) {
##     if(param_name == "c1") {
##       cost_first <- param_range[i]
##       cost_additional <- other_cost
##     } else {
##       cost_first <- other_cost
##       cost_additional <- param_range[i]
##     }
## 
##     #current total cost / budget check
##     results$total_cost[i] <- compute_design_cost(fixed_params$G,
##                                                  fixed_params$R,
##                                                  cost_first,
##                                                  cost_additional)
## 
##     if(results$total_cost[i] > budget) {
##       results$is_maximal[i] <- "over_budget"
##     } else {
##       # BUGFIX: budget boundary check
##       cost_plus_G <- compute_design_cost(fixed_params$G + 1,
##                                          fixed_params$R,
##                                          cost_first,
##                                          cost_additional)
##       cost_plus_R <- compute_design_cost(fixed_params$G,
##                                          fixed_params$R + 1,
##                                          cost_first,
##                                          cost_additional)
## 
##       results$is_maximal[i] <- if((cost_plus_G > budget) || (cost_plus_R > budget)) {
##         "maximal" } else { "non_maximal" }
##     }
## 
##     sim_results <- evaluate_design(
##       G = fixed_params$G,
##       R = fixed_params$R,
##       n_sims = n_sims,
##       alpha = fixed_params$alpha,
##       beta = fixed_params$beta,
##       gamma = fixed_params$gamma,
##       sigma = fixed_params$sigma,
##       family = family,
##       cl = cl
##     )
## 
##     results$var_estimate[i] <- sim_results$empirical_se^2
##     results$empirical_se[i] <- sim_results$empirical_se
##     results$model_se[i] <- sim_results$mean_model_se
##     results$coverage[i] <- sim_results$coverage
##   }
## 
##   return(results)
## }
## 
## 
## 


## ----eval = F---------------------------------------------------------------------------
## ######### Section: Simulation Runs #########
## 
## ####### Subsection: Univariate Runs #######
## 
## # G exploration
## res_G_gam_under_sig <- explore_parameter_effect(
##   param_name = "G",
##   param_range = seq(2, 100, by = 1),
##   fixed_params = list(
##     alpha = 1,
##     beta = 0.5,
##     gamma = 0.5,
##     sigma = 1
##   ),
##   budget = 1000,
##   cost_first = 100,
##   cost_additional = 1,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_G_gam_under_sig, "sim_data/res_G_gam_under_sig.rds")


## ---------------------------------------------------------------------------------------
res_G_gam_under_sig <- readRDS("sim_data/res_G_gam_under_sig.rds")


## ----eval = F---------------------------------------------------------------------------
## res_G_gam_over_sig <- explore_parameter_effect(
##   param_name = "G",
##   param_range = seq(2, 100, by = 1),
##   fixed_params = list(
##     alpha = 1,
##     beta = 0.5,
##     gamma = 2,
##     sigma = 1
##   ),
##   budget = 1000,
##   cost_first = 100,
##   cost_additional = 1,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_G_gam_over_sig, "sim_data/res_G_gam_over_sig.rds")


## ---------------------------------------------------------------------------------------
res_G_gam_over_sig <- readRDS("sim_data/res_G_gam_over_sig.rds")


## ----eval = F---------------------------------------------------------------------------
## # For R exploration
## res_R_gam_under_sig <- explore_parameter_effect(
##   param_name = "R",
##   param_range = seq(1, 240, by = 3),
##   fixed_params = list(
##     alpha = 1,
##     beta = 0.5,
##     gamma = 0.5,
##     sigma = 1
##   ),
##   budget = 1000,
##   cost_first = 50,
##   cost_additional = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_R_gam_under_sig, "sim_data/res_R_gam_under_sig.rds")

## ---------------------------------------------------------------------------------------
res_R_gam_under_sig <- readRDS("sim_data/res_R_gam_under_sig.rds")


## ----eval = F---------------------------------------------------------------------------
## # For R exploration
## res_R_gam_over_sig <- explore_parameter_effect(
##   param_name = "R",
##   param_range = seq(1, 240, by = 3),
##   fixed_params = list(
##     alpha = 1,
##     beta = 0.5,
##     gamma = 2,
##     sigma = 1
##   ),
##   budget = 1000,
##   cost_first = 50,
##   cost_additional = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_R_gam_over_sig, "sim_data/res_R_gam_over_sig.rds")


## ---------------------------------------------------------------------------------------
res_R_gam_over_sig <- readRDS("sim_data/res_R_gam_over_sig.rds")


## ----eval = F---------------------------------------------------------------------------
## # For gamma parameter (not used in report)
## res_gamma <- explore_model_parameter(
##   param_name = "gamma",
##   param_range = seq(0.1, 3, by = 0.2),
##   fixed_params = list(
##     G = 25,
##     R = 3,
##     alpha = 1,
##     beta = 0.5,
##     sigma = 0.5
##   ),
##   budget = 1000,
##   cost_first = 50,
##   cost_additional = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_gamma, "sim_data/res_gamma.rds")


## ---------------------------------------------------------------------------------------
res_gamma <- readRDS("sim_data/res_gamma.rds")


## ----eval = F---------------------------------------------------------------------------
## # Sigma parameter (not used in report)
## res_sigma <- explore_model_parameter(
##   param_name = "sigma",
##   param_range = seq(0.1, 2, by = 0.1),
##   fixed_params = list(
##     G = 20,
##     R = 5,
##     alpha = 1,
##     beta = 0.5,
##     gamma = 1,
##   ),
##   budget = 1000,
##   cost_first = 50,
##   cost_additional = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_sigma, "sim_data/res_sigma.rds")


## ---------------------------------------------------------------------------------------
res_sigma <- readRDS("sim_data/res_sigma.rds")


## ----eval = F---------------------------------------------------------------------------
## # alpha parameter  (not used in report)
## res_alpha <- explore_model_parameter(
##   param_name = "alpha",
##   param_range = seq(-2, 2, by = 0.2),
##   fixed_params = list(
##     G = 20,
##     R = 5,
##     beta = 0.5,
##     gamma = 1,
##     sigma = 0.5
##   ),
##   budget = 1000,
##   cost_first = 50,
##   cost_additional = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_alpha, "sim_data/res_alpha.rds")


## ---------------------------------------------------------------------------------------
res_alpha <- readRDS("sim_data/res_alpha.rds")


## ----eval = F---------------------------------------------------------------------------
## # Beta parameter (not used in report)
## res_beta <- explore_model_parameter(
##   param_name = "beta",
##   param_range = seq(0, 2, by = 0.1),
##   fixed_params = list(
##     G = 20,
##     R = 5,
##     alpha = 1,
##     gamma = 1,
##     sigma = 0.5
##   ),
##   budget = 1000,
##   cost_first = 50,
##   cost_additional = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_beta, "sim_data/res_beta.rds")

## ---------------------------------------------------------------------------------------
res_beta <- readRDS("sim_data/res_beta.rds")


## ----eval = F---------------------------------------------------------------------------
## # c1 parameter (not used in report)
## res_c1 <- explore_cost_parameter(
##   param_name = "c1",
##   param_range = seq(10, 200, by = 10),
##   fixed_params = list(
##     G = 20,
##     R = 5,
##     alpha = 1,
##     beta = 0.5,
##     gamma = 1,
##     sigma = 0.5
##   ),
##   budget = 1000,
##   other_cost = 2,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_c1, "sim_data/res_c1.rds")


## ---------------------------------------------------------------------------------------
res_c1 <- readRDS("sim_data/res_c1.rds")


## ----eval = F---------------------------------------------------------------------------
## # c2 parameter (not used in report)
## res_c2 <- explore_cost_parameter(
##   param_name = "c2",
##   param_range = seq(1, 20, by = 1),
##   fixed_params = list(
##     G = 20,
##     R = 5,
##     alpha = 1,
##     beta = 0.5,
##     gamma = 1,
##     sigma = 0.5
##   ),
##   budget = 1000,
##   other_cost = 50,
##   n_sims = 500,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(res_c2, "sim_data/res_c2.rds")

## ---------------------------------------------------------------------------------------
res_c2 <- readRDS("sim_data/res_c2.rds")


## ----eval = F---------------------------------------------------------------------------
## ######### Section: Simulation Runs #########
## 
## ####### Subsection: Factorial Runs #######
## 
## # Run factorial: normal
## resg <- sim_factorial(
##   budget = 1000,
##   sigma_range = seq(0.2, 2, by = 0.3),
##   gamma_range = seq(0.5, 3, by = 0.5),
##   c1_range = c(20, 50, 100),
##   c2_range = c(1, 10, 19),
##   n_sims = 1000,
##   family = "gaussian",
##   cl = cl
## )
## 
## #saveRDS(resg, "sim_data/results_g2.rds")


## ---------------------------------------------------------------------------------------
resg <- readRDS("sim_data/results_g2.rds")


## ----eval = F---------------------------------------------------------------------------
## # Run factorial: poisson
## resp <- sim_factorial(
##   budget = 1000,
##   gamma_range = seq(0.5, 3, by = 0.5),
##   c1_range = c(20, 50, 100),
##   c2_range = c(1, 10, 19),
##   n_sims = 500,
##   family = "poisson",
##   cl = cl
## )
## 
## #saveRDS(results_p, "sim_data/results_p3.rds")


## ---------------------------------------------------------------------------------------
resp <- readRDS("sim_data/results_p_fix2.rds")


## ---------------------------------------------------------------------------------------
######### Section: Simulation Plots #########

####### Subsection: Univariate Plot Functions #######

# NOTE: Normal case only, would need to alter label for Poisson
plot_G_exploration <- function(results, fixed_params, include_se = TRUE) {
#' Plot G Exploration
#'
#' Creates a plot showing empirical SE vs num of clusters (G) for design
#' exploration, highlighting budget constraints and maximal designs.
#'
#' @param results Data frame of results from G exploration
#' @param fixed_params List of fixed parameter values
#' @param include_se Logical, whether to include SE ribbons in the plot
#'
#' @return ggplot object showing empirical SE against num of clusters
  subtitle <- TeX(sprintf(
    "Fixed parameters: $\\alpha=%.2f$,
    $\\beta=%.2f$, $\\gamma=%.2f$, $\\sigma=%.2f$,
    $c_1=%.0f$, $c_2=%.0f$, $B=%.0f$",
    fixed_params$alpha, fixed_params$beta, fixed_params$gamma, fixed_params$sigma,
    results$cost_first, results$cost_additional, results$budget
  ))
  
  
  results$label <- sprintf("SE = %.3f\nG = %d\nR = %d", 
                           results$empirical_se,
                           results$param_value,
                           results$R_used)
  
  p <- ggplot(results, aes(x = param_value)) +
    geom_line(aes(y = empirical_se), color = "grey70") +
    geom_point(aes(y = empirical_se, color = is_maximal), size = 3) +
    scale_color_manual(values = c("black" = "black", 
                                  "non_maximal" = "black",
                                  "maximal" = "red", 
                                  "over_budget" = "grey50"),
                       name = "Budget Usage",
                       labels = c("non_maximal" = "Non-maximal",
                                  "maximal" = "Maximal", 
                                  "over_budget" = "Over Budget"))+
    labs(
      title = "Empirical SE vs Number of Clusters (G)",
      subtitle = subtitle,
      x = "Number of Clusters (G)",
      y = "Empirical Standard Error"
    ) +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, NA))
  
  if(include_se) {
    p <- p + 
      geom_ribbon(
        aes(ymin = pmax(0, empirical_se - 2*model_se),
            ymax = empirical_se + 2*model_se),
        alpha = 0.2
      )
  }
  
  results$to_label <- (seq_len(nrow(results)) - 1) %% 3 == 0
  
  p + geom_text_repel(
    data = subset(results, to_label), 
    aes(y = empirical_se, 
        label = label,
        color = is_maximal),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    force = 2,
    color = "blue",
    show.legend = FALSE
  )
}


# NOTE: Normal case only, would need to alter label for Poisson
plot_R_exploration <- function(results, fixed_params, include_se = TRUE) {
#' Plot R Exploration
#'
#' Creates a plot showing empirical SE vs replicates per cluster (R) for
#' design exploration, highlighting budget constraints and maximal designs.
#'
#' @param results Data frame of results from R exploration
#' @param fixed_params List of fixed parameter values
#' @param include_se Logical, whether to include SE ribbons in the plot
#'
#' @return ggplot object showing empirical SE against replicates per cluster
  subtitle <- TeX(sprintf(
    "Fixed parameters: $\\alpha=%.2f$,
    $\\beta=%.2f$, $\\gamma=%.2f$, $\\sigma=%.2f$,
    $c_1=%.0f$, $c_2=%.0f$, $B=%.0f$",
    fixed_params$alpha, fixed_params$beta, fixed_params$gamma, fixed_params$sigma,
    results$cost_first, results$cost_additional, results$budget
  ))
  
  
  results$label <- sprintf("SE = %.3f\nG = %d\nR = %d", 
                           results$empirical_se,
                           results$G_used,
                           results$param_value)
  
  p <- ggplot(results, aes(x = param_value)) +
    geom_line(aes(y = empirical_se), color = "grey70") +
    geom_point(aes(y = empirical_se, color = is_maximal), size = 3) +
    scale_color_manual(values = c("black" = "black", 
                                  "non_maximal" = "black",
                                  "maximal" = "red", 
                                  "over_budget" = "grey50"),
                       name = "Budget Usage",
                       labels = c("non_maximal" = "Non-maximal",
                                  "maximal" = "Maximal", 
                                  "over_budget" = "Over Budget"))+
    labs(
      title = "Empirical SE vs Replicates per Cluster (R)",
      subtitle = subtitle,
      x = "Replicates per Cluster (R)",
      y = "Empirical Standard Error"
    ) +
    theme_minimal() + theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, NA))
  
  if(include_se) {
    p <- p + 
      geom_ribbon(
        aes(ymin = pmax(0, empirical_se - 2*model_se),
            ymax = empirical_se + 2*model_se),
        alpha = 0.2
      )
  }
  
  results$to_label <- (seq_len(nrow(results)) - 1) %% 3 == 0
  
  p + geom_text_repel(
    data = subset(results, to_label), 
    aes(y = empirical_se, 
        label = label,
        color = is_maximal),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    force = 2,
    color = "blue",
    show.legend = FALSE
  )
}


plot_model_parameter <- function(results, param_name, fixed_params) {
  param_labels <- list(
    alpha = "Intercept (α)",
    beta = "Treatment Effect (β)",
    sigma = if(family == "gaussian") "Measurement Error SD (σ)" else "Expected Count (μ)",
    gamma = "Cluster Effect SD (γ)"
  )
  
  subtitle <- sprintf(
    "Fixed parameters: G=%d, R=%d, budget=%d, c1=%d, c2=%d", 
    results$G_used[1], results$R_used[1], 
    results$budget[1], results$cost_first[1], results$cost_additional[1]
  )
  
  results$label <- sprintf("SE = %.3f", results$empirical_se)
  
  ggplot(results, aes(x = param_value)) +
    geom_line(aes(y = empirical_se), color = "grey70") +
    geom_point(aes(y = empirical_se, color = is_maximal), size = 3) +
    geom_ribbon(aes(ymin = pmax(0, empirical_se - 2*model_se),
                    ymax = empirical_se + 2*model_se),
                alpha = 0.2) +
    scale_color_manual(values = c("maximal" = "red", 
                                  "non_maximal" = "black",
                                  "over_budget" = "grey50")) +
    labs(
      title = sprintf("Empirical SE vs %s", param_labels[[param_name]]),
      subtitle = subtitle,
      x = param_labels[[param_name]],
      y = "Empirical Standard Error"
    ) +
    theme_minimal() +
    scale_y_continuous(limits = c(0, NA)) +
    geom_text_repel(
      data = subset(results, seq_len(nrow(results)) %% 3 == 0),
      aes(y = empirical_se, label = label),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.2,
      force = 2,
      color = "blue",
      show.legend = FALSE
    )
}



# NOTE: Normal case only, would need to alter label for Poisson
plot_cost_parameter <- function(results, param_name, fixed_params) {
  subtitle <- sprintf(
    "Fixed parameters: G=%d, R=%d, α=%.2f, β=%.2f, γ=%.2f, σ=%.2f, B=%d",
    results$G_used[1], results$R_used[1],
    fixed_params$alpha, fixed_params$beta, fixed_params$gamma, fixed_params$sigma,
    results$budget[1]
  )
  
  results$label <- sprintf("SE = %.3f", results$empirical_se)
  
  ggplot(results, aes(x = param_value)) +
    geom_line(aes(y = empirical_se), color = "grey70") +
    geom_point(aes(y = empirical_se, color = is_maximal), size = 3) +
    geom_ribbon(aes(ymin = pmax(0, empirical_se - 2*model_se),
                    ymax = empirical_se + 2*model_se),
                alpha = 0.2) +
    scale_color_manual(values = c("maximal" = "red", 
                                  "non_maximal" = "black",
                                  "over_budget" = "grey50")) +
    labs(
      title = sprintf("Empirical SE vs %s", 
                      ifelse(param_name == "c1", "First Sample Cost", 
                             "Additional Sample Cost")),
      subtitle = subtitle,
      x = ifelse(param_name == "c1", "Cost of First Sample (c1)", 
                 "Cost of Additional Sample (c2)"),
      y = "Empirical Standard Error"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_y_continuous(limits = c(0, NA)) +
    geom_text_repel(
      data = subset(results, seq_len(nrow(results)) %% 3 == 0),
      aes(y = empirical_se, label = label),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.2,
      force = 2,
      color = "blue",
      show.legend = FALSE
    )
}



####### Subsection: Factorial Plot Function #######


create_factorial_analysis_plots <- function(
#' Create Factorial Analysis Plots
#'
#' Generates a set of plots from factorial analysis results, including
#' contour plots for optimal designs, ICC analysis, and cost efficiency.
#'
#' @param results_df Data frame containing factorial analysis results
#' @param budget Total budget for designs
#' @param family Model family, either "gaussian" or "poisson"
#'
#' @return A list containing three ggplot objects (contour plot, ICC plot,
#' efficiency plot) and a summary statistics table    
  
  results_df,
  budget = 1000,
  family = "gaussian") {
  
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(geomtextpath)
  library(ggrepel)
  
  results_df <- results_df %>%
  mutate(
    total_variance = if(family == "gaussian") 
                     gamma^2 + measurement_var^2 
                     else 
                     gamma^2 + log(measurement_var^2),
      cost_combination = paste0("c[1] == ", c1, " ~ c[2] == ", c2)
    ) %>%
    mutate(
      cost_combination = factor(
        cost_combination,
        levels = c("c[1] == 100 ~ c[2] == 1",
                   "c[1] == 100 ~ c[2] == 10",
                   "c[1] == 100 ~ c[2] == 19",
                   "c[1] == 50 ~ c[2] == 1", 
                   "c[1] == 50 ~ c[2] == 10", 
                   "c[1] == 50 ~ c[2] == 19", 
                   "c[1] == 20 ~ c[2] == 1",
                   "c[1] == 20 ~ c[2] == 10",
                   "c[1] == 20 ~ c[2] == 19") 
      )
    )
  
  # Contour plot
  contour_plot <- ggplot(results_df, 
                         aes(x = measurement_var, 
                             y = gamma)) +
    geom_contour_filled(
      aes(z = optimal_G),
      binwidth = 3
    ) +
    scale_fill_viridis_d(name = "Optimal G",
                         option = "magma",
                         guide = guide_legend(
                           position = "bottom",
                           direction = "horizontal", 
                           label.position = "bottom",
                           nrow = 1     
                         )) +
    geom_contour(
      aes(z = optimal_R, color = after_stat(level)),
      binwidth = 3,
      linetype = "dotted",
      linewidth = .15
    ) +
    scale_color_viridis_c(name = "Optimal R",
                          option = "magma",
                          direction = -1) +
    facet_wrap(~ cost_combination, labeller = label_parsed) + 
    labs(
      title = "Optimal Number of Clusters (G) & Replicates (R)",
      subtitle = "Contours show optimal G & R that minimizes empirical SE",
      x = if(family == "poisson") 
        expression(sqrt(mu) ~ "(Measurement Error)") 
      else 
        expression(sigma ~ "(Measurement Error SD)"),
      y = expression(gamma ~ "(Cluster Variation SD)")
    ) +
    theme_minimal()
  
  # ICC analysis plot
  icc_plot <- results_df %>%
    mutate(
      variance_cat = cut(total_variance, 
                         breaks = quantile(total_variance, probs = seq(0, 1, 0.25)),
                         labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"),
                         include.lowest = TRUE)
    ) %>%
    ggplot(aes(x = optimal_G, y = empirical_se, color = ICC)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    facet_wrap(~variance_cat) +
    scale_color_viridis_c(name = "ICC") +
    labs(
      title = "Effect of Number of Clusters on SE | Controlling for Total Variance",
      x = "Number of Clusters (G)",
      y = "Empirical SE"
    ) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 8),
          axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          legend.title = element_text(size = 8)
  )
  
  # Cost efficiency plot
  optimal_points <- results_df %>%
    group_by(cost_combination) %>%
    slice_min(empirical_se) %>%
    ungroup() %>%
    mutate(label = sprintf("SE = %.3f", empirical_se))
  
  efficiency_plot <- ggplot(results_df, 
                            aes(x = optimal_G, y = optimal_R, color = empirical_se)) +
    geom_point(size = 3, alpha = 0.7) + geom_line() +
    geom_point(data = optimal_points, size = 5, shape = 21, 
               color = "black", fill = "red", alpha = 0.7) +
    geom_text_repel(data = optimal_points, 
                    aes(label = label),
                    box.padding = 0.5,
                    point.padding = 0.5,
                    segment.color = 'grey50',
                    size = 3, max.overlaps = Inf) +
    facet_wrap(~ cost_combination, labeller = label_parsed) + 
    scale_color_viridis_c(name = "Empirical SE", option = "cividis") +
    labs(
      title = "Optimal Design Choices by Cost Combination",
      subtitle = "Red points indicate lowest SE combination",
      x = "Number of Clusters (G)",
      y = "Replicates per Cluster (R)"
    ) +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(size = 8), 
          plot.subtitle = element_text(size = 6),
          axis.title.x = element_text(size = 8), 
          axis.title.y = element_text(size = 8),
          legend.title = element_text(size = 8)
  )
  
  #c1 / c2 stats
  summary_stats <- results_df %>%
    group_by(c1, c2) %>%
    summarize(
      mean_G = mean(optimal_G),
      mean_R = mean(optimal_R),
      mean_SE = mean(empirical_se),
      min_SE = min(empirical_se),
      optimal_ICC = ICC[which.min(empirical_se)],
      optimal_total_var = total_variance[which.min(empirical_se)],
      .groups = 'drop'
    )
  
  return(list(
    contour_plot = contour_plot,
    icc_plot = icc_plot,
    efficiency_plot = efficiency_plot,
    summary_stats = summary_stats
  ))
}

# factorial summary table
factorial_summary <- function(results) {
  optimal_points <- results %>%
    group_by(c1, c2) %>%
    slice_min(empirical_se) %>%
    ungroup() %>%
    mutate(total_variance = gamma^2 + sigma^2) %>%
    arrange(empirical_se)
  
  summary_table <- optimal_points %>%
    select(
      c1, c2, cost_ratio, optimal_G, optimal_R, 
      empirical_se, ICC, total_variance
    ) %>%
    arrange(empirical_se)
  
  return(summary_table)
}






## ---------------------------------------------------------------------------------------
####### Subsection: Univariate Plots #######

p_G_gam_over_sig <- plot_G_exploration(res_G_gam_over_sig, fixed_params = list(
  alpha = 1,
  beta = 0.5,
  gamma = 2,
  sigma = 1
))
p_G_gam_under_sig <- plot_G_exploration(res_G_gam_under_sig, fixed_params = list(
  alpha = 1,
  beta = 0.5,
  gamma = 0.5,
  sigma = 1
))



## ----fig.height=3.8---------------------------------------------------------------------
#| fig-cap: "SE vs. G, Fixed Parameters ($\\gamma = 2$)"

p_G_gam_over_sig 


## ----fig.height=3.8---------------------------------------------------------------------
#| fig-cap: "SE vs. G, Fixed Parameters ($\\gamma = .5$)"

p_G_gam_under_sig


## ---------------------------------------------------------------------------------------

p_R_gam_over_sig <- plot_R_exploration(res_R_gam_over_sig, fixed_params = list(
  alpha = 1,
  beta = 0.5,
  gamma = 2,
  sigma = 1
))
p_R_gam_under_sig <- plot_R_exploration(res_R_gam_under_sig, fixed_params = list(
  alpha = 1,
  beta = 0.5,
  gamma = 0.5,
  sigma = 1
))



## ----fig.height=3.8---------------------------------------------------------------------
#| fig-cap: "SE vs. R, Fixed Parameters ($\\gamma = 2$)"

p_R_gam_over_sig 


## ----fig.height=3.8---------------------------------------------------------------------
#| fig-cap: "SE vs. R, Fixed Parameters ($\\gamma = .5$)"

p_R_gam_under_sig


## ----eval = F---------------------------------------------------------------------------
## p_gamma <- plot_model_parameter(res_gamma, "gamma", fixed_params = list(
##   G = 25,
##   R = 3,
##   alpha = 1,
##   beta = 0.5,
##   sigma = 0.5
## ))
## 
## p_sigma <- plot_model_parameter(res_sigma, "sigma", fixed_params = list(
##   G = 20, R = 5, alpha = 1, beta = 0.5, gamma = 1
## ))
## 
## p_alpha <- plot_model_parameter(res_alpha, "alpha", fixed_params = list(
##   beta = 0.5, gamma = 1, sigma = 0.5
## ))
## 
## p_beta <- plot_model_parameter(res_beta, "beta", fixed_params = list(
##   alpha = 1, gamma = 1, sigma = 0.5
## ))
## 
## p_c1 <- plot_cost_parameter(res_c1, "c1", fixed_params = list(
##   alpha = 1, beta = 0.5, gamma = 1, sigma = 0.5
## ))
## 
## p_c2 <- plot_cost_parameter(res_c2, "c2", fixed_params = list(
##   alpha = 1, beta = 0.5, gamma = 1, sigma = 0.5
## ))
## 


## ---------------------------------------------------------------------------------------
plots_g <- create_factorial_analysis_plots(resg$optimal_results)



## ----fig.width=12.5, fig.height=8.5-----------------------------------------------------
#| fig-cap: "Contour Plot: Normal Case"

plots_g$contour_plot


## ----fig.width=8, fig.height=2.8--------------------------------------------------------
#| fig-cap: "ICC Plot: Normal Case"

plots_g$icc_plot


## ----fig.width=8.2, fig.height=3.7------------------------------------------------------
#| fig-cap: "Efficiency Plot: Normal Case"

plots_g$efficiency_plot



## ---------------------------------------------------------------------------------------
plots_p <- create_factorial_analysis_plots(resp$optimal_results, family = "poisson")



## ----fig.width=12.5, fig.height=8.5-----------------------------------------------------
#| fig-cap: "Contour Plot: Poisson Case"

plots_p$contour_plot


## ----fig.width=8, fig.height=2.8--------------------------------------------------------
#| fig-cap: "ICC Plot: Poisson Case"

plots_p$icc_plot


## ----fig.width=8.2, fig.height=3.7------------------------------------------------------
#| fig-cap: "Efficiency Plot: Poisson Case"

plots_p$efficiency_plot



## ----ref.label=knitr::all_labels(), echo=TRUE, eval=FALSE-------------------------------
## NA

