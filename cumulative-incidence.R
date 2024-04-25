# Compute cumulative incidence estimates


# Loading ####

VERBOSE = TRUE
## Functions ####
library('survival')
library('prodlim')
library('mice')
library('parallel')
library('dplyr')
library('geeasy')
library('ggplot2')
library('stringr')

source('landmark_functions.R')
hardcoded_seed <- 7670

## Dataset ####
load('fake_data.Rdata')
load('imputed_values.Rdata')

df.etm <- within(fake_data.etm, {
  ti <- exit
  ageCancer <- as.numeric(as.character(ageCancer))
  event <- ifelse(to == "cens", 0, ifelse(to == "cardio", 2, 1))
  id2 <- id
  TypeSMN <- factor(sample(c("Important", "Other"), replace = TRUE, size = length(exit)))
})
df.etm_birth <- within(df.etm, {
  ti <- ti + ageCancer; entry <- entry + ageCancer; exit <- exit + ageCancer
})
df.etm.death <- mutate(df.etm, event = ifelse(event == 1, 2, ifelse(event == 2, 1, event)))
df.etm_birth.death <- mutate(df.etm_birth, event = ifelse(event == 1, 2, ifelse(event == 2, 1, event)))

# Generate pseudo-values ####

grid.diag <- c(15, 20, 25, 30, 35)
names(grid.diag) <- paste(grid.diag, 'since diagnosis')
grid.birth <- c(20, 25, 30, 35)
names(grid.birth) <- paste(grid.birth, 'since birth')
grid.cutoff <- seq(2, 20, by = 2)

{
  # cd is the event of interest
  grid.pseudos_diag <- lapply(
    grid.diag, landmark_pseudos,
    dat = df.etm, covars = c(), cutoff = grid.cutoff,
    f_filter = landmark_filter_fusion_cancer
  )
  grid.pseudos_birth <- lapply(
    grid.birth, landmark_pseudos,
    dat = df.etm_birth, covars = c(), cutoff = grid.cutoff,
    f_filter = landmark_filter_fusion_cancer
  )
  df.pseudos <- c(grid.pseudos_diag, grid.pseudos_birth)
  for (ind in seq_along(df.pseudos)) df.pseudos[[ind]][["name"]] <- names(df.pseudos)[ind]
}
{
  # death is the event of interest
  grid.pseudos_diag.death <- lapply(
    grid.diag, landmark_pseudos,
    dat = df.etm.death, covars = c(), cutoff = grid.cutoff,
    f_filter = landmark_filter_fusion_cancer
  )
  grid.pseudos_birth.death <- lapply(
    grid.birth, landmark_pseudos,
    dat = df.etm_birth.death, covars = c(), cutoff = grid.cutoff,
    f_filter = landmark_filter_fusion_cancer
  )
  df.pseudos.death <- c(grid.pseudos_diag.death, grid.pseudos_birth.death)
  for (ind in seq_along(df.pseudos.death)) df.pseudos.death[[ind]][["name"]] <- names(df.pseudos.death)[ind]
}


# Configuration ####

# Variables used to replace NAs with imputed values
id_missing.sj_rt_heart <- unique(
  subset(df.etm, is.na(sj_rt_heart), 'id2', drop = TRUE))
id_missing.sj_rt_brain <- unique(
  subset(df.etm, is.na(sj_rt_brain), 'id2', drop = TRUE))
id_missing.sj_rt_neck <- unique(
  subset(df.etm, is.na(sj_rt_neck), 'id2', drop = TRUE))


# Define models (ie covariates used) ####

#' Insert SMN covariate in the list of covariates
#'
#' This is useful because I have two versions of the SMN covariate,
#' and want to use the same adjusted models on both.
add_smn_in_covs <- function(x, prefix = "", suffix = "", var = "Has2K"){
  names(x) <- paste0(prefix, names(x), suffix)
  for (ind in seq_along(x)){
    x[[ind]] <- c(var, x[[ind]])
  }
  x
}

## When some covariates are imputed ####

# covars_config_raw <- list(
#   "Model adjusted for cumulative doses for RT (Gy) and CT (yes/no)" = c(
#     "CT", "sj_Gender", "sj_AgeDiag", "sj_rt_heart", "cutYearDiag"),
#   "Model adjusted for cumulative doses for RT (Gy) and CT (mg/m2)" = c(
#     "sj_Gender", "sj_AgeDiag", "sj_anthra", "sj_alkyl_noplatinum", "sj_platinum",
#     "sj_rt_brain", "sj_rt_neck", "sj_rt_heart", "cutYearDiag")
# )
covars_config_raw <- list("Some model" = c("CT", "cutYearDiag"))

# create lists storing covariates combinations, when using imputed data
covars_config_base <- add_smn_in_covs(covars_config_raw, var = "Has2K")
for (ind in seq_along(covars_config_base)){
  attr(covars_config_base[[ind]], "name") <- names(covars_config_base)[ind]
}
covars_config_split <- add_smn_in_covs(covars_config_raw, "type-only ", var = "f_type_2k")
for (ind in seq_along(covars_config_split)){
  attr(covars_config_split[[ind]], "name") <- names(covars_config_split)[ind]
}

## When no covariates are imputed ####

# covars_config_raw.null <- list(
#   "Univariable" = c(),
#   "Model adjusted for RT (yes/no) and CT (yes/no)" =
#     c('sj_Gender', 'RT', 'CT', "sj_AgeDiag", "cutYearDiag"),
#   "Model adjusted for RT (yes/no) and cumulative doses of CT (mg/m2)" = c(
#     'sj_Gender', 'RT', 'sj_AgeDiag', 'sj_anthra', 'sj_alkyl_noplatinum', 'sj_platinum', 'cutYearDiag')
# )
covars_config_raw.null <- list("Some model" = c("RT", "CT", "cutYearDiag"))
covars_config_base.null <- add_smn_in_covs(covars_config_raw.null, var = "Has2K")
for (ind in seq_along(covars_config_base.null)){
  attr(covars_config_base.null[[ind]], "name") <- names(covars_config_base.null)[ind]
}
covars_config_split.null <- add_smn_in_covs(covars_config_raw.null, "type-only ", var = "f_type_2k")
for (ind in seq_along(covars_config_split.null)){
  attr(covars_config_split.null[[ind]], "name") <- names(covars_config_split.null)[ind]
}



# Functions ####

#' Converts a mira object to a summarized version of it
#'
#' Because I am storing model results after each fit,
#' it is useful to first reduce its size to save disk space.
mira_geeglm2shortobj <- function(x, ...){
  x$analyses <- lapply(
    x$analyses,
    function(obj){
      res <- list("coef" = coef(obj),
                  "vcov" = vcov(obj),
                  "summary" = summary(obj),
                  "tidy" = tidy(obj),
                  "glance" = glance(obj),
                  "residuals" = residuals(obj)
                  )
      class(res) <- c("shortobj", "list")
      res
    }
  )
  x
}
coef.shortobj <- function(x,...)x$coef
vcov.shortobj <- function(x,...)x$vcov
summary.shortobj <- function(x,...)x$summary
tidy.shortobj <- function(x,...)x$tidy
glance.shortobj <- function(x,...)x$glance
residuals.shortobj <- function(x,...)x$residuals

#' Fit the model.
#' 
#' @param covar_config Covariates to include in the model.
#' @param imputed_values either imputed data with mice or NULL.
estimate_one_config_results <- function(covar_config, df.pseudo,
                                        imputed_values,
                                        return_fit = TRUE, verbose = TRUE){
  id_run <- round(runif(1), 5) * 1e5
  if (verbose)
    cat("Starting estimate_one_config_results() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  out <- if (is.null(imputed_values)) {
    estimate(df.pseudo, covar_config, return_fit = TRUE)
  } else with(
    imputed_values, {
      # Replace NAs from df.pseudo
      for (idx in id_missing.sj_rt_heart){
        new_value <- sj_rt_heart[id2 == idx & !is.na(sj_rt_heart)][1]
        df.pseudo[[1]][df.pseudo[[1]]$id2 == idx, 'sj_rt_heart'] <- new_value
      }
      for (idx in id_missing.sj_rt_brain){
        new_value <- sj_rt_brain[id2 == idx & !is.na(sj_rt_brain)][1]
        df.pseudo[[1]][df.pseudo[[1]]$id2 == idx, 'sj_rt_brain'] <- new_value
      }
      for (idx in id_missing.sj_rt_neck){
        new_value <- sj_rt_neck[id2 == idx & !is.na(sj_rt_neck)][1]
        df.pseudo[[1]][df.pseudo[[1]]$id2 == idx, 'sj_rt_neck'] <- new_value
      }
      # GEE
      estimate(df.pseudo, covar_config, return_fit = TRUE)
    }
  )
  
  if (verbose)
    cat("Finished estimate_one_config_results() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  out
}
#' Returns a gee object fitted using covar_config and df.pseudo.
#' 
#' @param covar_config A vector containing the names of the covariates to use.
#' @param df.pseudo A list returned by \link{calc_pseudo}.
get_one_config_results <- function(covar_config, df.pseudo, imputed_values,
                                   RdataFolder = NULL, config_name = NULL,
                                   verbose = TRUE, force_recompute = TRUE){
  id_run <- round(runif(1), 5) * 1e5
  if (verbose)
    cat("Starting get_one_config_results() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  if (!is.null(config_name))
    config_name <- paste(config_name, attr(covar_config, "name"), sep = "-")
  target_rdata <- if (is.null(RdataFolder) || is.null(config_name)) {
    NULL
  } else {
    file.path(RdataFolder, paste0(fs::path_sanitize(config_name, "-"), '.Rdata'))
  }
  
  # Iterate over each covariate configuration
  if (force_recompute || is.null(target_rdata) || !file.exists(target_rdata)){
    one_config_results <- estimate_one_config_results(covar_config, df.pseudo,
                                                      imputed_values, verbose = verbose)
    pooled_result <- if (is.null(imputed_values)) one_config_results else pool(one_config_results)
    one_config_results <- mira_geeglm2shortobj(one_config_results)
    if (!is.null(target_rdata)) {
      cat("Saving ", target_rdata, "\n", sep = "",
          file = "debug-ci-mice.log", append = TRUE)
      save(one_config_results, file = target_rdata)
    }
  } else { load(target_rdata) }
  
  if (verbose)
    cat("Finished get_one_config_results() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  
  pooled_result
}

#' Returns a list of gee objects generated by fitting each element of
#' covar_configs to df.pseudo.
#' 
#' @param df.pseudo A list returned by \link{calc_pseudo}.
#' @param covar_configs A list of vectors containing the names of variables
#' to include in the model. tpseudo should not be included.
get_one_landmark_results <- function(df.pseudo, covar_configs, imputed_values,
                                     RdataFolder = NULL, config_name = NULL,
                                     verbose = TRUE, force_recompute = TRUE){
  id_run <- round(runif(1), 5) * 1e5
  if (verbose)
    cat("Starting get_one_landmark_result() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  if (is.null(config_name))
    config_name <- df.pseudo$name
  
  one_landmark_result <- lapply(
    covar_configs,
    get_one_config_results,
    df.pseudo = df.pseudo,
    imputed_values = imputed_values,
    RdataFolder = RdataFolder,
    config_name = config_name,
    force_recompute = force_recompute
  )
  
  if (verbose)
    cat("Finished get_one_landmark_result() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  
  one_landmark_result
}

#' Returns a list of list of gee object.
#' 
#' This function uses parallel computing. It returns a list, where each element
#' correpsonds to one landmark time (ie. one element of df.pseudos). Each
#' element of this list is a list, where each cell corresponds to one set of
#' covariates (ie. one element of covar_configs). This sub-list elements are
#' gee objects returned by the geepack::geeglm function.
#' 
#' @param df.pseudos A list of results of calc_pseudo. You should name it,
#' because the order of the result is not guaranted to be the same.
#' @param covar_configs A list of vectors containing the names of variables
#' to include in the model. tpseudo should not be included.
#' @param imputed_values A mids object.
get_all_landmark_results <- function(df.pseudos, covar_configs, imputed_values,
                                     RdataFolder = NULL,
                                     force_recompute = TRUE){
  cat("Fresh get_all_landmark_result()\n",
      file = "debug-ci-mice.log", append = TRUE)
  if (!is.null(RdataFolder)) dir.create(RdataFolder, recursive = TRUE, showWarnings = FALSE)
  
  # Iterate over each grid of pseudo-values (ie. each landmark time)
  
  cores_2_use <- min(detectCores() - 1, length(df.pseudos))
  
  cl <- makeCluster(cores_2_use)
  clusterSetRNGStream(cl, 9956)
  
  clusterExport(cl, c("id_missing.sj_rt_heart",
                      "id_missing.sj_rt_brain",
                      "id_missing.sj_rt_neck",
                      "get_one_config_results",
                      "estimate_one_config_results",
                      "estimate",
                      "mira_geeglm2shortobj"))
  clusterEvalQ(cl, {library('geepack'); library('mice')})
  
  
  all_landmark_results <- parLapply(
    cl = cl,
    X = df.pseudos,
    fun = get_one_landmark_results,
    covar_configs = covar_configs,
    imputed_values = imputed_values,
    RdataFolder = RdataFolder,
    force_recompute = force_recompute
  )
  stopCluster(cl)
  all_landmark_results
}



# Execution ####

covars_config_split <- covars_config_split[
  startsWith(names(covars_config_split), 'type-only')]
covars_config_split.null <- covars_config_split.null[
  startsWith(names(covars_config_split.null), 'type-only')]

## SMN on CD. ####
#' Estimate the additive effect on the cumulative incidence of SMN on CD.
#' Use a fixed window (pseudo-values) and include platinum agent as its own
#'   category (instead of including in as alkylating agent.)

# compute and store a version where SMN is split into multiple categories
{
  simu_folder.cd_split <- 'fake-ci/cd/split/'
  all_landmark_results_det <- get_all_landmark_results(
    df.pseudos = df.pseudos,
    covar_configs = covars_config_split,
    imputed_values = imputed_values,
    RdataFolder = simu_folder.cd_split,
    force_recompute = FALSE
  )
  save(
    all_landmark_results_det,
    file = file.path(simu_folder.cd_split, "all_landmark_results.Rdata")
  )
  all_landmark_results_det.null <- get_all_landmark_results(
    df.pseudos = df.pseudos,
    covar_configs = covars_config_split.null,
    imputed_values = NULL,
    RdataFolder = simu_folder.cd_split,
    force_recompute = FALSE
  )
  save(
    all_landmark_results_det.null,
    file = file.path(simu_folder.cd_split, "all_landmark_results.null.Rdata")
  )
}

# compute and store a version where SMN is a binary status.
{
  simu_folder.cd_base <- 'fake-ci/cd/base/'
  all_landmark_results_det <- get_all_landmark_results(
    df.pseudos = df.pseudos,
    covar_configs = covars_config_base,
    imputed_values = imputed_values,
    RdataFolder = simu_folder.cd_base,
    force_recompute = FALSE
  )
  save(all_landmark_results_det, file = file.path(simu_folder.cd_base, "all_landmark_results.Rdata"))
  
  all_landmark_results_det.null <- get_all_landmark_results(
    df.pseudos = df.pseudos,
    covar_configs = covars_config_base.null,
    imputed_values = NULL,
    RdataFolder = simu_folder.cd_base,
    force_recompute = FALSE
  )
  save(all_landmark_results_det.null, file = file.path(simu_folder.cd_base, "all_landmark_results.null.Rdata"))
}


## SMN on DEATH ####
#' Estimate additive effect of SMN on death.
#'  Use a fixed window (pseudo-values) and include platinum agent as its own
#'  category (instead of including in as alkylating agent.)

# compute and store a version where SMN is split into multiple categories
{
  simu_folder.death_split <- 'fake-ci/death/split/'
  all_landmark_results_det <- get_all_landmark_results(
    df.pseudos = df.pseudos.death,
    covar_configs = covars_config_split,
    imputed_values = imputed_values,
    RdataFolder = simu_folder.death_split,
    force_recompute = FALSE
  )
  save(all_landmark_results_det, file = file.path(simu_folder.death_split, "all_landmark_results.Rdata"))
  
  all_landmark_results_det.null <- get_all_landmark_results(
    df.pseudos = df.pseudos.death,
    covar_configs = covars_config_split.null,
    imputed_values = NULL,
    RdataFolder = simu_folder.death_split,
    force_recompute = FALSE
  )
  save(all_landmark_results_det.null, file = file.path(simu_folder.death_split, "all_landmark_results.null.Rdata"))
}

# compute and store a version where SMN is a binary status.
{
  simu_folder.death_base <- 'fake-ci/death/base/'
  all_landmark_results_det <- get_all_landmark_results(
    df.pseudos = df.pseudos.death,
    covar_configs = covars_config_base,
    imputed_values = imputed_values,
    RdataFolder = simu_folder.death_base,
    force_recompute = FALSE
  )
  save(all_landmark_results_det, file = file.path(simu_folder.death_base, "all_landmark_results.Rdata"))
  
  all_landmark_results_det.null <- get_all_landmark_results(
    df.pseudos = df.pseudos.death,
    covar_configs = covars_config_base.null,
    imputed_values = NULL,
    RdataFolder = simu_folder.death_base,
    force_recompute = FALSE
  )
  save(all_landmark_results_det.null, file = file.path(simu_folder.death_base, "all_landmark_results.null.Rdata"))
}


