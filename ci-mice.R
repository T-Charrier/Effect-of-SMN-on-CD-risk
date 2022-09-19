# title: "Effect of 2k on Cardiovascular Risk, Cumulative Incidence Conditional
# on Survival at Time After First Diagnosis, a FCCSS Study."
# author: "Thibaud Charrier"


# Description

#' This .R file generates a .Rdata file and .csv files containing results
#' for the following analysis of the 'CD increase after SMN' paper:
#' Effect of SMN on cumulative incidence of CD, using years since diagnosis as a time-scale.  
#' Multiple Imputations is used.
#' For details about the methods, read the article.

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
library('tictoc')
source('landmark_functions.R')
hardcoded_seed <- 7670

## Dataset ####
load('_data/fake-data.Rdata')
df.etm_birth <- within(df.etm, {
  ti <- ti + ageCancer; entry <- entry + ageCancer; exit <- exit + ageCancer
})

# Configuration of models ####

# Variables used to replace NAs with imputed values
id_missing.sj_rt_heart <- unique(
  subset(df.etm, is.na(sj_rt_heart), 'id2', drop = TRUE))
id_missing.sj_rt_brain <- unique(
  subset(df.etm, is.na(sj_rt_brain), 'id2', drop = TRUE))
id_missing.sj_rt_neck <- unique(
  subset(df.etm, is.na(sj_rt_neck), 'id2', drop = TRUE))

# Set the covariates to use in the model. 
covar_configs <- list(
  "RT heart dose" = c("Has2K", "CT", "sj_Gender", "sj_AgeDiag", "sj_rt_heart"),
  "RT CT doses" = c("Has2K", "sj_Gender", "sj_AgeDiag", "sj_anthra",
                    "sj_alkyl_noplatinum", "sj_platinum", "sj_rt_heart"),
  "RT CT doses 2" = c("Has2K", "sj_Gender", "sj_AgeDiag", "sj_anthra",
                      "sj_alkyl_noplatinum", "sj_platinum",
                      "sj_rt_heart", "cutYearDiag"),
  "All doses" = c("Has2K", "sj_Gender", "sj_AgeDiag", "sj_anthra",
                  "sj_alkyl_noplatinum", "sj_platinum",
                  "sj_rt_brain", "sj_rt_neck", "sj_rt_heart"),
  "All doses 2" = c("Has2K", "sj_Gender", "sj_AgeDiag", "sj_anthra",
                    "sj_alkyl_noplatinum", "sj_platinum",
                    "sj_rt_brain", "sj_rt_neck", "sj_rt_heart", "cutYearDiag")
)
for (ind in seq_along(covar_configs)){
  attr(covar_configs[[ind]], "name") <- names(covar_configs)[ind]
}
covar_configs.null <- list(
  "Univariable" = c("Has2K"),
  "Basic Information" = c('Has2K', 'sj_Gender', 'RT', 'CT', "sj_AgeDiag"),
  "Basic Information 2" = c('Has2K', 'sj_Gender', 'RT', 'CT', 'sj_AgeDiag', 'cutYearDiag'),
  "CT doses" = c('Has2K', 'sj_Gender', 'RT', 'sj_AgeDiag', 'sj_anthra', 'sj_alkyl_noplatinum', 'sj_platinum')
)
for (ind in seq_along(covar_configs.null)){
  attr(covar_configs.null[[ind]], "name") <- names(covar_configs.null)[ind]
}

# Multiple Imputation ####

#' We want to impute missing values once per patient.
#' This data.frame allows just that.
mice_where <- with(df.etm, {
  mice_where <- matrix(0, nrow = nrow(df.etm), ncol = ncol(df.etm))
  colnames(mice_where) <- colnames(df.etm)
  first_occ <- df.etm %>% select(id2) %>% group_by(id2) %>%
    mutate(first_occ = seq(1,n()) == 1) %>% pull(first_occ)
  mice_where[is.na(sj_rt_heart) & first_occ, 'sj_rt_heart'] <- 1
  mice_where[is.na(sj_rt_brain) & first_occ, 'sj_rt_brain'] <- 1
  mice_where[is.na(sj_rt_neck) & first_occ, 'sj_rt_neck'] <- 1
  
  mice_where
})
colnames(mice_where) <- colnames(df.etm)

# generate mice results. (ie. multiple imputation)
imputed_values <- mice::mice(
  data = df.etm,
  m = 20,
  blocks = list(sj_rt_brain = 'sj_rt_brain', sj_rt_heart = 'sj_rt_heart',
                sj_rt_neck = 'sj_rt_neck'),
  formulas = list(
    sj_rt_brain = ~ cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) +
      CT * sj_AgeDiag + sj_anthra +
      sj_alkyl  + iccc_group_lab + sj_Gender +
      sj_rt_heart:(cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) + sj_AgeDiag) +
      sj_rt_neck:(cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) + sj_AgeDiag),
    sj_rt_heart = ~ cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) +
      CT * sj_AgeDiag + sj_anthra +
      sj_alkyl  + iccc_group_lab + sj_Gender +
      sj_rt_brain:(cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) + sj_AgeDiag) +
      sj_rt_neck:(cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) + sj_AgeDiag)
    , sj_rt_neck = ~ cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) +
                           CT * sj_AgeDiag + sj_anthra +
                           sj_alkyl  + iccc_group_lab + sj_Gender +
                           sj_rt_brain:(cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) + sj_AgeDiag) +
      sj_rt_heart:(cut(yearDiag, breaks = c(0, 1980, 1990, Inf)) + sj_AgeDiag)
  ),
  # This is temporary, and needs further work to make sure the imputations
  # are correct
  print = FALSE,
  # We do not want a verbose mice, we want an efficient one.
  ignore = with(df.etm, !(RT == "Yes" & from == 'init')),
  # from == init makes sure each patient has equal weight
  # RT == "Yes" makes sure we don't predict rt dose using patient unexposed
  # to radiotherapy
  where = mice_where,
  # mice_where makes sure each patient has only one imputed value.
  # This requires to fill the remaining NAs later in with.mice.
  maxit = 7,
  # We have a quick convergence, this lowers the computation cost
  seed = hardcoded_seed
)

# Generate pseudo-values ####
# choose the landmark times
grid.diag <- c(15, 20, 25, 30, 35)
names(grid.diag) <- paste(grid.diag, 'since diagnosis')
grid.birth <- c(20, 25, 30, 35)
names(grid.birth) <- paste(grid.birth, 'since birth')
# choose the pseudo-values grid points. They are years after landmark time.
grid.cutoff <- seq(2, 20, by = 2)
# use the time-scale "years since diagnosis"
grid.pseudos_diag <- lapply(
  grid.diag, landmark_pseudos,
  dat = df.etm, covars = c(), cutoff = grid.cutoff,
  f_filter = landmark_filter_fusion_cancer
)
# use the time-scale "attained age, years"
grid.pseudos_birth <- lapply(
  grid.birth, landmark_pseudos,
  dat = df.etm_birth, covars = c(), cutoff = grid.cutoff,
  f_filter = landmark_filter_fusion_cancer
)
df.pseudos <- c(grid.pseudos_diag, grid.pseudos_birth)
for (ind in seq_along(df.pseudos)){
  df.pseudos[[ind]][["name"]] <- names(df.pseudos)[ind]
}


# Estimation-related functions ####

#' Converts a mira object to a summarized version of it.
#' This is used to save disk space.
#' To disable it, use `mira_geeglm2shortobj <- identity`
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

#' 
estimate_one_config_results <- function(covar_config, df.pseudo,
                                        imputed_values,
                                        return_fit = TRUE){
  if (is.null(imputed_values)) {
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
}

#' Returns a gee object fitted using covar_config and df.pseudo.
#' 
#' @param covar_config A vector containing the names of the covariates to use.
#' @param df.pseudo A list returned by \link{calc_pseudo}.
#' @param verbose Boolean. Print to debug file.
get_one_config_results <- function(covar_config, df.pseudo, imputed_values,
                                   include_no_trt_config = FALSE,
                                   RdataFolder = NULL, config_name = NULL,
                                   verbose = TRUE){
  id_run <- round(runif(1), 5) * 1e5
  if (verbose)
    cat("Starting get_one_config_results() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  if (!is.null(config_name))
    config_name <- paste(config_name, attr(covar_config, "name"), sep = "-")
  
  # Iterate over each covariate configuration
  one_config_results <- estimate_one_config_results(covar_config,
                                                    df.pseudo,
                                                    imputed_values)
  pooled_result <- if (is.null(imputed_values)) one_config_results else pool(one_config_results)
  one_config_results <- mira_geeglm2shortobj(one_config_results)
  if (!is.null(RdataFolder) & !is.null(config_name))
    save(one_config_results, file = paste0(RdataFolder, config_name, '.Rdata'))
  
  # same thing but without trt
  if (include_no_trt_config &&
      length(covar_config.notrt <- covar_config[covar_config != 'Has2K']) > 0){
    one_config_results.notrt <- estimate_one_config_results(
      covar_config.notrt, df.pseudo, imputed_values)
    pooled_result.notrt <- if (is.null(imputed_values)) {
      one_config_results.notrt } else pool(one_config_results.notrt)
    one_config_results.notrt <- mira_geeglm2shortobj(one_config_results.notrt)
    if (!is.null(RdataFolder) & !is.null(config_name))
      save(one_config_results.notrt,
           file = paste0(RdataFolder, config_name, '-notrt.Rdata'))
  } else pooled_result.notrt <- NULL
  
  if (verbose)
    cat("Finished get_one_config_results() #", id_run, "\n", sep = "",
        file = "debug-ci-mice.log", append = TRUE)
  
  list(trt = pooled_result, notrt = pooled_result.notrt)
}

#' Returns a list of gee objects generated by fitting each element of
#' covar_configs to df.pseudo.
#' 
#' @param df.pseudo A list returned by \link{calc_pseudo}.
#' @param covar_configs A list of vectors containing the names of variables
#' to include in the model. tpseudo should not be included.
#' @param verbose Boolean. Print to debug file.
get_one_landmark_results <- function(df.pseudo, covar_configs, imputed_values,
                                     include_no_trt_config = FALSE,
                                     RdataFolder = NULL, config_name = NULL,
                                     verbose = TRUE){
  id_run <- round(runif(1), 5) * 1e5 # conflict isn't a strong issue.
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
    include_no_trt_config = include_no_trt_config,
    RdataFolder = RdataFolder,
    config_name = config_name
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
                                     include_no_trt_config = FALSE,
                                     RdataFolder = "MiceGeeglm/wip-model/"){
  cat("Fresh get_all_landmark_result()\n",
      file = "debug-ci-mice.log", append = TRUE)
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
    include_no_trt_config = include_no_trt_config,
    RdataFolder = RdataFolder
  )
  stopCluster(cl)
  all_landmark_results
}


# Execution ####

tmp_folder <- NULL # Set it to an existing folder to store results
                    # after they are computed. Can be useful for debug or
                    # to do the analysis in multiple times.
all_landmark_results_plat <- get_all_landmark_results(
  df.pseudos = df.pseudos, covar_configs = covar_configs,
  imputed_values = imputed_values, RdataFolder = tmp_folder,
  include_no_trt_config = TRUE
)
#' Remove comment of the next line to save results to disk.
#' Computation time can be high, so beware if you don't save them.
# save(all_landmark_results_plat, file = file.path(tmp_folder, "all_landmark_results.Rdata"))


all_landmark_results_plat.null <- get_all_landmark_results(
  df.pseudos = df.pseudos, covar_configs = covar_configs.null,
  imputed_values = NULL, RdataFolder = tmp_folder,
  include_no_trt_config = TRUE
)
#' Remove comment of the next line to save results to disk.
#' Computation time can be high, so beware if you don't save them.
# save(all_landmark_results_plat.null, file = file.path(tmp_folder, "all_landmark_results.null.Rdata"))


# Save Results ####
#' Save a short summarized version of the results.

clean_df <- function(x, ...) UseMethod("clean_df")
clean_df.default <- function(x, ...){
  data.frame(
    term = NA, estimate = NA, std.err = NA,
    p.value = NA, conf2.5 = NA, conf97.5 = NA
  )
}
clean_df.geeglm <- function(x, ...){
  res <- cbind(summary(x)$coefficients, confint(x, level = .95))
  res$term <- rownames(res)
  rownames(res) <- NULL
  transmute(
    res,
    term = term,
    estimate = Estimate,
    std.err = Std.err,
    p.value = `Pr(>|W|)`,
    conf2.5 = `2.5 %`,
    conf97.5 = `97.5 %`
  )
}
clean_df.mipo <- function(x, ...){
  res <- cbind(
    summary(x),
    { y <- mice:::confint.mipo(x);
    data.frame(conf2.5 = y[, '2.5 %'], conf97.5 = y[, '97.5 %']) },
    x$pooled[, c("m", "ubar", "b", "t", "dfcom", "riv", "lambda", "fmi")]
  )
  rownames(res) <- NULL
  transmute(
    res,
    term = term,
    estimate = estimate,
    std.err = std.error,
    p.value = p.value,
    conf2.5 = conf2.5,
    conf97.5 = conf97.5
  )
}

get_2k_info <- function(all_landmark_results){
  res <- lapply(all_landmark_results, function(one_landmark_result) {
    res <- do.call(rbind, lapply(one_landmark_result, function(x){
      obj <- clean_df(x$trt)
      select(filter(obj, term == 'Has2KTRUE'), -term)
    }))
    res$covariates <- rownames(res)
    rownames(res) <- NULL
    res
  })
  for (name in names(res)) {
    res[[name]]$time_scale <- ifelse(
      stringr::str_detect(name, 'since birth'),
      'Attained age', 'Time since diagnosis')
    res[[name]]$landmark_time <- as.integer(stringr::str_sub(name , 1, 2))
  }
  names(res) <- NULL
  res <- do.call(rbind, res)
  res
}

table_all_infos <- function(all_landmark_results, smn_version = TRUE){
  res <- lapply(all_landmark_results, function(one_landmark_result) {
    res <- do.call(rbind, lapply(one_landmark_result, function(x){
      xx <- if (smn_version) x$trt else x$notrt
      clean_df(xx)
    }))
    res$covariates <- stringr::str_remove(rownames(res), '\\.[0-9]+')
    rownames(res) <- NULL
    res
  })
  for (name in names(res)) {
    res[[name]]$time_scale <- ifelse(
      stringr::str_detect(name, 'since birth'),
      'Attained age', 'Time since diagnosis')
    res[[name]]$landmark_time <- as.integer(stringr::str_sub(name , 1, 2))
  }
  names(res) <- NULL
  filter(do.call(rbind, res), !is.na(term))
}

#' Separate the term-level part of the fit results in two columns
#' 
#' @param x a data.frame with column term.
clear_trt_info <- function(x){
  covars <- c(
    "RT", "CT", "sj_Gender",
    "sj_AgeDiag", "sj_anthra", "sj_alkyl_noplatinum", "sj_alkyl", "sj_platinum",
    "sj_rt_heart", "sj_rt_brain", "sj_rt_neck", "cutYearDiag"
  )
  covars_names <- c(
    'RT', 'CT', 'Sex',
    'Age at diagnosis', 'Anthracyclines doses',
    'Alkylating agent', 'Alkylating agent', 'Platinum agent',
    'Mean heart RT dose', 'Mean brain RT dose', 'Neck RT', 
    'Year of first cancer diagnosis'
  )
  x <- filter(x, !stringr::str_detect(term, 'Has2KTRUE'))
  covs <- rep(NA, nrow(x))
  labs <- rep(NA, nrow(x))
  for (x_ind in seq_along(x$term)){
    for (cov_ind in seq_along(covars)){
      if (stringr::str_detect(x$term[x_ind], covars[cov_ind])){
        covs[x_ind] <- covars_names[cov_ind]
        labs[x_ind] <- stringr::str_sub(x$term[x_ind],
                                        1 + stringr::str_length(covars[cov_ind]))
        break
      }
    }
  }
  x$term <- covs
  x$labs <- labs
  x
}

#' Summarize and save to .csv results of one or more results of get_all_landmark_results
#' 
#' Create .csv files containing summarized information of all the ... . Names
#' of ... argument are ignored. They are all written to the same file. Identical
#' configuration will not issue any error or warning, but you will not be
#' able to distinguish them - expect using the writing order, which is
#' preserved. Created files are 'Estimates of SMN effect.csv',
#' 'Estimates of TRT effect.csv' and 'Estimates of TRT effect without SMN.csv'.
#' The time effect (ie. baseline cumulative incidence function) are not stored.
#' 
#' @param folder folder to create .csv files in.
#' @param ... results of get_all_landamrk_results. Names are ignored.
save_infos <- function(folder, ...){
  table_2k_estimates <- do.call(rbind, lapply(list(...), get_2k_info))
  table_all_landmark_results <- do.call(rbind, lapply(list(...), table_all_infos)) %>%
    filter(!stringr::str_detect(term, 'Has2K'),
           !stringr::str_detect(term, 'as\\.factor\\(tpseudo\\)'))
  table_all_landmark_results_no_smn <- do.call(
    rbind, lapply(list(...), table_all_infos, smn_version = FALSE)) %>%
    filter(!stringr::str_detect(term, 'Has2K'),
           !stringr::str_detect(term, 'as\\.factor\\(tpseudo\\)'))
  
  write.csv(table_2k_estimates,
            file = file.path(folder, 'Estimates of SMN effect.csv'),
            row.names = FALSE)
  write.csv(clear_trt_info(table_all_landmark_results),
            file = file.path(folder, 'Estimates of TRT effect.csv'),
            row.names = FALSE)
  write.csv(clear_trt_info(table_all_landmark_results_no_smn),
            file = file.path(folder, 'Estimates of TRT effect without SMN.csv'),
            row.names = FALSE)
}

save_infos('Results/ci-mice-plat',
           all_landmark_results_plat, all_landmark_results_plat.null)
