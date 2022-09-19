# title: "Effect of 2k on Cardiovascular Risk, A Cause-Specific-Hazard Analysis"
# author: "Thibaud Charrier"


# Description

# This .R file generates a .Rdata file and .csv files containing results
# for the following analysis of the 'CD increase after SMN' paper:
# Effect of SMN on CSH of CDs using SMN as a time-dependent binary covariate.
# For details about the methods, read the article.


# Loading ...

library('survival')
library('mice')
library('dplyr')
library('ggplot2')
library('kableExtra')
load('_data/df-etm.Rdata')

folder <- 'Results/csh-filled-mice/'

#' For the sake of beauty, converts booleans into 'Yes'/'No' factors.
bool2yes <- function(x, yes = 'Yes', no = 'No', inverse.levels = FALSE){
  b.levels <- if (inverse.levels) c(no, yes) else c(yes, no)
  y <- factor(rep(NA, length(x)), levels = b.levels)
  y[x] <- yes
  y[!x] <- no
  y
}

df.cshpreg <- within(
  df.etm,
  {
    # Build event variable
    to2 <- as.character(to)
    to2[to == 'cancer+deces'] <- 'deces'
    to2[to == 'cancer+cardio'] <- 'cardio'
    to2[to == 'cancer'] <- 'cens'
    to2 <- factor(to2, c('cens', 'cardio', 'deces'))
    event <- to2
    
    
}) %>%
  filter(exit > 5) %>%
  mutate(entry = ifelse(entry < 5, 5, entry))


# Multiple Imputation ##### Multiple Imputation ####

#' We want to impute missing values once per patient.
#' This data.frame allows just that.
mice_where <- with(df.cshpreg, {
  mice_where <- matrix(0, nrow = nrow(df.cshpreg), ncol = ncol(df.cshpreg))
  colnames(mice_where) <- colnames(df.cshpreg)
  first_occ <- df.cshpreg %>% select(id2) %>% group_by(id2) %>%
    mutate(first_occ = seq(1,n()) == 1) %>% pull(first_occ)
  mice_where[is.na(sj_rt_heart) & first_occ, 'sj_rt_heart'] <- 1
  mice_where[is.na(sj_rt_brain) & first_occ, 'sj_rt_brain'] <- 1
  mice_where[is.na(sj_rt_neck) & first_occ, 'sj_rt_neck'] <- 1
  
  mice_where
})
colnames(mice_where) <- colnames(df.cshpreg)

# generate mice results. (ie. multiple imputation)
imputed_values <- mice::mice(
  data = df.cshpreg,
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

# Legacy code
# imputed_values$data <- within(
#   imputed_values$data,
#   {
#     to2 <- as.character(to)
#     to2[to == 'cancer+deces'] <- 'deces'
#     to2[to == 'cancer+cardio'] <- 'cardio'
#     to2[to == 'cancer'] <- 'cens'
#     to2 <- factor(to2, c('cens', 'cardio', 'deces'))
#     event <- to2
#   }
# )
# imputed_values$data <- left_join(
#   imputed_values$data,
#   group_by(select(df.cshpreg, ctr, numcent, sj_alkyl_noplatinum, sj_platinum),
#            ctr, numcent) %>% filter(1:n() == 1)
# )

#' We imputed missing values once per patient. Duplicate those values where needed.
multistate_fix_imp <- function(missing_vec, ids){
  loc_missing <- is.na(missing_vec)
  for (ind in which(loc_missing)){
    id <- ids[ind]
    value <- missing_vec[ids == id & !is.na(missing_vec)]
    value <- if (length(value) > 1 & mean(as.integer(value)) != as.integer(value[1])){
      warning("Found different values, defaulting to the first")
      value[1]
    } else  if (length(value) == 0){
      warning("Found no value, defaulting to NA")
      NA
    } else value
    missing_vec[ind] <- value
  }
  missing_vec
}


# Models (w/o Imputation )####

## cg: univariable ####
cshpreg.uni <- coxph(
  Surv(entry, exit, event) ~ Has2K,
  data = df.cshpreg, id = id2
)

## cg: Very Little Information ####
cshpreg.vl <- coxph(
  list(Surv(entry, exit, event) ~ Has2K, 1:2 ~ RT + CT + sj_Gender + sj_AgeDiag),
  data = df.cshpreg, id = id2
)

## cg: Very Little Information ####
cshpreg.vl_yd <- coxph(
  list(Surv(entry, exit, event) ~ Has2K, 1:2 ~ RT + CT + sj_Gender +
         sj_AgeDiag + cutYearDiag),
  data = df.cshpreg, id = id2
)

## cg: Chimio Doses ####
cshpreg.ct <- coxph(
  list(Surv(entry, exit, event) ~ Has2K,
       1:2 ~ RT + sj_Gender + sj_AgeDiag + sj_anthra +
         sj_alkyl_noplatinum + sj_platinum),
  data = df.cshpreg, id = id2
)


# Models (w/ Imputation) ####

## cg: RT doses @ Heart ####

cshpreg_pooled.rt <- with(
  imputed_values,
  {
    sj_rt_heart <- multistate_fix_imp(sj_rt_heart, id2)
    (cshpreg.rt <- coxph(
      list(Surv(entry, exit, event) ~ Has2K,
           1:2 ~ CT + sj_Gender + sj_AgeDiag + sj_rt_heart
      ),
      id = id2
    ))
  }
) %>% pool()

## cg: Rodrigue ####
cshpreg_pooled.ra <- with(
  imputed_values,
  {
    sj_rt_heart <- multistate_fix_imp(sj_rt_heart, id2)
    coxph(
      list(
        Surv(entry, exit, event) ~ Has2K,
        1:2 ~ sj_Gender + sj_AgeDiag + sj_anthra +
          sj_alkyl_noplatinum + sj_platinum + sj_rt_heart
      ),
      id = id2
    )
  }
) %>% pool()

## cg: Rodrigue ####
cshpreg_pooled.ra_yd <- with(
  imputed_values,
  {
    sj_rt_heart <- multistate_fix_imp(sj_rt_heart, id2)
    coxph(
      list(Surv(entry, exit, event) ~ Has2K,
           1:2 ~ sj_Gender + sj_AgeDiag + sj_anthra +
             sj_alkyl_noplatinum + sj_platinum + sj_rt_heart +
            cutYearDiag),
      id = id2
    )
  }
) %>% pool()

## cg: St.Jude Extension ####
cshpreg_pooled.sj <- with(
  imputed_values,
  {
    sj_rt_heart <- multistate_fix_imp(sj_rt_heart, id2)
    sj_rt_brain <- multistate_fix_imp(sj_rt_brain, id2)
    sj_rt_neck <- multistate_fix_imp(sj_rt_neck, id2)
    coxph(
      list(
        Surv(entry, exit, event) ~ Has2K,
        1:2 ~ Has2K + sj_Gender + sj_AgeDiag +  sj_anthra +
          sj_alkyl_noplatinum + sj_platinum +
          sj_rt_brain + sj_rt_heart + sj_rt_neck
      ),
      id = id2
    )
  }
) %>% pool()

## cg: St.Jude Extension ####
cshpreg_pooled.sj_yd <- with(
  imputed_values,
  {
    coxph(
      list(
        Surv(entry, exit, event) ~ Has2K,
        1:2 ~ Has2K + sj_Gender + sj_AgeDiag +  sj_anthra +
          sj_alkyl_noplatinum + sj_platinum +
          sj_rt_brain + sj_rt_heart + sj_rt_neck + cutYearDiag
      ),
      data = df.cshpreg,
      id = id2
    )
  }
) %>% pool()

# Summarize Results & Save to files ####


## SMN ####
summary.csh.covariates <- c(
  "Univariable",
  "Basic Information",
  "Basic Information 2",
  "CT doses",
  "RT heart dose",
  'RT CT doses',
  'RT CT doses 2',
  "All doses",
  'All doses 2'
)
summary.csh.list <- list(
  univariable = cshpreg.uni,
  basic_info = cshpreg.vl,
  basic_info_2 = cshpreg.vl_yd,
  ct_doses = cshpreg.ct,
  rt_doses = cshpreg_pooled.rt,
  rt_ct_doses = cshpreg_pooled.ra,
  rt_ct_doses_2 = cshpreg_pooled.ra_yd,
  st_jude = cshpreg_pooled.sj,
  st_jude_2 = cshpreg_pooled.sj_yd
)
names(summary.csh.list) <- summary.csh.covariates

#' Transform a fit object into a clean and short data.frame.
clean_df <- function(x, ...) UseMethod("clean_df")
clean_df.coxphms <- function(x, ...){
  res <- cbind(summary(x)$coefficients, confint(x, level = .95))
  res <- as.data.frame(res)
  res$term <- rownames(res)
  rownames(res) <- NULL
  transmute(
    res,
    term = term,
    estimate = exp(coef),
    std.err = `robust se`,
    p.value = `Pr(>|z|)`,
    conf2.5 = exp(`2.5 %`),
    conf97.5 = exp(`97.5 %`)
  )
}
clean_df.mipo <- function(x, ...){
  res <- cbind(
    summary(x),
    { y <- mice:::confint.mipo(x);
    data.frame(conf2.5 = y[, '2.5 %'], conf97.5 = y[, '97.5 %']) }
  )
  rownames(res) <- NULL
  transmute(
    res,
    term = term,
    estimate = exp(estimate),
    std.err = std.error,
    p.value = p.value,
    conf2.5 = exp(conf2.5),
    conf97.5 = exp(conf97.5)
  )
}

# cleanly summarize models results.
summary.csh.est <- do.call(rbind, lapply(summary.csh.list, clean_df))
summary.csh.est$covariates <- stringr::str_remove(rownames(summary.csh.est), '\\.[0-9]+')
rownames(summary.csh.est) <- NULL

summarized_smn_effect <- filter(summary.csh.est, term =='Has2KTRUE_1:2')

## Rename Treatment ####

#' Separate the variable-level info in two columns.
#' 
#' Warning: If you're using unsupported covariates names the results
#' will be misleading.
clear_trt_info <- function(x){
  covars <- c(
    "RT", "CT", "sj_Gender",
    "sj_AgeDiag", "sj_anthra", "sj_alkyl_noplatinum", "sj_platinum",
    "sj_rt_heart", "sj_rt_brain", "sj_rt_neck", "cutYearDiag"
  )
  covars_names <- c(
    'RT', 'CT', 'Sex',
    'Age at diagnosis', 'Anthracyclines doses',
    'Alkylating agent', 'Platinum agent',
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

summarized_treatment_infos <- clear_trt_info(
  filter(summary.csh.est, !stringr::str_detect(term, 'Has2KTRUE'))
)


