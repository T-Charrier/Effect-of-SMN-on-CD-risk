# Compute cause-specific-hazard estimates


# Loading ...

library('survival')
library('mice')
library('dplyr')
library('ggplot2')
library('kableExtra')
source('landmark_functions.R')

load('fake_data.Rdata')
load('imputed_values.Rdata')

# Final data preparation step ####

df.cshpreg <- within(
  fake_data.etm,
  {
    # Build event variable
    to2 <- as.character(to)
    to2[to == 'cancer+deces'] <- 'deces'
    to2[to == 'cancer+cardio'] <- 'cardio'
    to2[to == 'cancer'] <- 'cens'
    to2 <- factor(to2, c('cens', 'cardio', 'deces'))
    event <- to2
    Has2K = (from == "cancer")
    f_type_2k = ifelse(Has2K, TypeSMN, "Hello")
  }) %>%
  filter(exit > 5) %>%
  mutate(entry = ifelse(entry < 5, 5, entry))

# copy modifications on imputed data.frame
imputed_values$data <- within(
  imputed_values$data,
  {
    to2 <- as.character(to)
    to2[to == 'cancer+deces'] <- 'deces'
    to2[to == 'cancer+cardio'] <- 'cardio'
    to2[to == 'cancer'] <- 'cens'
    to2 <- factor(to2, c('cens', 'cardio', 'deces'))
    event <- to2
    Has2K = (from == "cancer")
    f_type_2k = ifelse(Has2K, TypeSMN, "Hello")
  }
)
# imputed_values$data <- left_join(
#   imputed_values$data,
#   group_by(select(df.cshpreg, ctr, numcent, sj_alkyl_noplatinum, sj_platinum),
#            ctr, numcent) %>% filter(1:n() == 1)
# )

#' Copy imputed data to fill non-imputed index
#' 
#' @param missing_vec F/T vector saying which rows where imputed
#' @param ids of patients.
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



# Analysis by occurrence of SMN ####

# Models no using imputed data

## Univariable ####
csh.univariable <- coxph(
  Surv(entry, exit, event) ~ Has2K,
  data = df.cshpreg, id = id2
)
## Model adjusted for RT (yes/no) and CT (yes/no) ####
csh.bin_bin <- coxph(
  Surv(entry, exit, event) ~ Has2K + RT + CT + sj_Gender +
    sj_AgeDiag + cutYearDiag,
  data = df.cshpreg, id = id2
)

## Model adjusted for RT (yes/no) and cumulative doses of CT (mg/m2) ####
csh.bin_cumdo <- coxph(
  Surv(entry, exit, event) ~ Has2K + RT + sj_Gender + sj_AgeDiag +
    sj_anthra + sj_alkyl_noplatinum + sj_platinum + cutYearDiag,
  data = df.cshpreg, id = id2
)


# Models using imputed data

## Model adjusted for cumulative doses for RT (Gy) and CT (yes/no) ####
csh.cumdo_bin <- with(
  imputed_values,
  {
    sj_rt_heart <- multistate_fix_imp(sj_rt_heart, id2)
    (cshpreg.rt <- coxph(
      Surv(entry, exit, event) ~ Has2K + CT + sj_Gender + sj_AgeDiag + sj_rt_heart + cutYearDiag,
      id = id2
    ))
  }
) %>% pool()

## Model adjusted for cumulative doses for RT (Gy) and CT (mg/m2) ####
csh.cumdo_cumdo <- with(
  imputed_values,
  {
    coxph(
      Surv(entry, exit, event) ~ Has2K + sj_Gender + sj_AgeDiag +  sj_anthra +
        sj_alkyl_noplatinum + sj_platinum +
        sj_rt_brain + sj_rt_heart + sj_rt_neck + cutYearDiag,
      data = df.cshpreg,
      id = id2
    )
  }
) %>% pool()



# Analysis by type of SMN ####

# Models no using imputed data

## Univariable ####
csh.type.univariable <- coxph(
  Surv(entry, exit, event) ~ f_type_2k,
  data = df.cshpreg, id = id2
)
## Model adjusted for RT (yes/no) and CT (yes/no) ####
csh.type.bin_bin <- coxph(
  Surv(entry, exit, event) ~ f_type_2k + RT + CT + sj_Gender +
    sj_AgeDiag + cutYearDiag,
  data = df.cshpreg, id = id2
)

## Model adjusted for RT (yes/no) and cumulative doses of CT (mg/m2) ####
csh.type.bin_cumdo <- coxph(
  Surv(entry, exit, event) ~ f_type_2k + RT + sj_Gender + sj_AgeDiag +
    sj_anthra + sj_alkyl_noplatinum + sj_platinum + cutYearDiag,
  data = df.cshpreg, id = id2
)


# Models using imputed data

## Model adjusted for cumulative doses for RT (Gy) and CT (yes/no) ####
csh.type.cumdo_bin <- with(
  imputed_values,
  {
    sj_rt_heart <- multistate_fix_imp(sj_rt_heart, id2)
    (cshpreg.rt <- coxph(
      Surv(entry, exit, event) ~ f_type_2k + CT + sj_Gender + sj_AgeDiag + sj_rt_heart + cutYearDiag,
      id = id2
    ))
  }
) %>% pool()

## Model adjusted for cumulative doses for RT (Gy) and CT (mg/m2) ####
csh.type.cumdo_cumdo <- with(
  imputed_values,
  {
    coxph(
      Surv(entry, exit, event) ~ f_type_2k + sj_Gender + sj_AgeDiag +  sj_anthra +
        sj_alkyl_noplatinum + sj_platinum +
        sj_rt_brain + sj_rt_heart + sj_rt_neck + cutYearDiag,
      data = df.cshpreg,
      id = id2
    )
  }
) %>% pool()

