#' Impute missing radiotherapy dosimetry data
#' 
#' Those imputed data do not need to be run at each analysis,
#' and should be (and are) shared between csh and ci analysis.
#' 
#' Do not run again without deleting .Rdata (or changing its name) first,
#' otherwise a stop() will occur

load("fake_data.Rdata")
df.etm <- fake_data.etm

library('dplyr')
library('mice')

hardcoded_seed <- 7670
storage_file <- "imputed_values.Rdata"

# do not impute data multiple time for same patient.
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

# actually impute the data
imputed_values <- mice::mice(
  data = df.etm,
  m = 3,
  blocks = list(sj_rt_brain = 'sj_rt_brain', sj_rt_heart = 'sj_rt_heart',
                sj_rt_neck = 'sj_rt_neck'),
  formulas = list(
    sj_rt_brain = ~ cutYearDiag +
      CT * sj_AgeDiag + sj_anthra +
      sj_alkyl_noplatinum + sj_platinum + iccc_group_lab + sj_Gender +
      sj_rt_heart:(cutYearDiag + sj_AgeDiag) +
      sj_rt_neck:(cutYearDiag + sj_AgeDiag),
    sj_rt_heart = ~ cutYearDiag +
      CT * sj_AgeDiag + sj_anthra +
      sj_alkyl_noplatinum + sj_platinum  + iccc_group_lab + sj_Gender +
      sj_rt_brain:(cutYearDiag + sj_AgeDiag) +
      sj_rt_neck:(cutYearDiag + sj_AgeDiag)
    , sj_rt_neck = ~ cutYearDiag +
      CT * sj_AgeDiag + sj_anthra +
      sj_alkyl_noplatinum + sj_platinum  + iccc_group_lab + sj_Gender +
      sj_rt_brain:(cutYearDiag + sj_AgeDiag) +
      sj_rt_heart:(cutYearDiag + sj_AgeDiag)
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
  maxit = 3,
  # We have a quick convergence, this lowers the computation cost
  seed = hardcoded_seed
)

if (file.exists(storage_file)) {
  stop('storage_file already exists, do not overwrite it for safety.')
} else save(imputed_values, file = storage_file)
