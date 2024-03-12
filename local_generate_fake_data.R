library('dplyr')



### completly random data

n = 1e4
fake_ass_covs <- data.frame(
  CT = factor(
    sample(c("Yes", "No"), n, TRUE, c(.6, .4)),
    levels = c("Yes", "No")
  ),
  RT = factor(
    sample(c("Yes", "No"), n, TRUE, c(.8, .4)),
    levels = c("Yes", "No")
  ),
  iccc_group_lab = factor(
    sample(c("First", "Second", "Third", "Fourth"), n, TRUE, c(.1, .4, .3, .2)),
    levels = c("First", "Second", "Third", "Fourth")
  ),
  
  sj_Gender = factor(
    sample(c("Purple", "Wyvern"), n, TRUE),
    levels = c("Purple", "Wyvern")
  ),
  sj_AgeDiag = factor(
    sample(c("<1 y", "[1, 5[", "[5, 10[", "[10, 15[", ">=15 y"), n, TRUE, c(.3, .15, .2, .2, .15)),
    levels = c("<1 y", "[1, 5[", "[5, 10[", "[10, 15[", ">=15 y")
  ),
  cutYearDiag = factor(
    sample(c("Long ago", "Not so long ago", "Quite recently", "Now"), n, TRUE, c(.1, .3, .5, .1)),
    levels = c("Long ago", "Not so long ago", "Quite recently", "Now")
  ),
  
  sj_rt_heart = "0", sj_rt_brain = "0", sj_rt_neck = "No",
  sj_anthra = "0", sj_alkyl_noplatinum = "No", sj_platinum = "No"
)

fake_ass_covs2 <- within(fake_ass_covs, {
  ageCancer <- c(0, 1, 5, 10, 15)[as.integer(sj_AgeDiag)] +
    runif(length(sj_AgeDiag), 0, 1) * c(1, 4, 5, 5, 3)[as.integer(sj_AgeDiag)]
  
  n_RT <- sum(RT == "Yes")
  sj_rt_heart[RT == "Yes"] <- sample(c("]0, 5]", "]5, 15]", "]15, 35]", ">=35 "), n_RT, TRUE, c(.5, .35, .1, .05))
  sj_rt_heart <- factor(sj_rt_heart, levels = c("0", "]0, 5]", "]5, 15]", "]15, 35]", ">=35 "))
  sj_rt_brain[RT == "Yes"] <- sample(c("]0, 5]", "]5, 15]", "]15, 35]", ">=35 "), n_RT, TRUE, c(.6, .25, .05, .1))
  sj_rt_brain <- factor(sj_rt_brain, levels = c("0", "]0, 5]", "]5, 15]", "]15, 35]", ">=35 "))
  sj_rt_neck[RT == "Yes"] <- sample(c("Yes", "No"), n_RT, TRUE, c(.85, .15))
  sj_rt_neck <- factor(sj_rt_neck, levels = c("Yes", "No"))
  
  n_CT = sum(CT == "Yes")
  sj_anthra[CT == "Yes"] <- sample(c("0", "]0, 100]", "]100, 250]", ">250"), n_CT, TRUE, c(.3, .3, .25, .15))
  sj_anthra <- factor(sj_anthra, levels = c("0", "]0, 100]", "]100, 250]", ">250"))
  sj_alkyl_noplatinum[CT == "Yes"] <- sample(c("Yes", "No"), n_CT, TRUE, c(.7, .3))
  sj_alkyl_noplatinum <- factor(sj_alkyl_noplatinum, levels = c("Yes", "No"))
  sj_platinum[CT == "Yes"] <- sample(c("Yes", "No"), n_CT, TRUE, c(.4, .6))
  sj_platinum <- factor(sj_platinum, levels = c("Yes", "No"))
})

# fake events

fake_model <- model.matrix(
  ~ sj_Gender + sj_AgeDiag + cutYearDiag + sj_rt_heart + sj_rt_brain +
    sj_rt_neck + sj_anthra + sj_alkyl_noplatinum + sj_platinum,
  data = fake_ass_covs2)

beta_cancer <- matrix(nrow = ncol(fake_model), ncol = 1, c(
  .23, .001, -.05, -.018, -.024, -.036, # up to age >= 15y
  -.025, -.047, -.084, .02, .029, .08, .124, # up to heart >= 35
  .01, .0478, .041, .055, -.03, # up to rt_neck
  .031, .047, .142, -.023, .0002
))
prob_cancer <- fake_model %*% beta_cancer %>% as.vector()
prob_cancer[prob_cancer < 0] <- 0

beta_cardio <- matrix(nrow = ncol(fake_model), ncol = 1, c(
  .208, .003, -.037, -.027, -.029, -.036, # up to age >= 15y
  -.0145, -.0347, -.074, .035, .039, .0575, .152, # up to heart >= 35
  .01474, .0363, .051, .064, -.05, # up to rt_neck
  .02759, .0347, .112, -.0413, .0002
))
prob_cardio <- fake_model %*% beta_cardio %>% as.vector()
prob_cardio[prob_cardio < 0] <- 0

beta_death <- matrix(nrow = ncol(fake_model), ncol = 1, c(
  .5, -.04, -.037, -.057, -.089, -.146, # up to age >= 15y
  -0.207, -.2347, -.274, .035, .039, .0575, .152, # up to heart >= 35
  .01474, .0363, .051, .094, -.05, # up to rt_neck
  0.148, 0.247, 0.29215, -.0413, .0102
))
prob_death <- fake_model %*% beta_death %>% as.vector()
prob_death[prob_death < 0] <- 0

has_cancer <- rep(NA, n)
for (ind in seq_len(n))
  has_cancer[ind] <- sample(c(TRUE, FALSE), 1, TRUE, c(prob_cancer[ind], 1 - prob_cancer[ind]))
has_cardio <- rep(NA, n)
for (ind in seq_len(n))
  has_cardio[ind] <- sample(c(TRUE, FALSE), 1, TRUE, c(prob_cardio[ind], 1 - prob_cardio[ind]))
has_died <- rep(NA, n)
for (ind in seq_len(n))
  has_died[ind] <- sample(c(TRUE, FALSE), 1, TRUE, c(prob_death[ind], 1 - prob_death[ind]))


t_cancer <- runif(n, 10, 40)
t_cancer[!has_cancer] <- 100
t_cardio <- runif(n, 18, 60)
t_cardio[!has_cardio] <- 100

t_death <- runif(n, 15, 80)
t_death[has_cancer & t_death < t_cancer] <- t_cancer[has_cancer & t_death < t_cancer] +
  runif(sum(has_cancer & t_death < t_cancer), 0, 20)
t_death[has_cardio & t_death < t_cardio] <- 100

fake_ass_times <- within(fake_ass_covs2, {
  t_cardio <- t_cardio
  t_cancer <- t_cancer
  t_death <- t_death
  cens_cardio = has_cardio
  cens_cancer = has_cancer
  cens_death = has_died
})

fake_ass_etm <- phdCharrier::etmprep_fixtra(
  c(NA, "t_cardio", "t_cancer", "t_death"), c(NA, "cens_cardio", "cens_cancer", "cens_death"),
  data = fake_ass_times,
  tra = matrix(c(F,T,T,T, F,F,F,F, F,T,F,T, F,F,F,F), ncol=4, byrow=TRUE),
  state.names = c("init", "cardio", "cancer", "death"),
  cens.name = "cens",
  keep = colnames(fake_ass_times)
)

fake_ass_etm %>% group_by(from, to) %>% count()

fake_data.etm <- fake_ass_etm
save(fake_data.etm, file = "fake_data.Rdata")

