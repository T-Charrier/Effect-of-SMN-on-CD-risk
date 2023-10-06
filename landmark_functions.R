#' Compute pseudo-values with Aalen-Joanson
#' 
#' @param data Dataset to be analyzed, with event=1 for the cause of interest
#' @param analysis MUST be "other"
#' @param usepseudo logical, should we use the package "pseudo" to compute
#' pseudo values ? Both option yield very similar results (<10e-12).
calc_pseudo <- function(data, analysis="other", cutoffs=NULL, poi=seq(0.3,0.9,0.1),
                      poievt = c(1,2), usepseudo = FALSE){
  # Calculate pseudo-values using prodlim package
  data <- data[order(data$ti),]
  if (is.null(poievt)) poievt <- c(1,2)
  if (is.null(cutoffs))
    cutoffs <- quantile(data$ti[data$event %in% poievt], probs=poi)
  b <- NULL
  # usepseudo is slower, but requires less code.
  # Implementations return negligible difference (<10e-12)
  if (usepseudo) { # use pseudo::pseudoci
    rawpseudos <- with(data, pseudo::pseudoci(ti, event, cutoffs))
    rawpseudos1 <- rawpseudos$pseudo$cause1
    for (j in 1:ncol(rawpseudos1)){
      b <- rbind(b, cbind(data, pseudo=rawpseudos1[,j], tpseudo=cutoffs[j],
                          id=1:nrow(data)))
    }
  } else { # use prodlim::jacknife.competing.risks
    pl <- prodlim::prodlim(Hist(ti, event) ~ 1, data=data)
    jack <- prodlim::jackknife.competing.risks(pl, times=cutoffs, cause=1)
    
    # Create an augmented dataset with one line per patient per time-point in grid and
    # corresponding pseudo-value
    for (j in 1:ncol(jack)){
      b <- rbind(b, cbind(data, pseudo=jack[,j], tpseudo=cutoffs[j], id=1:nrow(data)))
    }
  }
  b <- b[order(b$id),]
  
  # Returns a list with augmented dataset, grid of time-points and number of patients
  return(list(b,cutoffs,nrow(data), n.evt = sum(data$event == 1),
              n.evt.2k = sum(data$event == 1 & data$Has2K),
              n.2k = sum(data$Has2K)))
}


#' Return a dataset filtered according to a landmark strategy at time s,
#' and shift to a 2-state output using cancer information.
landmark_filter_fusion_cancer <- function(s, dat, covars = c()){
  df.2risk.s <- dat %>%
    group_by(id2) %>%
    mutate(entrycancer = sum(entry * (from == 'cancer')),
           everHad2k = sum(from == 'cancer')) %>%
    ungroup() %>%
    filter(to != 'cancer', ti > s) %>%
    mutate(from = factor('init'))
  df.2risk.s <- within(
    df.2risk.s, {
      event[to == 'cens'] <- 0
      Has2K <- (s >= entrycancer) & everHad2k
      f_type_2k <- ifelse(Has2K, TypeSMN, "Hello")
      entry <- 0
      tmp_entrycancer <- NULL
  })
  for (covar in covars){
    df.2risk.s <- df.2risk.s[!is.na(df.2risk.s[, covar]), ]
  }
  df.2risk.s
}


#' Call calc_pseudo on a dataset filtered on individuals surviving up to s,
#' and whose covariates are known. (ie. not-NA)
#'
#' Please note other  covariates can be added to the model. This is just a nice
#' easy way to filter out NA's pre-fitting
#' 
#' @param s landmark time
#' @param dat data.frame with covars and ti, event.
#' @param covars Remove rows for which one or more of those columns have missing values.
#' @param cutoff Grid of pseudo-values. s will be added to it.
landmark_pseudos <- function(s, dat, covars = c(), cutoff=NULL, poi=NULL, poievt=NULL, f_filter = NULL){
  dat.s <- f_filter(s, dat, covars)
  if (!is.null(cutoff)) cutoff <- cutoff + s
  dat.s.pseudo <- calc_pseudo(dat.s, poi=poi, cutoffs=cutoff, poievt=poievt)
  dat.s.pseudo[[1]]$s <- s
  dat.s.pseudo
}

#' Estimate coefficients using an additive model.
#' 
#' @param list A list as in the output from calc_pseudo.
estimate <- function(list, covar, return_fit = FALSE){
  b <- list[[1]]
  cutoffs <- list[[2]]
  n <- list[[3]]
  # Fit model using geese function from geepack package
  form <- as.formula(paste("pseudo~as.factor(tpseudo)+", paste(covar,collapse="+"), "-1"))
  # fit <- geepack::geese(form, data = b, id=id, jack=TRUE, scale.fix=TRUE,
  # family = gaussian(link = "identity"), corstr="independence")
  fit2 <- geepack::geeglm(
    form, data = b, id=id, scale.fix=TRUE,
    family = gaussian(link = "identity"), corstr="independence")
  if (return_fit) return(fit2)
  
  {
    # Get Variance Estimators and p-values using summary function
    # fit.summary <- summary(fit)$mean
    # res <- fit.summary[(length(cutoffs)+1):nrow(fit.summary),]
    fit2.summary <- summary(fit2)$coefficients
    fit2.confint <- confint(fit2)
    res <- cbind(fit2.summary, fit2.confint)
    colnames(res) <- c('estimate', 'san.se', 'wald', 'p', '2.5%', '97.5%')
    
    fit2.anova <- anova(fit2)
    colnames(fit2.anova) <- paste0("anova.", colnames(fit2.anova))
    rownames(fit2.anova) <- NULL
    anova_rows <- c(rep(1, length(unique(b$tpseudo))))
    for (ind in seq_along(covar)){
      anova_rows <- c(
        anova_rows, rep(ind + 1, length(unique(b[,covar[ind],drop=T])) - 1)
        )
    }
    res <- cbind(res, fit2.anova[anova_rows,])
    # res <- res[(length(cutoffs)+1):nrow(res),]
  }
  class(res) <- 'egee'
  return(res) 
}
