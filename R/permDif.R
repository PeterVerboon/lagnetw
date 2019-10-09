#' Test the difference of the network connectivity between two groups
#'
#' @param dat data frame
#' @param vars vector with the names of the variables which are centered. 
#' @param group dichotomous variable that indicates the two groups to be compared
#' @param subjnr identification number of the subjects
#' @param level1 variable, indicating level 1 variable, e.g. beeps. Must be specified.
#' @param level2 variable, indicating level 2 variable, e.g. days. Can be left empty. 
#' @param randomVars logical, indicating whether the variables should be included as random effects.
#' @param perms number of permutations.             
#' @param optim optimizer used in lmer, options: "bobyqa" or "Nelder_Mead", see lmerControl (lme4)           
#' @import ggplot2
#' @return Estimate of difference between groups wrt network connectivity with p value based on permutations.
#'          and permutation distribution. The p-values for differences of all individual paths and of summaries are given.
#' @export
#'
#' @examples
#' data("gratitude")
#' vars <- c("pa_1","pa_2","pa_3","na_1","na_2","na_3")
#' out <- permDif(dat=gratitude,vars=vars, group="wellBeing", subjnr="subjnr",
#' level1="beepno", level2 = "dayno", randomVars = FALSE, perms = 10) 
permDif <- function(dat, vars, group, subjnr, level1, level2 = NULL, randomVars = F, 
                   perms = 500, optim = "bobyqa") {
  
  res <- list(intermediate = list(),
              output = list());
  resall <- list()
  
  dat$group <- dat[,group] 
   if (is.null(level2)) {
     level2 <- "ones"
     dat$ones <- 1
   }
  
  dat1 <- lagnetw::centerESM(data = dat, subjnr = subjnr, addmeans = F, varnames = vars, center = "grand_mean")
  dat1 <- lagnetw::lagESM(data = dat1, subjnr = subjnr, level1 = level1, level2 = level2, lagn = 1, varnames = vars )
  
  k <- length(vars)  
  
  ### set up matrices to store coefficients in
  b.diff.obs <- rep(NA_real_, 4)
  names(b.diff.obs) <- c("total      ", "diagonal    ", "off-diag    ", "standard deviation")
  b1 <- matrix(NA, nrow=k, ncol=k)
  b2 <- matrix(NA, nrow=k, ncol=k)
  
  ### set up matrix to store the p values in (two definitions ) 
  p1.perm <-  rep(NA, length(b.diff.obs))
  p2.perm <-  matrix(NA, nrow = k, ncol = k)
  
  ### create names for lagged vars and random term
  varsp <- paste0(vars[1],"L",1)
  for (i in 2:length(vars)) {
    varsp <- paste0(varsp, " + ", vars[i],"L",1)
  }
 
  ## variables as random effect or intercept only
  if (randomVars == TRUE) {
    pred1 <- paste0("(", varsp," + (", varsp, "|",subjnr,"))")
  } else {
    pred1 <- paste0("(", varsp," + (", 1, "|",subjnr,"))")
  }
 
  ### first choose optimizer
  
  if (is.null(optim))  optim <- "bobyqa"
  if (!optim %in% c("Nelder_Mead", "bobyqa")) {
    warning("Optimizer not correctly specified. Default is used.")
    optim <- "bobyqa"
  }
  
  ## analyses of observed data
  for (i in 1:k) {
    
    ff <- stats::as.formula(paste0(vars[i], "~", pred1, sep=""))
    
    ### fit models
    res1 <- lme4::lmer(ff, data=subset(dat1, group == 1), REML=FALSE, 
                       control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE))
    res2 <- lme4::lmer(ff, data=subset(dat1, group == 2), REML=FALSE, 
                       control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE))
    
    ### store coefficients
    b1[i,] <- lme4::fixef(res1)[2:(k+1)] 
    b2[i,] <- lme4::fixef(res2)[2:(k+1)]
    
  }
  
  colnames(p2.perm) <- rownames(p2.perm) <- colnames(b1) <- rownames(b1) <- colnames(b2) <- rownames(b2) <- vars
  res$output$fixedEffects1 <- b1
  res$output$fixedEffects2 <- b2
  res$output$observedDifs <- difs <- b1 - b2
  
  ### Permutation analyses
  
  ### number of observations per person
  nobs.per.person  <- sapply(split(dat1$group, dat1$subjnr), length)
  
  ### group of each person
  group.per.person <- sapply(split(dat1$group, dat1$subjnr), function(x) x[1])
  
  ### repeatedly apply permfunc() function
  
  pb1 <- utils::txtProgressBar(min = 0, max = perms, style = 3)
  permres <- lapply(1:perms, permfunc, dat=dat1, pb = pb1, outnames=vars, 
                    pred=pred1, nobs.per.person=nobs.per.person, group.per.person=group.per.person)
  close(pb1)
  
  
  ### turn results into an array
  permres <- do.call(Map, c(rbind, permres))
  permres1 <- permres$FEsummary
  permres2 <- array(permres$FEdifferences,dim = c(6,perms,6))
  
  ### table with permutation based p-values (two definitions of the p-values)
  
  b.diff.obs[1] <- mean(abs(b1)) - mean(abs(b2))
  b.diff.obs[2] <- mean(abs(diag(b1))) - mean(abs(diag(b2)))
  
  ### set diagonal elements to NA (then can take the mean with na.rm=TRUE and get the mean of the off-diagonal elements)
  diag(b1) <- NA
  diag(b2) <- NA
  b.diff.obs[3] <- mean(abs(b1), na.rm=TRUE) - mean(abs(b2), na.rm=TRUE)
  b.diff.obs[4] <- stats::sd(b1, na.rm =TRUE) - stats::sd(b2, na.rm =TRUE)
  
  for (j in 1:length(b.diff.obs)) {
    p1.perm[j] <- 2*min(mean(permres1[,j] >= b.diff.obs[j], na.rm=TRUE), 
                            mean(permres1[,j] <= b.diff.obs[j], na.rm=TRUE))
  }
  
  
 #  
 #  ### table with model-based and permutation based p-values 
   for (i in 1:k) {
     for (j in 1:k) {
     p2.perm[i,j] <- 2*min(mean(permres2[i,,j] >= difs[i,j], na.rm=TRUE),
                                mean(permres2[i,,j] <= difs[i,j], na.rm=TRUE))
     }
   }
 # 
  ###  results
  
  res$output$permutations <- permres1
  res$output$pvalues.summary <- round(cbind("difference" = b.diff.obs, "p_value" = p1.perm), 4)
  res$output$pvalues.all <- round(p2.perm, 3)
  res$output$plimit_adjusted <- noquote(paste0("Bonferroni corrected alpha level: ", round(0.05/(k**2),4)))
  
  
  class(res) <- "permDif"
  return(res)
  
} # end function



