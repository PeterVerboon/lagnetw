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
#' @import ggplot2
#' @return Estimate of difference between groups wrt network connectivity with p value based on permutations.
#'            and permutation distribution
#' @export
#'
#' @examples
#' data("gratitude")
#' vars <- c("pa_1","pa_2","pa_3","na_1","na_2","na_3")
#' out <- conDif(dat=gratitude,vars=vars, group="wellBeing", subjnr="subjnr",
#' level1="beepno", level2 = "dayno", randomVars = F, perms = 100) 
conDif <- function(dat, vars, group, subjnr, level1, level2 = NULL, randomVars = F, perms = 500) {
  
  res <- list(intermediate = list(),
              output = list());
  
  dat$group <- dat[,group] 
   if (is.null(level2)) {
     level2 <- "ones"
     dat$ones <- 1
   }
  
  dat1 <- lagnetw::centerESM(data = dat, subjnr = subjnr, addmeans = F, varnames = vars, center = "grand_mean")
  dat1 <- lagnetw::lagESM(data = dat1, subjnr = subjnr, level1 = level1, level2 = level2, lagn = 1, varnames = vars )
  
  k <- length(vars)  
  
  ### set up matrices to store coefficients in
  b1 <- matrix(NA, nrow=k, ncol=k)
  b2 <- matrix(NA, nrow=k, ncol=k)
  
  ### create names for lagged vars and random terms
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
 
  ## analyses of observed data
  for (i in 1:k) {
    
    ff <- stats::as.formula(paste0(vars[i], "~", pred1, sep=""))
    
    ### fit models
    res1 <- lme4::lmer(ff, data=subset(dat1, group == 1), REML=FALSE)
    res2 <- lme4::lmer(ff, data=subset(dat1, group == 2), REML=FALSE)
    
    ### store coefficients
    b1[i,] <- lme4::fixef(res1)[2:(k+1)]
    b2[i,] <- lme4::fixef(res2)[2:(k+1)]
    
  }
  
  res$output$fixedEffects1 <- b1
  res$output$fixedEffects2 <- b2
  
  ### Permutation analyses
  
  ### number of observations per person
  nobs.per.person  <- sapply(split(dat1$group, dat1$subjnr), length)
  
  ### group of each person
  group.per.person <- sapply(split(dat1$group, dat1$subjnr), function(x) x[1])
  
  ### repeatedly apply permfunc() function
  
  pb1 <- utils::txtProgressBar(min = 0, max = perms, style = 3)
  permres <- lapply(1:perms, permfunc, dat=dat1, pb = pb1, outnames=vars, pred=pred1, nobs.per.person=nobs.per.person, group.per.person=group.per.person)
  close(pb1)
  
  ### turn results into a matrix
  permres <- do.call(rbind, permres)
  
  
  ### table with permutation based p-values (two definitions of the p-values)
  
  b.diff.obs <- rep(NA_real_, 4)
  names(b.diff.obs) <- c("total diff", "diagonal diff", "offdiag diff", "SD diff")
  
  b.diff.obs[1] <- mean(abs(b1)) - mean(abs(b2))
  b.diff.obs[2] <- mean(abs(diag(b1))) - mean(abs(diag(b2)))
  
  ### set diagonal elements to NA (then can take the mean with na.rm=TRUE and get the mean of the off-diagonal elements)
  diag(b1) <- NA
  diag(b2) <- NA
  b.diff.obs[3] <- mean(abs(b1), na.rm=TRUE) - mean(abs(b2), na.rm=TRUE)
  b.diff.obs[4] <- stats::sd(b1, na.rm =TRUE) - stats::sd(b2, na.rm =TRUE)
  
  p.perm.def1 <- p.perm.def2 <- rep(NA, length(b.diff.obs))
  for (j in 1:length(b.diff.obs)) {
    p.perm.def1[j] <- 2*min(mean(permres[,j] >= b.diff.obs[j], na.rm=TRUE), mean(permres[,j] <= b.diff.obs[j], na.rm=TRUE))
    p.perm.def2[j] <- min(1, 2*ifelse(b.diff.obs[j] > 0, mean(permres[,j] >= b.diff.obs[j], na.rm=TRUE), mean(permres[,j] <= b.diff.obs[j], na.rm=TRUE)))
  }
  
  ###  results
  
  outp <- round(cbind(b.diff.obs, "p-perm.def1"=p.perm.def1, "p-perm.def2"=p.perm.def2), 4)
  
  res$output$resObserved1 <- res1
  res$output$resObserved2 <- res2
  res$output$pvals <- outp
  res$output$permutations <- permres <- data.frame(permres)
  
  colnames(b1) <- rownames(b1) <- colnames(b2) <- rownames(b2) <- vars

  
  
  class(res) <- "conDif"
  return(res)
  
} # end function



