#' Test the difference of the network connectivity between two groups
#'
#' @param dat data frame
#' @param vars vector with the names of the dependent variables. 
#' @param covs vector with the names of the covariates. 
#' @param group dichotomous variable that indicates the two groups to be compared
#' @param subjnr identification number of the subjects
#' @param randomVars vector, indicating which variables should be included as random effects. 
#'                   If "all" then all fixed effects are taken. If "null" only intercept is used as random effect.
#' @param subset  subset of predictor variables which are compared in summary statistics. If null then result is also null.                  
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
#' randomVars = FALSE, perms = 10) 
permDif <- function(dat, vars, covs = NULL, group, subjnr, randomVars = NULL, subset = NULL,
                   perms = 500, optim = "bobyqa") {
  
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  resall <- list()
  
  dat$group <- dat[,group] 

    k <- length(vars)  
    if (is.null(covs)) p <- k
    if (!is.null(covs)) p <- k + length(covs)
  
  ### set up matrices to store coefficients in
  b.diff.obs <- rep(NA_real_, 5)
  names(b.diff.obs) <- c("total      ", "diagonal    ", "off-diag    ", "subset  ", "standard deviation")
  b1 <- matrix(NA, nrow=k, ncol=p)
  b2 <- matrix(NA, nrow=k, ncol=p)
  
  ### set up matrix to store the p values in (two definitions ) 
  p1.perm <-  rep(NA, length(b.diff.obs))
  p2.perm <-  matrix(NA, nrow = k, ncol = p)
  
  ### create names for lagged vars and random term 
  
  pnames <- varsp <- paste0(vars[1],"L",1)
    for (i in 2:length(vars)) {
        pnames <- c(pnames, paste0(vars[i],"L",1))
        varsp <- paste0(varsp, " + ", vars[i],"L",1)
    }
   if (!is.null(covs)) {
    varsp <- c(varsp, covs)
    pnames <- c(pnames, covs)
  }
  
  ## variables as random effect or intercept only
  if (randomVars == "all") {
    pred1 <- paste0((paste0(varsp, collapse = " + "))," + (",(paste0(varsp, collapse = " + ")), "|",subjnr,")")
  } 
  if (is.null(randomVars)) {
    pred1 <- paste0((paste0(varsp, collapse = " + "))," + (",1, "|",subjnr,")")
  }
  if (!is.null(randomVars) & !(randomVars == "all")) {
    if (sum(!(randomVars %in% names(dat))) > 0) {return("At least one random term is not present in the data")}
    pred1 <- paste0((paste0(varsp, collapse = " + "))," + (",randomVars, "|",subjnr,")")
  }
 
  ### first choose optimizer
  
  if (is.null(optim))  optim <- "bobyqa"
  if (!optim %in% c("Nelder_Mead", "bobyqa")) {
    warning("Optimizer not correctly specified. Default is used.")
    optim <- "bobyqa"
  }
  
  res$intermediate$formula <- stats::as.formula(paste0(vars[1], "~", pred1, sep=""))
  
  ## analyses of observed data
  for (i in 1:k) {
   
    ff <- stats::as.formula(paste0(vars[i], "~", pred1, sep=""))
    
    ### fit models
    res1 <- lme4::lmer(ff, data=subset(datx, group == 1), REML=FALSE, 
                       control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE))
    res2 <- lme4::lmer(ff, data=subset(datx, group == 2), REML=FALSE, 
                       control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE))
    
    ### store coefficients
    b1[i,] <- lme4::fixef(res1)[2:(p+1)] 
    b2[i,] <- lme4::fixef(res2)[2:(p+1)]
    
  }
  
  rownames(p2.perm) <- rownames(b1) <- rownames(b2) <- vars
  colnames(p2.perm) <- colnames(b1) <- colnames(b2) <- pnames
  res$output$fixedEffects1 <- round(b1,4)
  res$output$fixedEffects2 <- round(b2,4)
  difs <- b1 - b2
  res$output$observedDifs <- round(difs, 4) 
 
   
  ## store differences in observed statistics
  
  b.diff.obs[1] <- mean(abs(b1)) - mean(abs(b2))
  b.diff.obs[2] <- mean(abs(diag(b1))) - mean(abs(diag(b2)))
  
  ## next only when a subset of the predictors is requested
  if (!is.null(subset)) {
    s <- (vars %in% subset)*(1:k)
    b.diff.obs[4] <- mean(abs(b1[,s]), na.rm=TRUE) - mean(abs(b2[,s]), na.rm=TRUE)
    res$intermediate$subset1 <- b1[,s]
    res$intermediate$subset2 <- b2[,s]
    res$intermediate$dif.subset <- abs(b1[, s])- abs(b2[, s])
  }
  
  ## set diagonal elements to NA (then can take the mean with na.rm=TRUE and get the mean of the off-diagonal elements)
  diag(b1) <- NA
  diag(b2) <- NA
  b.diff.obs[3] <- mean(abs(b1), na.rm=TRUE) - mean(abs(b2), na.rm=TRUE)
  
  ## and finally the differences in standard deviations
  b.diff.obs[5] <- stats::sd(b1, na.rm =TRUE) - stats::sd(b2, na.rm =TRUE)
  
  
  
  ### Permutation analyses
  
  ### number of observations per person
  nobs.per.person  <- sapply(split(dat$group, dat$subjnr), length)
  
  ### group of each person
  group.per.person <- sapply(split(dat$group, dat$subjnr), function(x) x[1])
  
  ### repeatedly apply permfunc() function
  
  pb1 <- utils::txtProgressBar(min = 0, max = perms, style = 3)
  permres <- lapply(1:perms, permfunc, dat=dat, pb = pb1, outnames=vars, 
                    pred=pred1, subset = subset, nobs.per.person=nobs.per.person, group.per.person=group.per.person)
  close(pb1)
  
  
  ### turn results into an array
  permres <- do.call(Map, c(rbind, permres))
  permres1 <- permres$FEsummary
  permres2 <- array(permres$FEdifferences,dim = c(k,perms,k))
  
  ### table with permutation based p-values 
  
  for (j in 1:length(b.diff.obs)) {
    p1.perm[j] <- 2*min(mean(permres1[,j] >= b.diff.obs[j], na.rm=TRUE), 
                            mean(permres1[,j] <= b.diff.obs[j], na.rm=TRUE))
  }
  
  
 #  
 #  ### table with permutation based p-values 
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
  res$output$pvalues.all <- round(p2.perm, 3)[,1:k]
  res$output$plimit_adjusted <- noquote(paste0("Bonferroni corrected alpha level: ", round(0.05/(k**2),4)))
  
  
  class(res) <- "permDif"
  return(res)
  
} # end function



