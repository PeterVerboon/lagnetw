#' Test the difference of the network connectivity between two groups
#' 
#' The function uses permutations to obtain p-values of the connectivity differences, 
#' used in permDif function.
#'
#' @param perms number of permutations
#' @param dat data frame 
#' @param pb name progressbar
#' @param outnames vector with the names of the variables which are centered.
#' @param pred formula with names of predictor variables for lmer.
#' @param subset  subset of predictor variables which are compared in summary statistics. If null then result is also null.                  
#' @param nobs.per.person vector with number of observations per person             
#' @param group.per.person vector with group number for each person
#' @param type type of analyses: lagged ("lagged") or contemporaneous predictors ("contemp")
#' @param optim optimizer used in lmer, options: "bobyqa" or "Nelder_Mead", see lmerControl (lme4)           
             
#'
#' @return matrix (4 x 3) with total differences, diagonal differences, off-diagonal differences,
#'         and differences between standard deviations, with their p-values (2 definitions).
#' @export
#'
permfunc <- function(perms, dat, pb, outnames, pred,pnames, subset = NULL, 
                     nobs.per.person, group.per.person, type = "lagged", optim = "bobyqa") {

  utils::setTxtProgressBar(pb, perms)

   ### set up matrices to store coefficients in
    k <- length(outnames)

    b1.perm <- matrix(NA, nrow=k, ncol=k)
    b2.perm <- matrix(NA, nrow=k, ncol=k)
   
  
   ### reshuffle group variable
   dat$group <- rep(sample(group.per.person), times=nobs.per.person)

   datx <- dat
   
   ### fit models to reshuffled data

   for (i in 1:k) {

      ff <- stats::as.formula(paste0(outnames[i], "~", pred, sep=""))
      
     if (type == "contemp") {
        datx <- dat
        datx[,pnames[i]] <- rnorm(dim(datx)[1], 0, .5)
     }
      

      ### fit models 
      res1.perm <- try(lme4::lmer(ff, data=subset(datx, datx$group == 1), REML=FALSE, 
                                  control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE)), 
                       silent = TRUE) 
      res2.perm <- try(lme4::lmer(ff, data=subset(datx, datx$group == 2), REML=FALSE, 
                                  control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE)), 
                       silent = TRUE)
      
      ### if one of the models doesn't converge, return NA; otherwise store coefficients

      if (inherits(res1.perm, "try-error") | inherits(res2.perm, "try-error")) {
         return(NA)
      } else {
         b1.perm[i,] <- lme4::fixef(res1.perm)[2:(k+1)]               # the covariates are ignored here
         b2.perm[i,] <- lme4::fixef(res2.perm)[2:(k+1)]
         }
   
   }

   ### if we get this far, then all models converged and we can return the relevant statistics
   b1.permD <- b1.perm
   b2.permD <- b2.perm
   diag(b1.permD) <- NA
   diag(b2.permD) <- NA
   sav <- rep(NA_real_, 5)
   sav[1] <- c(mean(abs(b1.perm)) - mean(abs(b2.perm)))
   sav[2] <- c(mean(abs(diag(b1.perm))) - mean(abs(diag(b2.perm))))
   sav[3] <- mean(abs(b1.permD), na.rm=TRUE) - mean(abs(b2.permD), na.rm=TRUE)
   if (!is.null(subset)) {
      s <- (outnames %in% subset)*(1:k)
      sav[4] <- mean(abs(b1.perm[,s]), na.rm=TRUE) - mean(abs(b2.perm[,s]), na.rm=TRUE)
   } 
   sav[5] <- stats::sd(b1.perm, na.rm =TRUE) - stats::sd(b2.perm, na.rm =TRUE)
   names(sav) <- c("total ", "diagonal ", "off-diag ", "subset  ", "SD ")
   
   
   difFE <- b1.perm - b2.perm
   
  
   output <- list("FEsummary" = sav, "FEdifferences" = difFE, "perm1" = b1.perm, "perm2" = b2.perm)
   
  #   perm1 = k x k matrices with fixed effect in group 1 times permutations
   #  perm2 = k x k matrices with fixed effect in group 2 times permutations
   #  FEdifferences = k x k matrices with differences between fixed effect in group 1 and 2 times permutations
   #  FEsummary = vector with 4 summary differences between fixed effect in group 1 and 2 times permutations
   
   return(output)

} 









