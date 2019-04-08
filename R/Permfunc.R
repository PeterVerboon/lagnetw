#' Function to test the difference of the network connectivity between two groups
#'
#' @param perms number of permutations
#' @param dat data frame 
#' @param pb name progressbar
#' @param outnames vector with the names of the variables which are centered.
#' @param pred formula with names of predictor variables for lmer.
#' @param nobs.per.person vector with number of observations per person             
#' @param group.per.person vector with group number for each person             
#'
#' @return matrix (4 x 3) with total differences, diagonal differences, off-diagonal differences,
#'         and differences between standard deviations, with their p-values (2 definitions).
#' @export
#'
permfunc <- function(perms, dat, pb, outnames, pred, nobs.per.person, group.per.person) {

  setTxtProgressBar(pb, perms)

   ### set up matrices to store coefficients in
    k <- length(outnames)
    
    b1.perm <- matrix(NA, nrow=k, ncol=k)
    b2.perm <- matrix(NA, nrow=k, ncol=k)
   
  
   ### reshuffle group variable
   dat$group <- rep(sample(group.per.person), times=nobs.per.person)

   
   ### fit models to reshuffled data

   for (i in 1:k) {

      ff <- as.formula(paste0(outnames[i], "~", pred, sep=""))

      ### fit models
      res1.perm <- try(lme4::lmer(ff, data=subset(dat, group == 1), REML=FALSE), silent = TRUE) 
      res2.perm <- try(lme4::lmer(ff, data=subset(dat, group == 2), REML=FALSE), silent = TRUE)
      
      ### if one of the models doesn't converge, return NA; otherwise store coefficients

      if (inherits(res1.perm, "try-error") | inherits(res2.perm, "try-error")) {
         return(NA)
      } else {
         b1.perm[i,] <- lme4::fixef(res1.perm)[2:(k+1)]
         b2.perm[i,] <- lme4::fixef(res2.perm)[2:(k+1)]
         }
  
   }

   ### if we get this far, then all models converged and we can return the relevant statistics

   sav <- rep(NA_real_, 4)
   sav[1] <- c(mean(abs(b1.perm)) - mean(abs(b2.perm)))
   sav[2] <- c(mean(abs(diag(b1.perm))) - mean(abs(diag(b2.perm))))
   diag(b1.perm) <- NA
   diag(b2.perm) <- NA
   sav[3] <- mean(abs(b1.perm), na.rm=TRUE) - mean(abs(b2.perm), na.rm=TRUE)
   sav[4] <- sd(b1.perm, na.rm =TRUE) - sd(b2.perm, na.rm =TRUE)
   names(sav) <- c("total diff", "diagonal diff", "offdiag diff", "SD diff")
   return(sav)

} 





