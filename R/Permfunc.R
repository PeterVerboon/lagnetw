#' Test the difference of the network connectivity between two groups
#' 
#' The function uses permutations to obtain p-values of the connectivity differences.
#'
#' @param perms number of permutations
#' @param dat data frame 
#' @param pb name progressbar
#' @param outnames vector with the names of the variables which are centered.
#' @param pred formula with names of predictor variables for lmer.
#' @param nobs.per.person vector with number of observations per person             
#' @param group.per.person vector with group number for each person
#' @param optim optimizer used in lmer, options: "bobyqa" or "Nelder_Mead", see lmerControl (lme4)           
             
#'
#' @return matrix (4 x 3) with total differences, diagonal differences, off-diagonal differences,
#'         and differences between standard deviations, with their p-values (2 definitions).
#' @export
#'
permfunc <- function(perms, dat, pb, outnames, pred, nobs.per.person, group.per.person, optim = "bobyqa") {

  utils::setTxtProgressBar(pb, perms)

   ### set up matrices to store coefficients in
    k <- length(outnames)
    
    b1.perm <- matrix(NA, nrow=k, ncol=k)
    b2.perm <- matrix(NA, nrow=k, ncol=k)
   
  
   ### reshuffle group variable
   dat$group <- rep(sample(group.per.person), times=nobs.per.person)

   
   ### fit models to reshuffled data

   for (i in 1:k) {

      ff <- stats::as.formula(paste0(outnames[i], "~", pred, sep=""))

      ### fit models 
      res1.perm <- try(lme4::lmer(ff, data=subset(dat, dat$group == 1), REML=FALSE, 
                                  control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE)), 
                       silent = TRUE) 
      res2.perm <- try(lme4::lmer(ff, data=subset(dat, dat$group == 2), REML=FALSE, 
                                  control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE)), 
                       silent = TRUE)
      
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
   sav[4] <- stats::sd(b1.perm, na.rm =TRUE) - stats::sd(b2.perm, na.rm =TRUE)
   names(sav) <- c("total diff", "diagonal diff", "offdiag diff", "SD diff")
   return(sav)

} 




#' Test the difference of the network paths between two groups
#' 
#' The function uses permutations to obtain p-values of the path differences.
#'
#' @param iter number of permutations
#' @param dat data frame 
#' @param pred vector with names of predictor variables
#' @param dv name dependent variable
#' @param group dichotomous variable that indicates the two groups to be compared
#' @param subjnr identification number of the subjects
#' @param nobs.per.person vector with number of observations per person             
#' @param group.per.person vector with group number for each person 
#' @param is.conti boolean to indicate if dependent variable is continuous or dichotomous (NOT ACTIVE)            
#'
#' @return matrix (4 x 3) with total differences, diagonal differences, off-diagonal differences,
#'         and differences between standard deviations, with their p-values (2 definitions).
#' @export
#'
permfunc2 <- function(iter, dat, pred, dv, group, subjnr, nobs.per.person, group.per.person, is.conti = TRUE) {

  
  ### reshuffle outcome variable within subjects
  ## dat$outcome <- unlist(sapply(split(dat[,dv], dat[,subjnr]), sample))
  
  ### reshuffle group variable
  dat$group <- rep(sample(group.per.person), times=nobs.per.person)
  
  ff <- stats::as.formula(paste0(dv," ~ ", pred, sep=""))
  
  ### fit model for both groups to data with reshuffled outcome
  if (is.conti) 
  {res.perm1 <- try(nlme::lme(ff, random = ~ 1 | subjnr, data=subset(dat, group == 1), na.action = stats::na.omit, 
                        control=list(opt="optim")), silent=TRUE)
  res.perm2 <- try(nlme::lme(ff, random = ~ 1 | subjnr, data=subset(dat, group == 2), na.action = stats::na.omit, 
                       control=list(opt="optim")), silent=TRUE)
  } else {
    res.perm1 <- try(lme4::glmer(outcome ~ pred, random = ~ 1 | subjnr, data=subset(dat, group = 1), family = stats::binomial), silent=TRUE)
    res.perm2 <- try(lme4::glmer(outcome ~ pred, random = ~ 1 | subjnr, data=subset(dat, group = 2), family = stats::binomial), silent=TRUE)
  }
  
  ### if model doesn't converge, return NA; otherwise return coefficients
  if (is.conti) {
  if (inherits(res.perm1, "try-error") | inherits(res.perm2, "try-error")) {
    return(difFixEffects <- NA)
     } else {
       return(difFixEffects <- nlme::fixef(res.perm1) - nlme::fixef(res.perm2))
       }
  } else {
  if (inherits(res.perm1, "try-error") | inherits(res.perm2, "try-error")) {
    return(difFixEffects <- NA)
     } else {
       return(difFixEffects <- lme4::fixef(res.perm1) - lme4::fixef(res.perm2))
       }
  }
  
}  # end permfunc2



# a <- permfunc2(iter=1,dat=dat1, pred = varsp, dv = outname,subjnr = subjnr,
#               group = group,nobs.per.person = nobs.per.person, group.per.person=group.per.person)




