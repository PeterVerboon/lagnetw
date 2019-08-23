
#' Test the difference of the network path strenghts between two groups
#' 
#' Function that uses resampling method (permutations) to obtain distribion and p-values for
#' the differences of the path strengths of a network between two groups.
#'
#' @param dat data frame
#' @param vars variables in network 
#' @param group dichotomous variable that indicates the two groups to be compared
#' @param subjnr identification number of the subjects
#' @param level variable indicating level (e.g. days or beeps)
#' @param randomVars boolean indicating whether the variables are random effects in the model
#' @param perms number of permutations
#'
#' @return matrix (k x k) with observed differences of the paths ( k = number of variables in network), 
#'         2 matrices (k x k) with the p-values (2 definitions) obtained from permutation distribution for the differences.
#' @export
#'
#' @examples 
#' data("gratitude")
#' vars <- c("pa_1","pa_2","pa_3","na_1","na_2","na_3")
#' sigConnect(dat = gratitude, vars = vars , group="wellBeing", subjnr ="subjnr", 
#' level = "beepno", randomVars = FALSE, perms = 10)

sigConnect <- function(dat, vars, group, subjnr, level, randomVars = F,  perms = 500) {
  
  result <- list(input = as.list(environment()),
                 intermediate = list(),
                 output = list());
  
  
dat1 <- lagnetw::centerESM(data = dat, subjnr = subjnr, addmeans = F, varnames = vars, center = "grand_mean")
dat1 <- lagnetw::lagESM(data = dat1, subjnr = subjnr,  level1 = level, lagn = 1, varnames = vars )

dat1$group <- dat1[,group]

k <- length(vars)  

### set outcomes to loop through
outnames <- vars

### create names for lagged vars and random terms
varsp <- paste0(vars[1],"L",1)
for (i in 2:k) {
  varsp <- paste0(varsp, " + ", vars[i],"L",1)
}
## variables as random effect or intercept only
if (randomVars == TRUE) {
  randomPred1 <- paste0( varsp, "|",subjnr)
} else {
  randomPred1 <- stats::as.formula(paste0("~", 1, "|", subjnr))
}


### set up matrix to store the differences between the coefficients in
difs <- matrix(NA, nrow=k, ncol=k+1)

### set up matrix to store the mean differences (across permutations) for each dependent variable in
meanDifs <- matrix(NA, nrow=k, ncol=k+1)

### set up matrix to store the p values in (two definitions ) 
p.perm.def1 <- p.perm.def2 <- matrix(NA, nrow = k, ncol = k+1)

### number of observations per person
nobs.per.person  <- sapply(split(dat1$group, dat1$subjnr), length)

### group of each person
group.per.person <- sapply(split(dat1$group, dat1$subjnr), function(x) x[1])

### loop over groups and variables 
   i <- 0
   for (outname in outnames) {
     i <- i + 1
         
      cat("Outcome:", outname, "\n")

      ff <- stats::as.formula(paste0(outname, "~", varsp, sep=""))

      is.conti <- TRUE  ## not yet dichotomous variables
      
      ### fit model with actual data
      if (is.conti) {
         res1 <- nlme::lme(ff, random = randomPred1, data = subset(dat1, group == 1), na.action= stats::na.omit, control=list(opt="optim"))
         res2 <- nlme::lme(ff, random = randomPred1, data = subset(dat1, group == 2), na.action= stats::na.omit, control=list(opt="optim"))
      } else {
         res1 <- lme4::glmer(ff, random = randomPred1, data = subset(dat1, group == 1), family= stats::binomial)
         res2 <- lme4::glmer(ff, random = randomPred1, data = subset(dat1, group == 2), family= stats::binomial)
      }

         difs[i,] <- lme4::fixef(res1) - lme4::fixef(res2)
         
      
      ### repeatedly apply permfunc() function
      permres <- lapply(1:perms, permfunc2, dat=dat1, pred = varsp, dv = outname,subjnr = subjnr,
                        group = group,nobs.per.person = nobs.per.person, group.per.person=group.per.person, is.conti=TRUE)
     
      ### turn results into a matrix
      permres <- do.call(rbind, permres)
      
      meanDifs[i,] <- apply(permres, 2, mean)
     
      ### table with model-based and permutation based p-values (two definitions of the p-values)
      for (j in 1:(k+1)) {
         p.perm.def1[i,j] <- 2*min(mean(permres[,j] >= difs[i,j], na.rm=TRUE), 
                                 mean(permres[,j] <= difs[i,j], na.rm=TRUE))
         p.perm.def2[i,j] <- min(1, 2*ifelse(difs[i,j] > 0, 
                                           mean(permres[,j] >= difs[i,j], na.rm=TRUE), 
                                           mean(permres[,j] <= difs[i,j], na.rm=TRUE)))
      }
      
   }
   
   colnames(difs) <- colnames(meanDifs) <- c("intercept", paste0(vars,"L1") )
   rownames(difs) <- rownames(meanDifs) <- vars
   colnames(p.perm.def1) <- colnames(p.perm.def2) <- c("intercept", paste0(vars,"L1") )
   rownames(p.perm.def1) <- rownames(p.perm.def2) <- vars
   
  result$output$observedDifs <- difs
  result$output$meanDifs <- meanDifs
  result$output$pvalues.def1 <- p.perm.def1
  result$output$pvalues.def2 <- p.perm.def2
  
  return(result)
}



