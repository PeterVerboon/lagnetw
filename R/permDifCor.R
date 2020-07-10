#' Test the difference of the network connectivity between two groups based on the correlations
#'
#' Note that these variables must exist in the data.
#' 
#' @param dat data frame
#' @param vars vector with the names of the dependent variables. 
#' @param group dichotomous variable that indicates the two groups to be compared
#' @param perms number of permutations.             
#' @return Estimate of difference between groups wrt correlation based network connectivity with p value based on permutations.
#'          The p-values for differences of all individual correlations and of summaries are given.
#'          Also p-values for the differences of the centrality measures outDegree and closeness and betweenness are given.
#' @export
#'
#' @examples
#' data("gratitude")
#' vars <- c("pa_1","pa_2","pa_3","na_1","na_2","na_3")
#' out <- permDifCor(dat=gratitude,vars=vars, group="wellBeing", perms = 10)
 
permDifCor <- function(dat, vars,  group, perms = 1000) {
  
  res <- list(input = as.list(environment()),
              intermediate = list(),
              output = list());
  
  dat$group <- dat[,group]
  s <- dat[,group] == min(dat[,group])
  dat$group[s]  <- 1
  dat$group[!s]  <- 2

    k <- length(vars)  
    m <- k*(k-1)/2

  ### set up matrix to store the p values in (two definitions ) 
 
  obsCentrality1 <-  obsCentrality2 <- matrix(NA, nrow=3, ncol = k)
  colnames(obsCentrality1) <- colnames(obsCentrality2) <- vars
  
  difPerm <- perm1 <- perm2 <-  matrix(NA, nrow = perms, ncol = m)  ## correlations (network weights)
  meanDifPerm <-  rep(NA, perms)
  permCentrality1outDegree <- permCentrality2outDegree <- matrix(NA, nrow=perms, ncol=k)
  permCentrality1closeness <- permCentrality2closeness <- matrix(NA, nrow=perms, ncol=k)
  permCentrality1between <- permCentrality2between <- matrix(NA, nrow=perms, ncol=k)
  colnames(permCentrality1outDegree) <- colnames(permCentrality2outDegree) <- vars
  colnames(permCentrality1closeness) <- colnames(permCentrality2closeness) <- vars
  colnames(permCentrality1between) <- colnames(permCentrality2between) <- vars
 
  ## observed networks
  
  dat1 <- subset(dat, dat$group == 1)
  dat2 <- subset(dat, dat$group == 2)
  ## analyses of observed data
  res1 <- qgraph::qgraph(cor(dat1[,vars], use = "complete.obs"), graph = "cor", DoNotPlot = TRUE)
  b1 <- res1$Edgelist$weight
  res2 <- qgraph::qgraph(cor(dat2[,vars], use = "complete.obs"), graph = "cor", DoNotPlot = TRUE)
  b2 <- res2$Edgelist$weight
  
  difObs <- b1 - b2
  meanDifObs <- mean(abs(b1 - b2))
  res$output$observedDifs <- round(difObs, 4) 
  res$output$observedMeanDif <- round(meanDifObs, 4) 
  res$output$network1 <- res1 
  res$output$network2 <- res2
 
  obsCentrality1[1,] <- (qgraph::centrality(res1))$OutDegree
  obsCentrality2[1,] <- (qgraph::centrality(res2))$OutDegree
  obsCentrality1[2,] <- (qgraph::centrality(res1))$Closeness
  obsCentrality2[2,] <- (qgraph::centrality(res2))$Closeness
  obsCentrality1[3,] <- (qgraph::centrality(res1))$Betweenness
  obsCentrality2[3,] <- (qgraph::centrality(res2))$Betweenness
   
  
  ### Permutation analyses
  
  ### number of observations per person
  nobs.per.person  <- sapply(split(dat$group, dat$subjnr), length)
  
  ### group of each person
  group.per.person <- sapply(split(dat$group, dat$subjnr), function(x) x[1])
  
 
   ### repeatedly apply permutations function
   ### by reshuffling group variable
  
   for (i in 1:perms) {
     datx <- dat 
     datx$group <- rep(sample(group.per.person), times=nobs.per.person)
     datx1 <- subset(datx, datx$group == 1)
     datx2 <- subset(datx, datx$group == 2)
     
     resx1 <- qgraph::qgraph(cor(datx1[,vars], use = "complete.obs"), graph = "cor", DoNotPlot = TRUE)
     perm1[i,] <- resx1$Edgelist$weight
     resx2 <- qgraph::qgraph(cor(datx2[,vars], use = "complete.obs"), graph = "cor", DoNotPlot = TRUE)
     perm2[i,] <- resx2$Edgelist$weight
     difPerm[i,] <- resx1$Edgelist$weight - resx2$Edgelist$weight
     meanDifPerm[i] <- mean(abs(resx1$Edgelist$weight - resx2$Edgelist$weight))
     permCentrality1outDegree[i,] <- (qgraph::centrality(resx1))$OutDegree
     permCentrality2outDegree[i,] <- (qgraph::centrality(resx2))$OutDegree
     permCentrality1closeness[i,] <- (qgraph::centrality(resx1))$Closeness
     permCentrality2closeness[i,] <- (qgraph::centrality(resx2))$Closeness
     permCentrality1between[i,] <- (qgraph::centrality(resx1))$Betweenness
     permCentrality2between[i,] <- (qgraph::centrality(resx2))$Betweenness
   }

  res$intermediate$permutation1 <- perm1
  res$intermediate$permutation2 <- perm2
  res$intermediate$difPermutation <- difPerm
  res$intermediate$meanDifPermutation <- meanDifPerm
  
 
    ### permutation based p-values 
  
    a <- rbind(apply((difPerm >= difObs),2,mean), apply((difPerm <= difObs),2,mean))
    res$output$pDif <- 2*apply(a,2,min)
    res$output$pMeanDif <- 2*min(mean(meanDifPerm >= meanDifObs, na.rm=TRUE), mean(meanDifPerm <= meanDifObs, na.rm=TRUE))
  
  
  #  ### permutation based p-values for centrality measures

    a <- rbind(apply((permCentrality1outDegree >= obsCentrality1[1,]),2,mean), 
               apply((permCentrality1outDegree <= obsCentrality1[1,]),2,mean))
    res$output$pOutDegree1 <- 2*apply(a,2,min)
    a <- rbind(apply((permCentrality2outDegree >= obsCentrality2[1,]),2,mean), 
               apply((permCentrality2outDegree <= obsCentrality2[1,]),2,mean))
    res$output$pOutDegree2 <- 2*apply(a,2,min)
    
    a <- rbind(apply((permCentrality1closeness >= obsCentrality1[2,]),2,mean), 
               apply((permCentrality1closeness <= obsCentrality1[2,]),2,mean))
    res$output$pCloseness1 <- 2*apply(a,2,min)
    a <- rbind(apply((permCentrality2closeness >= obsCentrality2[2,]),2,mean), 
               apply((permCentrality2closeness <= obsCentrality2[2,]),2,mean))
    res$output$pCloseness2 <- 2*apply(a,2,min)
    
    a <- rbind(apply((permCentrality1between >= obsCentrality1[3,]),2,mean), 
               apply((permCentrality1between <= obsCentrality1[3,]),2,mean))
    a <- 2*apply(a,2,min)
    a[1 < a] <- 1
    res$output$pBetweenness1 <- a
    a <- rbind(apply((permCentrality2between >= obsCentrality2[3,]),2,mean), 
               apply((permCentrality2between <= obsCentrality2[3,]),2,mean))
    a <- 2*apply(a,2,min)
    a[1 < a] <- 1
    res$output$pBetweenness2 <-  a
  
 
  ###  results
  
  class(res) <- "permDifCor"
  return(res)
  
} # end function



