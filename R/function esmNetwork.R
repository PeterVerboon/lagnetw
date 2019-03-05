
#' Function to compute network of ESM variables using lags in multilevel analysis
#'
#' @param dat data frame 
#' @param subjnr variable name, indicating the subjects
#' @param level1 variable name, indicating level 1 variable, e.g. beeps
#' @param level2 variable name, indicating level 2 variable, e.g. days
#' @param vars vector with the names of the variables for which the network is computed
#' @param covs covariates in the analyses, which are not plotted in the network
#' @param randomAll logical indicates whether all variables should be used as random effects
#' @param randomVars vector of variable names, used as random effects
#' @param groups variable used to label groups in the network figure
#' @param lagn number of lags used in the network
#' @param labs labels used in the network plot
#' @param solid effect size above which lines are shown as solid (default = .10)
#' @param plimit p-value under which lines are shown (default = .05)
#' @param titlePlot title for the plot
#'
#' @return a qgraph object (network) that can be plotted
#' @export
#'
#' @examples
#' data("DataNews")
#' res <- esmNetwork(dat = DataNews, subjnr="subjnr",level1 = "beepnr",
#'                   level2= "daynr", vars = vars, labs = labs, lagn = 1)
#'        
esmNetwork <- function(dat, subjnr, level1, level2 = NULL,  vars, covs = NULL, 
                       randomAll = FALSE, randomVars = NULL, 
                       groups = NULL, lagn=1, layout = "spring", labs=NULL, 
                       solid = .10, plimit = .05, titlePlot="Figure"){
  
  result <- list(input = as.list(environment()),
                 intermediate = list(),
                 output = list());
  
  result$intermediate$dataName <- as.character(deparse(substitute(dat)));
  
  nvars <- length(vars)                # number of variables involved in the network analyses
  npred <- length(covs) + nvars        # number of predictors involved in the analyses
  
  result$intermediate$numberOfVars <- nvars
  result$intermediate$numberOfPreds <- npred
  
  if (nvars < 3) { stop("Number of variables in the network should be more than 2") }
  
  if (is.null(level2)) {dat$level2 <- 1; level2 <- "level2"}
  
   dat1 <- dat[,c(subjnr,level2,level1, covs,vars)]
  
   if (is.null(labs)) { labs <- vars}

  # Vector of predictor names (lagged variables)
   
   varsp <- paste0(vars[1],"L",lagn)
   for (i in 2:length(vars)) {
     varsp <- paste0(varsp, " + ", vars[i],"L",lagn)
   }
   
  # Vector of predictor names with random effect (lagged variables)
   
   if (!is.null(randomVars)) { 
     vrandom <- 1 
       for (i in (seq_along(randomVars))) {
           vrandom <- paste0(vrandom, "+", randomVars[i],"L",lagn)
       }
     }

   covs2 <- ifelse(is.null(covs), "", paste0(paste0(covs,collapse = " + "), " + "))
   

   ## check random effects and build formula
   
   if (randomAll) {
     vrandom <- varsp
     pred1 <- paste0("(",covs2, varsp," + (", varsp, "|",subjnr,"))")
   } else {
     if (!exists("vrandom")) vrandom <- 1
     pred1 <- paste0("(",covs2, varsp," + (", vrandom, "|",subjnr,"))") 
   }
   
 
  ## Construct lagged variables
   
  dat2 <- LagESM(dat1, 
                 subjnr=subjnr,
                 level2=level2,
                 level1=level1, 
                 lagn=lagn, 
                 varnames=vars)

  
  ### run MLA for all variables in network
  
  model1 <- list()
  
  for (j in 1:nvars) {
    ff=as.formula(paste(vars[j],"~", pred1, sep="")); 
    model1[[j]]<-lme4::lmer(ff, data=dat2, REML=FALSE)
    print(j)
  }
  
  
  ##  inferring the coefficients for the network from the fitted model1
  
  coef1=data.frame(matrix(unlist(lapply(model1,lme4::fixef),use.names=FALSE),
                          byrow=TRUE, 
                          ncol=(npred+1))) 
  colnames(coef1)=names(lme4::fixef(model1[[1]]))
  rownames(coef1)=vars
  
  se.coef1=data.frame(matrix(unlist(lapply(model1,arm::se.fixef),use.names=FALSE),
                             byrow=TRUE,
                             ncol=(npred+1))) 
  colnames(se.coef1)=names(lme4::fixef(model1[[1]]))
  rownames(se.coef1)=vars

  E <- cbind(from=rep(1:nvars,each=nvars),
             to=rep(1:nvars,nvars),
             weigth=unlist(coef1[,(2+npred-nvars):(npred+1)]))
  pvals <- 2*(1-pnorm(abs(unlist(coef1[,(2+npred-nvars):(npred+1)]/se.coef1[,(2+npred-nvars):(npred+1)]) )))
  edge.color <- addTrans(ifelse(E[,3] > 0, "green3", "red3"), ifelse(pvals < plimit, 255, 0))
  
  result$intermediate$pvals <- pvals
  result$intermediate$edgeColor <- edge.color
  
  G1 <- qgraph::qgraph(E,fade=FALSE,
                       groups = groups,
                       layout=layout,
                       labels=labs, 
                       lty=ifelse(E[,3] > solid, 1, 5),
                       edge.labels=F,
                       edge.color=edge.color, 
                       title=titlePlot)
 
   ## centrality measures
  
  C <- qgraph::centrality(E, alpha = 1, posfun = abs, all.shortest.paths = FALSE)
  C <- data.frame(cbind(C$Betweenness, C$Closeness, C$OutDegree, C$InDegree))
  names(C) <- c("Betweenness", "Closeness"," outDegree", " inDegree")
  row.names(C) <- vars
  
  ##  In "VV" the individual differences are taken from the fitted model1, each link now indicates 
  ##  the amount of variability across subjects
  ##  Assumed that there are 1 (intercept) + nvars random effecten
  
  if (randomAll) { 
  
     VV <- sqrt(t(matrix(unlist(lapply(model1,function(x){VV=diag(lme4::VarCorr(x)[[1]][2:(nvars+1),2:(nvars+1)])})),nvars,nvars)))
  
     E2 <- cbind(from=rep(1:nvars,each=nvars),to=rep(1:nvars,nvars),weigth=as.vector(VV))
     
     edge.color <- addTrans("blue", ifelse(E[,3]>.095, 255, 20))
     G2 <- qgraph::qgraph(E2,fade=FALSE,
                          layout="circular",
                          labels=labs,
                          lty=ifelse(E[,3] > solid, 1, 5),
                          edge.labels=F,
                          edge.color=edge.color,
                          title = "Path variability")
     
     (result$output <- list(plotData = E, network = G1, centralityMeasures = C, pathVariability = G2, model = model1))
  
  } else {
        (result$output <- list(plotData = E, network = G1, centralityMeasures = C,  model = model1))
    }
  
 
  class(result) <- "esmNetwork"
  return(result)
  
  
  
}  # end function




