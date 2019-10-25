
#' Compute network of ESM variables using lags 
#' 
#' Based om lagged variables from an ESM study a network is constructed.
#'
#' @param dat data frame 
#' @param subjnr variable name, indicating the subjects
#' @param level1 variable name, indicating level 1 variable, e.g. beeps
#' @param level2 variable name, indicating level 2 variable, e.g. days
#' @param vars vector with the names of the variables for which the network is computed
#' @param covs covariates in the analyses, which are not plotted in the network. Covariates are not lagged and not centered.
#' @param randomAll logical indicates whether all variables should be used as random effects
#' @param randomVars vector of variable names, used as random effects
#' @param randomIcept logical indicating whether there is a random intercept (TRUE) or not (FALSE)
#' @param fixedIcept logical indicating whether there is a fixed intercept (TRUE) or not (FALSE)
#' @param groups variable used to label groups of variables in the network figure
#' @param lagn number of lags used in the network
#' @param centerType vector of same length as vars that indicates for each variable which centering must be,
#'               applied, either "person" or "grand_mean" or a number. If only "person" or "grand_mean" is specified (default), 
#'               all variables are person or grand-mean centered, respectively.
#'               If centerType is NULL no variables are centered.
#' @param optim optimizer used in lmer, options: "bobyqa" or "Nelder_Mead", see lmerControl (lme4)           
#' @param layout layout specification of plot (see options in qgraph)
#' @param labs labels used in the network plot
#' @param solid effect size above which lines are shown as solid (default = .10)
#' @param plimit p-value under which lines are shown (default = .05)
#' @param titlePlot title for the plot
#'
#' @return a qgraph object (network) that can be plotted and the centrality measures
#' @export
#'
#' @examples
#' data("news")
#' vars <- c("Fearful","Hopeful","Anxious","Down","Irritated","Relaxed","Insecure")
#' labs <- c("FF","HO","ANX","DOW", "IRR","REL","INS")
#' res <- esmNetwork(dat = news, subjnr="subjnr",level1 = "beepnr",
#'                   level2= "daynr", vars = vars, labs = labs, lagn = 1)
#'        
esmNetwork <- function(dat, subjnr, level1, level2 = NULL,  vars, covs = NULL, 
                       randomAll = FALSE, randomVars = NULL, 
                       randomIcept = TRUE, fixedIcept = TRUE,
                       groups = NULL, lagn=1, centerType = "person",
                       optim = "bobyqa",
                       layout = "spring", labs=NULL, 
                       solid = .10, plimit = .05, titlePlot="Figure"){
  
  result <- list(input = as.list(environment()),
                 intermediate = list(),
                 output = list());
  
  result$intermediate$dataName <- as.character(deparse(substitute(dat)));
  
  
  if (!is.null(centerType)) {
    if (sum(!centerType %in% c("person","grand_mean") ) > 0) 
      return(cat("centerType is not correctly specified","\n"))
    if (length(vars) != length(centerType) & length(centerType) > 1)  
      return(cat("Number of variables in varnames does not equal length of center","\n")) 
  }
  
  if(randomAll == FALSE & is.null(randomVars) == TRUE & randomIcept == FALSE) {
     return(cat("Analysis not run, because no random effects were specified.","\n")) }
  
  nvars <- length(vars)                # number of variables involved in the network analyses
  npred <- length(covs) + nvars        # number of predictors involved in the analyses
  allpred <- c(vars, covs)
  
  result$intermediate$numberOfVars <- nvars
  result$intermediate$numberOfPreds <- npred
  
  if (nvars < 3) { stop("Number of variables in the network should be more than 2","\n") }

  if (is.null(level2)) {dat$level2 <- 1; level2 <- "level2"}
  
  ## Construct lagged variables
  
  dat1 <- lagESM(dat, 
                 subjnr=subjnr,
                 level2=level2,
                 level1=level1, 
                 lagn=lagn, 
                 varnames=vars)
  
  result$intermediate$workData1 <- dat1

    # center data
 
  if (!is.null(centerType)) {
     dat2 <- centerESM(data = dat1,
                       subjnr = subjnr,
                       addmeans = TRUE,
                       varnames = vars,
                       center = centerType)
  } else {
     dat2 <- dat1
  }

  
  result$intermediate$workData2 <- dat2
  
  
  if (is.null(labs)) { labs <- vars}

  # Vector of predictor names (lagged variables)
   
   varsp <- paste0(vars[1],"L",lagn)
   for (i in 2:length(vars)) {
     varsp <- paste0(varsp, " + ", vars[i],"L",lagn)
   }
   
  # Vector of predictor names with random effect (lagged variables)
   
   if (!is.null(randomVars)) { 
     vrandom <- ifelse(randomIcept, 1, 0)
       for (i in (seq_along(randomVars))) {
          if (randomVars[i] %in% covs) {
            vrandom <- paste0(vrandom, " + ", randomVars[i])
          } else {
            vrandom <- paste0(vrandom, " + ", randomVars[i],"L",lagn)
          }
          
       }
     }

   covs2 <- ifelse(is.null(covs), "", paste0(paste0(covs,collapse = " + "), " + "))
   

   ## check random effects and build formula
   
   if (randomAll) {
     vrandom <- varsp
     pred1 <- paste0("(",covs2, varsp," + (", varsp, " | ",subjnr,"))")
   } else {
     if (!exists("vrandom")) vrandom <- 1
     pred1 <- paste0("(",covs2, varsp," + (", vrandom, " | ",subjnr,"))") 
   }
   
   if (fixedIcept == FALSE) pred1 <- paste0("-1 + ", pred1) 
 

  
  ### run MLA for all variables in network
  
  ### first choose optimizer
  
  if (is.null(optim))  optim <- "bobyqa"
  if (!optim %in% c("Nelder_Mead", "bobyqa")) {
    warning("Optimizer not correctly specified. Default is used.")
    optim <- "bobyqa"
  }
  
  
  model1 <- list()
  
  for (j in 1:nvars) {
    ff = stats::as.formula(paste(vars[j],"~", pred1, sep="")); 
    model1[[j]]<-lme4::lmer(ff, data=dat2, REML=FALSE, 
                            control = lme4::lmerControl(optimizer = optim, calc.derivs = FALSE))
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
  pvals <- 2*(1- stats::pnorm(abs(unlist(coef1[,(2+npred-nvars):(npred+1)]/se.coef1[,(2+npred-nvars):(npred+1)]) )))
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
  coef2 <- coef1[-1] - diag(as.matrix(coef1[-1]))
  Ostrength <- (apply(abs(coef2[,(npred-nvars+1):npred]), 2, sum))
  Istrength <- (apply(abs(coef2[,(npred-nvars+1):npred]), 1, sum))
  C <- qgraph::centrality(E, alpha = 1, posfun = abs, all.shortest.paths = FALSE)
  C <- round(data.frame(cbind(Ostrength,Istrength, C$Betweenness, C$Closeness, C$OutDegree, C$InDegree)),3)
  names(C) <- c("OutStrength","InStrength", "Betweenness", "Closeness"," outDegree", " inDegree")
  row.names(C) <- vars
  
  coef1 <- round(coef1, 4)
  
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
     
     (result$output <- list(plotData = E, network = G1, centralityMeasures = C, pathVariability = G2,
                            coef = coef1, model = model1))
  
  } else {
        (result$output <- list(plotData = E, network = G1, centralityMeasures = C,  
                               coef = coef1, model = model1))
    }
  
 
  
  class(result) <- "esmNetwork"
  return(result)
  
  
  
}  # end function




