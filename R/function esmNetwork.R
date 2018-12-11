
#' Function to compute network of ESM variables using lags in multilevel analysis
#'
#' @param dat data frame 
#' @param subjnr variable, indicating subjects
#' @param daynr variable, indicating days
#' @param beepnr variable, indicating beeps
#' @param vars vector with the names of the variables for which the network is computed
#' @param covs covariates in the analyses, which are not plotted in the network
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
#' res <- esmNetwork(dat = DataNews, subjnr="subjnr", daynr= "daynr",
#'        beepnr = "beepnr", vars = vars, labs = labs, lagn = 1)
#'        
esmNetwork <- function(dat, subjnr, daynr = NULL, beepnr, vars, covs=NULL,
                       lagn=1, labs=NULL, solid = .10, plimit = .05, titlePlot="Figure"){
  
   if (is.null(daynr)) {dat$daynr <- 1; daynr <- "daynr"}
  
   dat1 <- dat[,c(subjnr,daynr,beepnr, covs,vars)]
  
   if (is.null(labs)) { labs <- vars}

    # Vector of predictor names (lagged variables)
   
   varsp <- paste0(vars[1],"L",lagn)
   for (i in 2:length(vars)) {
     varsp <- paste0(varsp, " + ", vars[i],"L",lagn)
   }
   
   covs2 <- ifelse(is.null(covs), "", paste0(covs," + "))
   
   # all predictors, only random intercept
   pred1 <- paste0("(",covs2, varsp," + (1|",subjnr,"))")
   
   nvars = length(vars)                 # number of variables involved in the network analyses
   npred = length(covs) + nvars        # number of predictors involved in the analyses
   
 
  ### Construct lagged variables
   
  dat2 <- LagESM(dat1, subjnr=subjnr,daynr=daynr,beepnr=beepnr, lagn=lagn, 
                  varnames=vars)
  
  model1=list()
  
  ### run MLA for all variables in network
  
  for (j in 1:nvars) {
    ff=as.formula(paste(vars[j],"~", pred1, sep="")); 
    model1[[j]]<-lme4::lmer(ff, data=dat2, REML=FALSE)
    print(j)
  }
  
  
  ###  inferring the coefficients or connection strengths for the network from the fitted model1
  
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
  
  G <- qgraph::qgraph(E,fade=FALSE,
                layout="spring",
                labels=labs, 
                lty=ifelse(E[,3]>solid,1,5),
                edge.labels=F,
                edge.color=edge.color, 
                title=titlePlot)
  
  return(G)
  
  
}  # end function



