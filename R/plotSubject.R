
#' Network plot of an individual 
#'
#' @param res result of class "esmNetwork" 
#' @param subjectID identification number that indicates the subject to be plotted
#' @param layout layout of the plot (see documentation of qgraph)
#' @param solid effect size above which lines are shown as solid (default = .10)
#' @param plimit p-value under which lines are shown (default = .05)
#'
#' @return an individual qgraph plot (network) 
#' @export
#'
#' @examples
#' plotSubject(res = res, subjectID = 31623)
#' 
plotSubject <- function(res, subjectID, layout = "circular", plimit = .05, solid = .10 ) {
 
  dat <- res$input$dat
  subjnr <- res$input$subjnr
  vars <- res$input$vars
  labs <- res$input$labs
  randomVars <- res$input$randomVars
  randomAll <- res$input$randomAll
  nvars <- res$intermediate$numberOfVars
  npred <- res$intermediate$numberOfPreds
  edge.color <- res$intermediate$edgeColor
  model <- res$output$model
  
  if (!subjectID %in% dat[,subjnr]) stop("subject identification is not in the data")
  if(is.null(randomVars) & !randomAll) warning("Only fixed (slope) effects in model: individual plots are superfluous")
  
  ncov <- npred - nvars
  nrand <- length(randomVars)
 
  nsubj <- length(unique(dat[,subjnr]))
  select <- c(1:nsubj)[subjectID == (unique(dat[,subjnr]))]

## The coefficients for the individual are computed
  
  ind <- array(0,c(nvars,nvars))
  
  if (randomAll) {
    for (j in 1:nvars){
        ind[,j] = as.numeric((lme4::fixef(model[[j]])[(2+ncov):(npred+1)] + 
                              lme4::ranef(model[[j]])[[1]][select,2:(nvars+1)]))
    }
    # add zero's in the random effects for the variables with no random effects
  } else {  
    test <- vars %in% randomVars
    if (nrand == 0) a <- lme4::ranef(model[[j]])[[1]][select,1]
    for (j in 1:nvars) {
      if (nrand > 0) a <- lme4::ranef(model[[j]])[[1]][select,2:(nrand+1)] 
        b <- rep(0,length(test))
        b[test] <- a
        ind[,j] = as.numeric((lme4::fixef(model[[j]])[(2+ncov):(npred+1)] + 
                              unlist(b)))
    }
  }
  
## The individual network is constructed

jj = rep(1:nvars,each=nvars)
jk = rep(1:nvars,nvars)
E3 = data.frame(from = jk, to = jj, weight = as.vector(ind))


G3 <- qgraph::qgraph(E3, 
                     layout = layout,
                     labels = labs,
                     lty = ifelse(E3[,3] > solid, 1, 5),
                     edge.labels = F,
                     title = paste0("Subject ", subjectID),
                     title.cex = .7)

plot(G3)

}

