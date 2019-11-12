
#' Function to plot the results of the permDif function for network connectivity 
#' differences between two groups
#' 
#' There are five lot types. The type number indicates what figure is plotted, 
#'            1 = mean differences total
#'            2 = mean differences auto-regression, 
#'            3 = mean differences, excluding auto regression,
#'            4 = mean differences, subset of predictors
#'            5 = mean differences of SD's 
#'
#' @param x result of permDif 
#' @param ... additional parameters
#' @return Four or five (when there subsets) figures are build and returned to the plot window
#' @export
#'
plot.permDif <- function(x, ...) {
  
  plotall <- list()

  for (type in c(1:5)) { 
    
  a <- x$output$permutations[,type]
  est <- x$output$pvalues.summary[type,1]
  lab <- colnames(x$output$permutations)[type]
  
  if (!is.na(est))  {
  
df <- with(stats::density(a), data.frame(x, y))
meanEst <- mean(df$x)
cutoff1 <- stats::quantile(a,0.025)
cutoff2 <- stats::quantile(a,0.975)

ylim1 <- colMeans((df[(min(abs(df$x - est)) == abs(df$x - est)),])[2])
ylim2 <- colMeans((df[(min(abs(df$x - meanEst)) == abs(df$x - meanEst)),])[2])


shade1 <- rbind(c(df[1, "x"], 0), subset(df ,x<=cutoff1 ),c(cutoff1,0) )
shade2 <- rbind(c(cutoff2,0), subset(df ,x>=cutoff2 ), c(df[nrow(df), "x"], 0))

p <- ggplot(data = df, aes(x = df$x, y = df$y)) + geom_line()
p <- p + geom_polygon(data = shade1, aes(shade1$x,shade1$y), fill = "grey")
p <- p + geom_polygon(data = shade2, aes(shade2$x,shade2$y), fill = "grey")
p <- p + geom_segment(aes(x=meanEst, y=0, xend = meanEst, yend=ylim2),color="black", size=.1)
p <- p + geom_segment(aes(x=est,     y=0, xend = est,     yend=ylim1),color="blue", linetype="dashed", size=.2)
p <- p + theme_bw() + xlab("") + ylab("")
p <- p + ggtitle(paste("Permutation distribution",lab," with observed estimate"))
plotall[[type]] <- p

cat(paste("building plottype: ", type, "\n"))

print(plotall[[type]]) 

  } else {
    cat(paste("No plot available for", lab, "\n"))
  }
  

}


}
