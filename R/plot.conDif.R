
#' Function to plot the results of the testCF function for network connectivity 
#' differences between two groups
#'
#' @param x result of conDif 
#' @param type number that indicates what output to plot, 
#'            1 = mean differences total
#'            2 = mean differences auto-regression, 
#'            3 = mean differences, excluding auto regression,
#'            4 = mean differences of SD's 
#' @export
#'
plot.conDif <- function(x, type=1) {
  
  a <- x$output$permutations[,type]
  est <- x$output$pvals[type,1]
  lab <- colnames(x$output$permutations)[type]
  
df <- with(stats::density(a), data.frame(x, y))
meanEst <- mean(df$x)
cutoff1 <- stats::quantile(a,0.025)
cutoff2 <- stats::quantile(a,0.975)

ylim1 <- colMeans((df[(min(abs(df$x - est)) == abs(df$x - est)),])[2])
ylim2 <- colMeans((df[(min(abs(df$x - meanEst)) == abs(df$x - meanEst)),])[2])

p <- ggplot(data = df, aes(x = df$x, y = df$y)) + geom_line()
p <- p + geom_segment(aes(x=est, y=0, xend = est, yend=ylim1),color="blue", linetype="dashed", size=.5)
p <- p + geom_segment(aes(x=meanEst, y=0, xend = meanEst, yend=ylim2),color="black", size=.5)
p <- p + geom_ribbon(data=subset(df ,x<=cutoff1 ),aes(ymax= df$y,ymin=0),
                     fill="gray10",colour=NA,alpha=0.5)
p <- p + geom_ribbon(data=subset(df ,x>=cutoff2 ),aes(ymax= df$y,ymin=0),
                     fill="gray10",colour=NA,alpha=0.5)
p <- p + ggtitle(paste("Permutation distribution",lab," with observed estimate"))

 return(p)

}
