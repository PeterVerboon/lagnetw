
#' Function to plot the results of the testCF function for network connectivity 
#' differences between two groups
#'
#' @param x result of testCF 
#' @param type number that indicates what output to plot, 
#'            1 = mean differences total
#'            2 = mean differences auto-regression, 
#'            3 = mean differences, excluding auto regression,
#'            4 = mean differences of SD's 
#' @export
#'
plot.testCF <- function(x, type=1) {
  
  a <- x$output$perms[,type]
  outp <- x$output$pval[type,1]
  
df <- with(density(a), data.frame(x, y))
meanEst <- mean(df$x)
est <- outp[1,1]
cutoff1 <- quantile(a,0.025)
cutoff2 <- quantile(a,0.975)

ylim1 <- as.numeric((df[(min(abs(df$x - est)) == abs(df$x - est)),])[2])
ylim2 <- as.numeric((df[(min(abs(df$x - meanEst)) == abs(df$x - meanEst)),])[2])

p <- ggplot(data = df, aes(x = x, y = y)) + geom_line()
p <- p + geom_segment(aes(x=est, y=0, xend = est, yend=ylim1),color="blue", linetype="dashed", size=.5)
p <- p + geom_segment(aes(x=meanEst, y=0, xend = meanEst, yend=ylim2),color="black", size=.5)
p <- p + geom_ribbon(data=subset(df ,x<=cutoff1 ),aes(ymax=y,ymin=0),
                     fill="gray10",colour=NA,alpha=0.5)
p <- p + geom_ribbon(data=subset(df ,x>=cutoff2 ),aes(ymax=y,ymin=0),
                     fill="gray10",colour=NA,alpha=0.5)
p <- p + ggtitle("Permutation distribution with observed estimate")

 return(p)

}
