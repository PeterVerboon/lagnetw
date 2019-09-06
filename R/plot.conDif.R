
#' Function to plot the results of the conDif function for network connectivity 
#' differences between two groups
#' 
#' There are four lot types. The type number indicates what figure is plotted, 
#'            1 = mean differences total
#'            2 = mean differences auto-regression, 
#'            3 = mean differences, excluding auto regression,
#'            4 = mean differences of SD's 
#'
#' @param x result of conDif 
#' @return Four figures are build with ggplot and returned to the plot window
#' @export
#'
plot.conDif <- function(x,...) {
  
  plotall <- list()

  for (type in c(1:4)) { 
    
  a <- x$output$permutations[,type]
  est <- x$output$pvals[type,1]
  lab <- colnames(x$output$permutations)[type]
  
df <- with(stats::density(a), data.frame(x, y))
meanEst <- mean(df$x)
cutoff1 <- stats::quantile(a,0.025)
cutoff2 <- stats::quantile(a,0.975)

ylim1 <- colMeans((df[(min(abs(df$x - est)) == abs(df$x - est)),])[2])
ylim2 <- colMeans((df[(min(abs(df$x - meanEst)) == abs(df$x - meanEst)),])[2])


shade1 <- rbind(c(df[1, "x"], 0), subset(df ,x<=cutoff1 ),c(cutoff1,0) )
shade2 <- rbind(c(cutoff2,0), subset(df ,x>=cutoff2 ), c(df[nrow(df), "x"], 0))

p <- ggplot(data = df, aes(x = df$x, y = df$y)) + geom_line()
p <- p + geom_segment(aes(x=est, y=0, xend = est, yend=ylim1),color="blue", linetype="dashed", size=.2)
p <- p + geom_segment(aes(x=meanEst, y=0, xend = meanEst, yend=ylim2),color="black", size=.2)
p <- p + geom_polygon(data = shade1, aes(x, y), fill = "grey")
p <- p + geom_polygon(data = shade2, aes(x, y), fill = "grey")
p <- p + theme_bw() + xlab("") + ylab("")
p <- p + ggtitle(paste("Permutation distribution",lab," with observed estimate"))
plotall[[type]] <- p

cat(paste("building plottype: ", type, "\n"))

print(plotall[[type]]) 

}


}
