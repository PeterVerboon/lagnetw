
#' This function is used to optimally display the networks with qgraph.
#' Works with either color and trans a vector of equal length,
#' or one of the two of length 1.
#'
#' @param color color specification for the arrows in the plot 
#' @param trans transparency: an integer between 0 and 255, 
#'              0 being fully transparent and 255 being fully visable
#'
#' @return a vector with specificatiosn used in qgraph
#' @export
#'
#' @examples
#' specs <- addTrans(color = "green", trans=255)
addTrans <- function(color, trans)
{
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) 
    stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(grDevices::col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}
