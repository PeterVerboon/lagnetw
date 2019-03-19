
#' Function to compute lagged n variables in ESM data and add it to data
#'
#' @param dat data frame 
#' @param subjnr variable, indicating subjects
#' @param level2 variable, indicating level 2 variable, e.g. days. Can be left empty. 
#' @param level1 variable, indicating level 1 variable, e.g. beeps. Must be specified.
#' @param varnames vector with the names of the variables which are centered.
#' @param center vector of same length as varnames that indicates for each variable which centering must be,
#'               applied, either "person" or "grand_mean" or a number. If only "person" is specified (default), 
#'               all variables are person centered. If center is NULL no variables are centered.
#' @param addmeans logical indicating whether a vector of person means should be added to the data 
#'               (works only if person centering is requested for that variable).             
#'
#' @return data frame with centered variables, and if requested, person means with suffix "_means".
#' @export
#'
#' @examples
#' data("newsData")
#' res <- centerESM(data = DataNews, subjnr="subjnr",
#'        addmeans = TRUE, 
#'        varnames = c("Lonely", "Relaxed", "Anxious"), 
#'        center = c("person", "person", "person"))

centerESM <- function(data, subjnr="subjnr", addmeans = TRUE, varnames = NULL, center = NULL) {
  
if (is.null(varnames)) return(cat("Argument varnames is not specified in function centerESM","\n"))

if (is.null(center))  return(cat("Warning: no centering is requested","\n"))
if (length(center) == 1 && center == "person")  center <- rep("person", length(varnames))
if (length(center) == 1 && center == "grand_mean")  center <- rep("grand_mean", length(varnames))
   
if (!is.null(center))
   if (length(varnames) != length(center))  
     return(cat("Number of variables in varnames does not equal length of center","\n")) 

  if(!subjnr %in% names(data)) 
    return(cat("Name of subjnr is not correctly specified in function centerESM","\n")) 

dat1 <- data


# centering of variables

  for (i in seq_along(varnames)) 
  {
    xx <- varnames[i]
    if (is.numeric(dat1[,xx])) {
      
      if (center[i] == "person") {
        pmean <- ave(dat1[,xx], dat1[,subjnr], FUN = mean, na.rm=TRUE)
        dat1[,xx] <- dat1[,xx] - pmean
        if (addmeans == TRUE) dat1[,paste0(varnames[i],"_means")] <- pmean
      }
      if (center[i] == "grand_mean") dat1[,xx] <- dat1[,xx] - mean(dat1[,xx], na.rm=TRUE)
      if (is.numeric(center[i])) dat1[,xx] <- dat1[,xx] - center[i]
    }
  }

 return(dat1)

}   # end function


#TEST

datC <- centerESM(data = DataNews, subjnr="Participant", 
                         addmeans = TRUE, 
                         varnames = c("Lonely", "Relaxed", "Anxious"), 
                         center = "person")
                 









