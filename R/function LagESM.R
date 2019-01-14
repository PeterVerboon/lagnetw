
#' Function to compute lagged n variables in ESM data and add it to data
#'
#' @param dat data frame 
#' @param subjnr variable, indicating subjects
#' @param level2 variable, indicating level 2 variable, e.g. days
#' @param level1 variable, indicating level 1 variable, e.g. beeps
#' @param lagn number of lags
#' @param varnames vector with the names of the variables for which lagged versions are computed
#'
#' @return data frame with lagged variables added, with suffix L1, L2, L3, respectively
#' @export
#'
#' @examples
#' data("newsData")
#' res <- lagESM(dat = newsData, subjnr="subjnr", level2= "daynr",
#'        level1 = "beepnr", lagn = 1, varnames = vars)

LagESM <- function(dat, subjnr="subjnr", level2 = NULL, level1, lagn=1, varnames) {
 
if (lagn > 3) {print("number of lags should not exceed 3"); return() }
   
b <- dat[order(dat[,subjnr],dat[,level2],dat[,level1]),]

L <- length(varnames)

  for (i in 1:L) {   
      newname <- paste(varnames[i],"L",lagn, sep="")
      b <-  DataCombine::slide(b, Var = varnames[i], slideBy = -lagn, NewVar = newname, reminder = FALSE)
  }

return(b)

}  








