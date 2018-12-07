
#' Function to compute lagged n variables in ESM data and add it to data
#'
#' @param dat data frame 
#' @param subjnr variable, indicating subjects
#' @param daynr variable, indicating days
#' @param beepnr variable, indicating beeps
#' @param lagn number of lags
#' @param varnames vector with the names of the variables for which lagged versions are computed
#'
#' @return data frame with lagged variables added, with suffix L1, L2, L3, respectively
#' @export
#'
#' @examples
#' data("newsData")
#' res <- lagESM(dat = newsData, subjnr="subjnr", daynr= "daynr",
#'        beepnr = "beepnr", lagn = 1, varnames = vars)

LagESM <- function(dat, subjnr="subjnr", daynr="daynr", beepnr="beepnr", lagn=1, varnames) {
 

if (lagn > 3) {print("number of lags should not exceed 3"); return() }
   
## add additional beeps at the end of each day with missings
 
  
vdaynr <- rep(sort(unique(dat[,daynr])), length(unique(dat[,subjnr])));  
vsubjnr <- rep(unique(dat[,subjnr]), each=length(unique(dat[,daynr])) )

a <- data.frame(cbind(vsubjnr,vdaynr))

a2 <- NULL; a3 <- NULL

                a1 <- a;   a1[,beepnr] <- max(unique(dat[,beepnr])) + 1 
if (lagn > 1)  {a2 <- a;   a2[,beepnr] <- max(unique(dat[,beepnr])) + 2 }
if (lagn == 3) {a3 <- a;   a3[,beepnr] <- max(unique(dat[,beepnr])) + 3 }
                
               
a <- rbind(a1,a2,a3)  ;                            
a[,c(4:dim(dat)[2])] <- NA
names(a) <- names(dat)

b <- as.data.frame(rbind(dat, a))

# s1 <- b[,subjnr]; s2 <- b[,daynr]; s3 <- b[,beepnr]
# b <- b[order(s1,s2,s3),]

b <- b[order(b[,subjnr],b[,daynr],b[,beepnr]),]


## add lagged variables 

L <- length(varnames)

for (i in 1:L) {   
  newname <- paste(varnames[i],"L",lagn, sep="")
  b[,newname] <- lag(b[,c(varnames[i])], n=lagn)
}


## remove additional beeps

b <- b[b[,beepnr] <= (max(unique(dat[,beepnr]))),]

return(b)

}   # end function








