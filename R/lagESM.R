
#' Function to compute lagged n variables in ESM data and add it to data
#'
#' @param data data frame 
#' @param subjnr variable, indicating subjects
#' @param level2 variable, indicating level 2 variable, e.g. days. Can be left empty. 
#' @param level1 variable, indicating level 1 variable, e.g. beeps. Must be specified.
#' @param lagn number of lags
#' @param varnames vector with the names of the variables for which lagged versions are computed
#'
#' @return data frame with lagged variables added, with suffix L1, L2, L3, respectively
#' @export
#'
#' @examples
#' data("news")
#' vars <- c("Fearful","Hopeful","Anxious","Down","Inspiring","Irritated","Relaxed","Insecure")
#' res <- lagESM(data = news, subjnr="subjnr", level2= "daynr",
#'        level1 = "beepnr", lagn = 1, varnames = vars)

lagESM <- function(data, subjnr="subjnr", level2 = NULL, level1 = "beepnr", lagn=1, varnames = NULL) {
  
  
if (lagn > 3)  return("Number of lags in function lagESM should not exceed 3","\n") 
if (is.null(level1))  
  return(cat("Level 1 is not correctly specified in function lagESM","\n")) 
if (is.null(varnames))  
  return(cat("Argument varnames is not specified in function lagESM","\n")) 
  
  if(!subjnr %in% names(data)) 
    return(cat("Name of subjnr is not correctly specified in function lagESM","\n")) 
  if(!level1 %in% names(data)) 
    return(cat("Name of level1 is not correctly specified in function lagESM","\n")) 
  if(!is.null(level2)) {
     if(!level2 %in% names(data)) 
       return(cat("Name of level2 is not correctly specified in function lagESM","\n")) 
     if(level1 == level2)
       return(cat("Level1 must be different from level2 in function lagESM","\n")) 
  }
if (is.null(level2)) {
  level2 <- "one"
  data$one <- 1
}
  
## add additional beeps at the end of each day with missings

rest <- names(data)[!names(data) %in% c(subjnr,level2, level1, varnames)]
dat <- data[,c(subjnr,level2, level1, varnames, rest)]

l2 <- unique(dat[,level2])[!is.na(unique(dat[,level2]))]
vdaynr <- rep(sort(l2), length(unique(dat[,subjnr])));  
vsubjnr <- rep(unique(dat[,subjnr]), each=length(l2))

a <- data.frame(cbind(vsubjnr,vdaynr))

a2 <- NULL; a3 <- NULL
l1 <- unique(dat[,level1])[!is.na(unique(dat[,level1]))]

                a1 <- a;   a1[,level1] <- max(l1) + 1 
if (lagn > 1)  {a2 <- a;   a2[,level1] <- max(l1) + 2 }
if (lagn == 3) {a3 <- a;   a3[,level1] <- max(l1) + 3 }
                
               
a <- rbind(a1,a2,a3)  ;                            
a[,c(4:dim(dat)[2])] <- NA
names(a) <- names(dat)

b <- as.data.frame(rbind(dat, a))

b <- b[order(b[,subjnr],b[,level2],b[,level1]),]


## add lagged variables 

for (i in seq_along(varnames)) {   
  newname <- paste(varnames[i],"L",lagn, sep="")
  b <-  DataCombine::slide(b, Var = varnames[i], slideBy = -lagn, 
                           NewVar = newname,  reminder = FALSE)
}

## remove additional beeps

b <- b[b[,level1] <= max(l1),]

return(b)

}   # end function








