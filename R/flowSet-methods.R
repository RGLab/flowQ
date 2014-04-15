setMethod("outliers", signature("flowSet"),
          function(x, test, stat, ...){
              
              if(length(x) < 4)
                stop("Not sufficient number of samples to perform the test")
              dat <- fsApply(x, each_col, stat, na.rm=TRUE)

              
              if(test == "Grubb") 
                  ans <- outlierG(x, stat, type = 11, opposite = FALSE, two.sided = FALSE)
              if(test == "Dixon")
                ans <- outlierD(x, stat, type = 11, opposite = FALSE, two.sided = FALSE)
              if(test == "X2")
                ans <- outlierX2(x, stat, variance, opposite = FALSE)

              return(ans)
          })



outlierD <- function(x, stat, type = 11, opposite = FALSE, two.sided = FALSE){

    if(length(x) < 4 || length(x) > 30) stop("Sample size must be in range 3-30")
    dat <- fsApply(x, each_col, stat, na.rm=TRUE)

    ans <- apply(dat, 2, function(y){
        dixon.test(y, type = type, opposite = opposite, two.sided = two.sided)}
                 )

    list(dat, ans)
}

outlierG <- function(x, stat, type = 11, opposite = FALSE, two.sided = FALSE){

    if(length(x) < 4) stop("Not sufficient number of samples to perform the test")
    dat <- fsApply(x, each_col, stat, na.rm=TRUE)

    ans <- apply(dat, 2, function(y)
                 grubbs.test(y, type = type, opposite = opposite, two.sided = two.sided))

    list(dat, ans)
}

outlierX2 <- function(x, stat, variance, opposite = FALSE){

    if(length(x) < 4) stop("Not sufficient number of samples to perform the test")
    dat <- fsApply(x, each_col, stat, na.rm=TRUE)

    ans <- apply(dat, 2, function(y)
                 chisq.out.test(y, variance = var(y), opposite = opposite))

    list(dat, ans)
}



##outlierG(GvHD, median)
##outlierD(GvHD, median)
##outlierX2(GvHD, median)
