## re-evaluate a process for new thresholds or cutoff values 
evaluateProcess <- function(process, thresh, ...)
{
    switch(process@type,
           "time line"=
       {
           efun <- function(x, c,two.sided="missing", absolute.value="missing"){
               qaScore <- flowViz:::computeQAScore(x@details$raw, c)
               for(i in 1:length(x@frameAggregators)){
                   x@frameAggregators[[i]]@passed <- qaScore[i]==0
                   x@frameAggregators[[i]]@x <- qaScore[i] 
               }
               nfail <- !sapply(x@frameAggregators, slot, "passed")
               val <- if(sum(nfail)==1) factor(2) else factor(0)
               if(sum(nfail)==0)
                   val <- factor(1)
               x@summaryAggregator <- discreteAggregator(val)
               return(x)
           }
           process@frameProcesses <- lapply(process@frameProcesses, efun,
                                            thresh) 
       },
           "margin events"=
       {
           efun <- function(x, c,two.sided="missing", absolute.value="missing"){
               
               sums <- x@details$events
               m <- x@details$m
               s <- x@details$s
               for(i in 1:length(x@frameAggregators)){
                   x@frameAggregators[[i]]@passed <-
                       sums[i] <= m[i]+s[i]*c & sums[i] >= m[i]-s[i]*c
                   x@frameAggregators[[i]]@x <- sums[i]
               }
               nfail <- !sapply(x@frameAggregators, slot, "passed")
               val <- if(sum(nfail)==1) factor(2) else factor(0)
               if(sum(nfail)==0)
                   val <- factor(1)
               x@summaryAggregator <- discreteAggregator(val)
               return(x)
           }
           process@frameProcesses <- lapply(process@frameProcesses, efun,
                                            thresh)    
       },

           "time flow"=
       {
           efun <- function(x, c,two.sided="missing", absolute.value="missing"){
               qaScore <- flowViz:::computeQAScore(x@details$raw, c)
               x@summaryAggregator@passed <- qaScore==0
               return(x)
           }
           process@frameProcesses <- lapply(process@frameProcesses, efun,
                                            thresh)    
       },
           "cell number"=
       {
           efun <- function(x, c, two.sided, absolute.value){
               x@details$absolute.value <- if(!missing(absolute.value)) absolute.value
               else NULL
               if(is.null(x@details$absolute.value))
               {
                   m <- x@details$mean
                   s <- x@details$sd
                   val <- x@summaryAggregator@x
                   if(!missing(two.sided))
                       x@details$two.sided <- two.sided
                   summary <- if(!x@details$two.sided) m-val else abs(m-val)
                   x@summaryAggregator@passed <- summary < c*s 
               }
               else
               {
                    x@summaryAggregator@passed <-
                        x@summaryAggregator@x > x@details$absolute.value
               }
               return(x)
           }
           process@frameProcesses <- lapply(process@frameProcesses, efun,
                                            thresh, ...)    
       },
           
           stop("Don't know how to deal with process of type '",  process@type,
                "'", call.=FALSE)
           )
    return(process)
}
