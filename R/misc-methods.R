## ===========================================================================
## writeLines for data frames
## ---------------------------------------------------------------------------
## write method to create HTML output
setMethod("writeLines", signature("data.frame", "file", "missing","missing"),
          function(text, con){
              if(length(text))
              {  
                  modes <- sapply(text, mode)
                  for(i in which(modes == "numeric"))
                       text[,i] <- as.character(round(as.numeric(text[,i]), 3))
                  
                  th <- paste("<th>",
                              colnames(text), "</th>", sep="")
                  th[ncol(text)] <- gsub("<th", "<th class=\"right\"", th[ncol(text)])
                  th <- paste("<tr>\n",
                              paste(th, sep="", collapse="\n"), "\n</tr>\n", sep="",
                              collapse="\n")
                  td <- apply(text, 1, function(x) paste("<td>", x, "</td>", sep=""))
                  if(!is.matrix(td))
                      td <- as.matrix(td)
                  td[ncol(text), ] <-
                      gsub("<td", "<td class=\"right\"", td[ncol(text)],)
                  td <- paste(apply(td,2,function(x)
                                    paste("<tr>\n",
                                          paste(x, collapse="\n", sep=""), "\n</tr>", sep="",
                                          collapse="/n")), collapse="\n", sep="")
                  paste("<tr>\n", td, "</tr>\n", sep="",
                        collapse="/n")
                  writeLines(paste("<table class=\"QAParameter\" align=\"center\">\n",
                                   th, td, "\n</table>", sep=""), con)
              }
          })




## display details about qaProcessSummary and write HTML object
setMethod("show", signature("qaProcessSummary"),
          function(object)
      {
          cat("qaProcess summary for", ncol(object@summary)-1,
              "parameters and", nrow(object@summary),
              "common samples\n") 
      })


## write method to create HTML output
setMethod("writeLines", signature("qaProcessSummary", "file", "missing","missing"),
          function(text, con)
      {

          samples <- rownames(text@summary)
          channels <- colnames(text@summary)
          ordNames <- names(sort(rowSums(sapply(text@ranges, c)),decreasing=T))
          
          panels <- length(text@panels)
          ## we want a summary of all failed qaChecks across panels
          sumRanges <- text@ranges[[1]]
          if(panels>1)
              for(r in 2:panels)
                  sumRanges <- sumRanges + text@ranges[[r]]
          ## The bounding table: columns are panels, rows are samples
          writeLines("<table class=\"qaPanelBound\">\n<tr>\n<td>\n", con)
          writeLines(paste("<table class=\"qaPanelSummary\">\n<tr",
                           "class=\"even\">\n<th class=\"left\">\n</th>"), con)
          if(panels>1)
              writeLines("<th class=\"left\"><span>Summary</span></th>", con)
          writeLines(paste(paste("<th class=\"right\"><a href=\"index", 1:panels, ".html\">",
                                 text@pnams, "\n<br><span>",
                                 shortNames(names(text@panels), n=12), "</span></a>",
                                 "\n</th>", sep="", collapse="\n"), "\n</tr>", sep=""), con)
          ## iterate over rows (samples)
          for(s in seq_along(samples))
          {
              thisSamp <- samples[s]
              class <- ifelse((s %% 2)==0, "even", "odd") 
              ## first the summary over all pannels
              writeLines(paste("\n<tr class=\"", class, "\">\n<th class=\"left\">",
                               thisSamp, "</th>", sep=""), con)
              if(panels>1)
              {
                  writeLines("\n<td class=\"sum\">", con) 
             
                  writeLines(htmlBarplot(text@summary[s,][ordNames], sumRanges[ordNames]), con)
                  writeLines("</td>", con)
              }
              ## now iterate over each pannel
              for(p in seq_len(panels)){
                  writeLines("\n<td class=\"bars\">", con)
                  writeLines(htmlBarplot(text@panels[[p]][s,][ordNames], text@ranges[[p]][ordNames],
                                         p, match(thisSamp,
                                                  text@mapping[[p]][,"sample"]),
                                         FALSE), con)
                  writeLines("</td>", con)
              }
              writeLines("</tr>", con)
          }
          writeLines("</table>\n</td>\n</tr>\n</table>", con)
      })


## create a HTML code "barplot" for a single sample
htmlBarplot <- function(data, ranges, panel=NULL, frame=NULL,
                        rownames=TRUE, ignoreChannel="time")
{
    sel <- match(tolower(ignoreChannel), tolower(names(data)))
    if(!is.na(sel)){
        data <- data[-sel]
        ranges <- ranges[-sel]
    }
    ld <- length(data)
    maxRange <- max(ranges)
    width <- ifelse(rownames, maxRange*20 + max(nchar(names(data)))*10, 0)
    out <- paste("\n<table class=\"qaBar\" align=\"center\"",
                 ifelse(is.null(panel), "", paste("onclick=\"link2Panel(", panel,
                                                  ", ", frame, ")\"", sep="")), ">", sep="")
    rows <- character(ld)
    for(i in seq_len(ld))
    {
        failed <- ifelse(is.na(data[i]), 0, data[i])
        passed <- ranges[i] - failed
        empty <- maxRange - ranges[i]
        lclass <- if(rownames) "" else "link "
        bclass <- if(ranges[i]==0) rep("", maxRange) else
        if(maxRange==1 && ranges[i]==1) "single " else
        c("left ", rep("", maxRange-2), "right ")
        td <-  c(rep(sprintf("<td class=\"%%s%spassed\"></td>\n", lclass), passed),
                 rep(sprintf("<td class=\"%%s%sfailed\"></td>\n", lclass), failed),
                 rep("<td class=\"%sempty\"></td>\n", empty))
        rn <- ifelse(rownames, paste("<th>", names(data)[i], "</th>", collapse="",
                                    sep=""), "")
        rows[i] <- sprintf("<tr>\n%s%s</tr>", rn,
                           paste(mapply(sprintf, td, bclass), sep="", collapse=""))
    }
    out <- c(out, rows, "</table>")
    return(out)
}

## truncate character vectors to 'n' chars for pretty names plotting
shortNames <- function(x, n=13)
{
    for(i in seq_along(x))
    {
        if(nchar(x[i])>n)
            x[i] <- paste(sapply(x[i], substring, 1,n), "...", sep="")
        if(x[i]=="")
            x[i] <- "&nbsp;"
    }
    return(x)
}

## formats the output from qaProcess objects to a tab delimited format 
## suitable for writing to a text file

txtFormatQAObject<- function(qp){

	switch(qp@type,
		"Density" =,
		"ECDF" =,
		"KLD" =,
		{	
			vals <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@x
								})	
								unlist(tm)		
					})))
			passed <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@passed
								})	
								unlist(tm)		
					})))	
			tbl <- combineAltCols(qp,vals,passed)	
			thresh <- qp@frameProcesses[[1]]@details$absolute.value
			if(is.null(thresh)){
				alpha <- qp@frameProcesses[[1]]@details$alpha
				tbl <- cbind(tbl,"Alpha" = matrix( rep(alpha,nrow(vals)),ncol=1))
			}else{
				tbl <- cbind(tbl,"Threshold" = matrix( rep(thresh,nrow(vals)),ncol=1)) 			
			}
			sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
					   ))	   				
			colnames(sumry) <- "Summary"	
			tbl <- cbind(tbl, "Summary" = sumry)		
		},
		"2DStat" = {
		
			channels <- names(qp@frameProcesses[[1]]@frameAggregators)
			thresh <- qp@frameProcesses[[1]]@details$thresh
			tempVals <- sapply(channels,function(x){
					qp@frameProcesses[[1]]@details$stat[[x]]
			})
			alqLen <- length(qp@frameProcesses[[2]]@details$stat[[1]][[1]])
			chnl <- sapply(channels,function(x){
				paste(x, "_Alq", seq_len(alqLen),sep="")
			})
			
			tmpS <- list()
			vals <- sapply(channels,function(x){
					tmp <- data.frame(tempVals[,x],stringsAsFactors=FALSE,check.names=F)
					rownames(tmp) <- chnl[,x]
					tmpS[[x]] <<- t(tmp)				
			})
			vals <- tmpS
			
			passed <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@passed
								})	
								unlist(tm)		
					})))
			colnames(passed) <- channels
			Threshold <- matrix(rep(thresh,length(qp@frameProcesses)),ncol=1)
			colnames(Threshold) <- "Threshold"
			
			sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
					   ))	   				
			colnames(sumry) <- "Summary"	
			newVal <- matrix( rep(paste(qp@type,qp@name,sep="_"),nrow(passed)),ncol=1)
			colnames(newVal) <- "Process"

			sapply(channels,function(x){
				psd <- passed[,x,drop=F]
				colnames(psd) <- paste(colnames(psd),"Passed",sep="_")
				tmp <- data.frame(vals[[x]],psd,check.names=F)
				newVal <<- cbind(newVal,tmp)
			})
			
			tbl <- cbind(newVal,Threshold)
			
	},
	"BoundaryEvents"={
	
			channels <- names(qp@frameProcesses[[1]]@frameAggregators)
			sampNames <- names(qp@frameProcesses)
			thresh <- qp@frameProcesses[[1]]@details$thresh
			tempVals <- sapply(channels,function(x){
					qp@frameProcesses[[1]]@details$bPerc[[x]]
			})
			alqLen <- length(qp@frameProcesses[[1]]@details$bPerc[[1]][[1]])
			chnl <- sapply(channels,function(x){
				paste(x, "_Alq", seq_len(alqLen),sep="")
			})
			
			tmpS <- list()
			vals <- sapply(channels,function(x){
				   
					tmp <- data.frame(tempVals[,x],stringsAsFactors=FALSE,check.names=F)
					rownames(tmp) <- chnl[,x]
					colnames(tmp) <- rownames(tempVals)
					tmpS[[x]] <<- t(tmp)					
			})
			vals <- tmpS
		
			
			passed <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@passed
								})	
								unlist(tm)		
					})))
			colnames(passed) <- channels
			Threshold <- matrix(rep(thresh,length(qp@frameProcesses)),ncol=1)
			colnames(Threshold) <- "Threshold"
			
			sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
					   ))	   				
			colnames(sumry) <- "Summary"	
			newVal <- matrix( rep(paste(qp@type,qp@name,sep="_"),nrow(passed)),ncol=1)
			colnames(newVal) <- "Process"
			
			sapply(channels,function(x){
			  
				psd <- passed[,x,drop=F]
				colnames(psd) <- paste(colnames(psd),"Passed",sep="_")
				tmp <- data.frame(vals[[x]],psd,check.names=F)
				newVal <<- cbind(newVal,tmp)
			})
			tbl <- cbind(newVal,Threshold)			

			},
	"margin events" = {
			sampNames <- names(qp@frameProcesses)
			channels <- names(qp@frameProcesses[[1]]@frameAggregators)
			vals <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@x
								})	
								unlist(tm)		
					})))
			passed <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@passed
								})	
								unlist(tm)		
					})))	
			tbl <- combineAltCols(qp,vals,passed)	
			thresh <- qp@frameProcesses[[1]]@details$absolute.value
			if(is.null(thresh)){
					cFactor <- qp@frameProcesses[[1]]@details$cFactor
					m <- qp@frameProcesses[[1]]@details$m
					s <- qp@frameProcesses[[1]]@details$s
					
					if(!is.null(m ) || !is.null(s)){
						threshLow <- m - s*cFactor
						threshHigh <- m + s*cFactor
					}else{
						threshLow <- threshHigh <- rep(NA,length(sampNames))
					}
					threshLow <- matrix(rep(threshLow,length(sampNames)),ncol=length(channels),byrow=T)
					colnames(threshLow) <- paste("threshLow",channels,sep="_")
					threshHigh <- matrix(rep(threshHigh,length(sampNames)),ncol=length(channels),byrow=T)
					colnames(threshHigh) <- paste("threshHigh",channels,sep="_")
					thresh <- cbind(threshLow,threshHigh)
					
			
			}else{
					thresh <- matrix(rep(qp@frameProcesses[[1]]@details$absolute.value,length(sampNames)),ncol=1)		
					colnames(thresh) <- "Threshold"			
			}
			sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
								 ))	   
			colnames(sumry) <- "Summary"
			tbl <- cbind(tbl,thresh,sumry)	
	},
	"cell number"= {
			sampNames <- names(qp@frameProcesses)
 			absolute.value <- qp@frameProcesses[[1]]@details$absolute.value
            co <- qp@frameProcesses[[1]]@details$co  

	     	thresh <- matrix(rep(co,length(sampNames),ncol=1))
			colnames(thresh) <- "Threshold"
			rownames(thresh) <- sampNames
			vals <- t(data.frame(lapply(qp@frameProcesses,function(x){
									x@details$qaScore
							   })
							  ))
            if(is.null(absolute.value))
				colnames(vals) <- "CalcValues"
            else
				colnames(vals) <- "CellNumber"
			sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
								 ))	 
			colnames(sumry) <- "Summary"
			passed <- sumry
			#passed[,1] <- NA
			colnames(passed) <- "Passed"
			
			process<- matrix( rep(paste(qp@type,qp@name,sep="_"),nrow(vals)),ncol=1)
			colnames(process) <- "Process"
			tbl <- cbind(data.frame(process),vals,thresh,passed,sumry)
	},
	"time line" = { ##qaProcess.timeline  ## anything less than or equal to zero is pass
			channels <- names(qp@frameProcesses[[1]]@frameAggregators)
			vals <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
													y@x
								})	
								unlist(tm)		
					})))
			colnames(vals) <- channels
			passed <- t(data.frame(lapply(qp@frameProcesses,function(x){
									tm <-lapply(x@frameAggregators,function(y){
												y@passed
								})	
								unlist(tm)		
					})))	
			colnames(passed) <- channels
			tbl <- combineAltCols(qp,vals,passed)	
			
			sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
					   ))	   
			colnames(sumry) <- "Summary"
			tbl <- cbind(tbl,sumry)	
	},
	"time flow" = { 
	## qaProcess.timeflow 
	## NA is filled in for passed field 
	## anything less than or equal to zero is pass
			sampNames <- names(qp@frameProcesses)
			vals <- t(data.frame(lapply(qp@frameProcesses,function(x){
										x@details$qaScore
							   })
							  ))
			colnames(vals) <- "TimeFlow"

			thresh <- matrix(rep(0,length(sampNames),ncol=1))
			colnames(thresh) <- "Threshold"
			rownames(thresh) <- sampNames	
	       	sumry <- t(data.frame(lapply(qp@frameProcesses, function(x){
										x@summaryAggregator@passed
								  })
								 ))	 
			colnames(sumry) <- "Summary"
			
			passed <- sumry
			#passed[,1] <- NA
			colnames(passed) <- "Passed"
			process<- matrix( rep(paste(qp@type,qp@name,sep="_"),nrow(vals)),ncol=1)
			colnames(process) <- "Process"
			tbl <- cbind(data.frame(process),vals,passed,thresh,sumry)			
	}
	)
	tbl
	
}

### helper function used by txtFormatQAObject to combine results of a particular
### parameter together

combineAltCols <- function(qp,vals,passed){
	parms <- colnames(vals)
	if(!all(colnames(passed) %in% parms))
		stop("Column names dont match, combine failed")
		
	newVal <- matrix( rep(paste(qp@type,qp@name,sep="_"),nrow(vals)),ncol=1)
	colnames(newVal) <- "Process"

	sapply(parms,function(x){
		tmp <- data.frame(vals[,x,drop=F],passed[,x,drop=F],check.names=F)
		colnames(tmp) <- c(x,paste(x,"Passed",sep="_"))
		newVal <<- cbind(newVal,tmp)
	})
	newVal
}

