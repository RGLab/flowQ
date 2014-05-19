## Create quasi-random guids. This is only based on the time stamp,
## not on MAC address.
##This function has been updated to reflect changes in the format.hexmode function in R
# guid <- function()
    # as.character(format.hexmode(as.integer(Sys.time())/runif(1)*proc.time()["elapsed"]))

guid <- function(len=10){
       ltrs <- c(LETTERS,letters)
       paste(c(sample(ltrs,1),sample(c(ltrs,0:9),len-1,replace=TRUE)),collapse="")
}


.mySd <- function(x, na.rm=FALSE){
    return(if(is.matrix(x)) apply(x, 2, sd, na.rm=na.rm) else sd(x, na.rm=na.rm))
}


## Check for the class of object x and its length and cast error if wrong
checkClass <- function(x, class, length=NULL, verbose=FALSE,
                       mandatory=TRUE)
{
    if(mandatory && missing(x))
        stop("Argument '", substitute(x), "' missing with no default",
             call.=verbose)
    msg <- paste("'", substitute(x), "' must be object of class ",
                 paste("'", class, "'", sep="", collapse=" or "), sep="")
    fail <- !any(sapply(class, function(c, y) is(y, c), x))
    if(!is.null(length) && length(x) != length)
    {
        if(!is.null(x))
        {
            fail <- TRUE
            msg <- paste(msg, "of length", length)
        }
    }
    if(fail) stop(msg, call.=verbose) else invisible(NULL)     
}

locateParameter<-function(flowList, parm, flowSetIndx, flowFrameIndx)
{   
    if(length(parm)!=1)
            stop("Only one parameter is to be specified")
    if(!is.null(flowList[[flowSetIndx]][[flowFrameIndx]]))
    {
    	temp <- pData(parameters(flowList[[flowSetIndx]][[flowFrameIndx]]))
        mIndx <- which(parm ==  as.character(temp[,"desc"]))
     	if(length(mIndx)>1)
	        stop("Multiple channels in a flowFrame stained for the same cell type")
        if(length(mIndx)==0)
			return(NA)
		else 
			return( as.character(temp[,"name"])[mIndx])
    }else{
	    return(NA)
    }
}


locateDuplicatedParameters<-function(flowList)
{
    alqLen<- length(flowList)
    names <-data.frame()
    for( i in seq_len(alqLen))
    {   
	    if(!all(fsApply(flowList[[i]],nrow)==1)){
			temp <- pData(parameters(flowList[[i]][[1]]))	
			mIndx <- is.na(temp["desc"])
			temp["desc"][mIndx] <- temp["name"][mIndx]
			names<-rbind(names,temp[1:nrow(temp)-1,"desc",drop=FALSE])
		}
    }
    dupes<- names[duplicated(names[,1]),1] 
    dIndx <- !duplicated(dupes)
    dupes<-dupes[dIndx]
	if(length(dupes)==0)
	  message("No duplicate parameters found")
	else
    return(as.character(dupes))	
}



preProcessFlowList <- function(flowList){
	missingPats  <- lapply(flowList, sampleNames)
	totNames <- unique(unlist(missingPats))

	newList <-list()
	for ( i in 1: length(flowList)){

		toAdd <- totNames[!totNames %in% missingPats[[i]]]
		temp <- flowList[[i]]	
		cols <- colnames(temp)
		m <- matrix(0, ncol=length(cols))
		colnames(m) <- cols
		mframe <- flowFrame(m)
		parameters(mframe) <- parameters(temp[[1]])

		tempList  <- as(temp,"list")
		for ( j in toAdd){
  			tempList[[j]] <- mframe
		}	
		newList[[i]] <- as(tempList, "flowSet")
		rownames(pData(newList[[i]]))<- c(rownames(pData(temp)), toAdd)
	}
	return(newList)
}


normalizeSets <- function(flowList,dupes,peaks=NULL)
{
    patientID<-sampleNames(flowList[[1]])
    alqLen<-length(flowList)
    for (cellType  in dupes)
    {
        for(i in patientID)
        {
            inList<-list()   
            alqTable <- rep(FALSE,alqLen)
            frameIndx<-1
            for(j in seq_len(alqLen))
            {
                parm <- locateParameter(flowList,cellType,j,i)
                if(length(parm)!=0 && !is.na(parm))
                {
                    value<-flowList[[j]][[i]][,parm]
                    colnames(value) <- cellType
                    inList[[frameIndx]] <- value            
                    frameIndx <- frameIndx + 1
                    alqTable[j] <- TRUE
                } ## end if(length ...
            } ## end for(j ...
                  
            inSet <- flowSet(inList)
            norm1 <- normalization(normFunction=function(x, parameters, ...)
                                   warpSet(x, parameters,peakNr=peaks,...),
                                   parameters=as.character(cellType),
                                   arguments=list(grouping=NULL),
                                   normalizationId="Norm")
            nData <- normalize(inSet,norm1)
            frameIndx<-1
            for(m in which(alqTable==T))
            {
                parm <- locateParameter(flowList,cellType,m,i)
                exprs(flowList[[m]]@frames[[i]])[,parm] <-
                    exprs(nData[[frameIndx]])
                frameIndx <- frameIndx +1
            } ## end for(m ...
        } ## end for(i ...
    } ## end for(celltype ...
    return(flowList)
}



.mybw.nrd0 <- function (x) 
{
    if (length(x) < 2L) 
        stop("need at least 2 data points")
    hi <- .mySd(x)
    if (!(lo <- min(hi, IQR(x)/1.34))) 
        (lo <- hi) || (lo <- abs(x[1L])) || (lo <- 1)
    0.9 * lo * length(x)^(-0.2)
}



## QA process indicating too many events on the margins in comparison to
## average number of margin events for a particular channel. The 'grouping'
## argument can be used to compare within groups. 'cFactor' is the cutoff
## value, essentially this is a factor in the confidence intervall defined
## by the standard deviation of the average number of margin events for a
## particular channel.
qaProcess.marginevents <- function(set, channels=NULL,side="both", grouping=NULL, outdir,
                                   cFactor=2, absolute.value=NULL,
                                   name="margin events",
                                   sum.dimensions=NULL, det.dimensions=c(3,1),
                                   pdf=TRUE,
                                   ...)
{
    ## some sanity checking
    checkClass(set, "flowSet")
    checkClass(outdir, "character", 1)
    checkClass(cFactor, "numeric", 1)
    checkClass(absolute.value, c("numeric", "NULL"), 1)
    if(!is.null(grouping))
        if(!is.character(grouping) || ! grouping %in% colnames(pData(set)))
            stop("'grouping' must be a character indicating one of the ",
                 "phenotypic variables in 'set'")
    cn <- colnames(set)
    parms <- if(is.null(channels)) setdiff(cn, flowCore:::findTimeChannel(set)) else
    if(!all(channels %in% cn)) stop("Invalid channel(s)") else channels
	
    if(!is.null(sum.dimensions))
          checkClass(sum.dimensions, "numeric")
    checkClass(det.dimensions, "numeric")
    checkClass(name, "character", 1)
    checkClass(pdf, "logical", 1)
    frameIDs <- sampleNames(set)
    ls <- length(set)
    lp <- length(parms)
    
    ## count events on the margins using an expression filter
    ranges <- range(set[[1]], parms)
    perc <- matrix(ncol=ls, nrow=lp, dimnames=list(parms, frameIDs))
    cat("computing margin events...")
	perc <- t(sapply(parms, function(x) 
        { 	
			bf <- boundaryFilter(x,side=side)
			tol <- bf@tolerance
			fsApply(set,function(y){
				dat <- exprs(y)[,x]
				ranges <- range(y)[,x]
				if(length(dat) ==0){
					res <- NULL
				}else{
					res <- 	switch(bf@side, 
							both={(dat >  (ranges[1] + tol)) & (dat < (ranges[2] - tol))},
							upper = dat < (ranges[2] - tol), lower = dat >  (ranges[1] + tol)
							)
					res[is.na(res)] <- TRUE	
				}
				trueCount <- sum(res==T)
				count <- length(dat)
				q <- 1- trueCount/count									
			})
		}))		
	colnames(perc) <- frameIDs
    cat("\n")
    perc[,frameIDs[fsApply(set,nrow)==1]] <- NA
	
    ## create summary plot
    require("lattice")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    sfile <- file.path(idir, "summary.pdf")

    #pdf(file=sfile, width=sum.dimensions[1], height=sum.dimensions[2])
	if(is.null(sum.dimensions))
		pdf(file=sfile)
	else 
		pdf(file=sfile, width=sum.dimensions[1], height=sum.dimensions[2])
    col.regions=colorRampPalette(c("white",  "darkblue"))(256)
    print(levelplot(t(perc[1:lp,]*100), scales = list(x = list(rot = 90)),
                    xlab="", ylab="", main="% margin events",
                    col.regions=col.regions))
	dev.off()
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)

    ## deal with groups if there are any
    frameProcesses <- list()
	grps <- if(!is.null(grouping)) pData(set)[,grouping] else rep(1,ls)
    grpsi <- split(1:ls, grps, drop=TRUE)
    sset <- split(set, grps)
    ## create graphs and aggregators for each frame (and each channel)
    cat("creating frame plots...")
    for(i in 1:length(set)){
        fnames <- NULL
        ## this will hold the aggregators for all channels
        agTmp <- aggregatorList()
        thisGrp <-
            if(is.null(grouping)) "1" else as.character(pData(set)[i,grouping])
        mv <- sv <- NULL
        for(j in 1:length(parms))
        {
            ## the frame and parameter specific density plots
            tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
                                          gsub("\\..*$", "", parms[j]), ".pdf",
                                          sep=""))
			pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
			par(mar=c(1,0,1,0))
			set <- sset[[thisGrp]]
           	passed <- TRUE
			if(nrow(set[[i]]) >1){
				bw <- .mybw.nrd0(exprs(set[[i]][,parms[j]]))
				dens <-density(exprs(set[[i]][,parms[j]]), bw=bw)
				rng <- extendrange(range(set[[i]][,parms[j]])[c("min","max"),],f=0.05)
				plot(dens,xlim=rng, type="n", axes=FALSE, ann=FALSE)
	     		polygon(dens,col="gray", border="blue")
				percTmp <- perc[j,grps==as.numeric(thisGrp)]
				m <- mean(percTmp[!is.na(percTmp)])
				s <- .mySd(percTmp[!is.na(percTmp)])
				## test whether the particular frame and channel passes the check
				## and use a rangeAggregator to store that information
				if(!is.na(perc[j,i])){
					passed <- if(is.null(absolute.value)){
									perc[j,i] <= m+s*cFactor & perc[j,i] >= m-s*cFactor
							}else{
									perc[j,i] <= absolute.value
							}
			    }
				mv <- c(mv, m)
				sv <- c(sv,s)
			}else{
				plot(1,1,axes=F,ann=F,cex=0)
				text(1,1,"MISSING",col="red")
			}
			dev.off()
			fnames <- c(fnames, tfile)
            
            agTmp[[j]] <- new("rangeAggregator", passed=passed, x=perc[j,i], min=0, max=1)
            cat(".")
        }

        ## summarize the results for separate channels
        names(agTmp) <- parms
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1 && length(nfail)>2) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        ba <- new("discreteAggregator", x=val)
        ## bundle up the graphs for all channels for this particular frame
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir, width=150, pdf=pdf)
        fid <- frameIDs[i]
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                summaryAggregator=ba,
                                                frameAggregators=agTmp,
                                                frameGraphs=fGraphs,
                                                details=list(events=perc[,i],
                                                mevents=rowMeans(perc),
                                                m=mv, s=sv,
                                                absolute.value=absolute.value,
												cFactor=cFactor))
    }
    
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name,
               type="margin events", summaryGraph=sgraph,
               frameProcesses=frameProcesses))
}    

## QA process indicating strange patterns in signal intensity over time.
## This is done for all channels at once now, with drilldown to the separate
## channel results. The cutoff is the variance cutoff used directly in function
## 'timeLinePlot'
qaProcess.timeline <- function(set, channels=NULL, outdir, cutoff=1,
                               name="time line",
                               sum.dimensions=NULL, det.dimensions=c(3,2),
                               pdf=TRUE, ...)
{
    ## some sanity checking
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    ls <- length(set)
    if(is.null(channels)){
        parms <- setdiff(colnames(set[[1]]), c("time", "Time"))
	}else{
        if(!all(channels %in% colnames(set[[1]])))
            stop("Invalid channel(s)")
        parms <- channels
    }
    lp <- length(parms)
    if(!is.character(outdir) || length(outdir)!=1)
        stop("'outdir' must be a valid file path")
    if(!is.numeric(cutoff) || length(cutoff)!=1)
        stop("'cutoff' must be numeric scalar")
	if(!is.null(sum.dimensions))
          checkClass(sum.dimensions, "numeric")
	checkClass(det.dimensions, "numeric")
    
    
    ## create summary plots for each channel
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    sfiles <- NULL
    summary <- vector(lp, mode="list")
    for(j in seq_len(lp))
    {
        sfile <- file.path(idir, paste("summary_", j, ".pdf", sep=""))
		if(is.null(sum.dimensions))
			pdf(file=sfile)
		else
			pdf(file=sfile,width=sum.dimensions[1],height=sum.dimensions[2])
    
        binSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
        if( !all(fsApply(set,nrow) ==1) ){
			summary[[j]] <- timeLinePlot(set[fsApply(set,nrow) !=1], parms[j], binSize=binSize,
                                     varCut=cutoff)
		}else{
			plot(1,1,axes=F,ann=F,cex=0)
			text(1,1,"MISSING",col="red")
		}
	    summary[[j]][sampleNames(set)[fsApply(set,nrow)==1 ]] <- 0
	    dev.off()
        sfiles <- c(sfiles, sfile)
        cat(".")
    }  ## make s
	
	
    sampID <- names(summary[[1]])
    ##glue together the summary graphs and create a qaGraph object
    sfile <- paste(idir, "summary.pdf", sep="/")
	system(paste("convert ", " -density 240x240 +append ",
		  paste(sfiles, collapse=" ")," ",sfile, sep=""))
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=max(350,150*lp), pdf=pdf)

    ## create graphs and aggregators for each frame and channel
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    for(i in seq_len(ls))
    {
        fnames <- NULL
        agTmp <- aggregatorList()
        for(j in seq_len(lp))
        {
            tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
                                          gsub("\\..*$", "", parms[j]), ".pdf",
                                          sep=""))
            pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
			if(nrow(set[[i]]) >1){
				timeLinePlot(set[i], parms[j],
                         main=paste("score=", signif(summary[[j]][i], 4), sep=""),
                         cex.main=2, binSize=binSize, varCut=cutoff)
			    agTmp[[j]] <- new("numericAggregator", passed=summary[[j]][[sampID[i]]]<=0,
                              x=as.numeric(summary[[j]][[sampID[i]]]))
			}else{
				plot(1,1,axes=F,ann=F,cex=0)
				text(1,1,"MISSING",col="red")
				agTmp[[j]] <- new("numericAggregator", passed=TRUE,
                              x=as.numeric(NA))
			}
            dev.off()
            fnames <- c(fnames, tfile)
            
            cat(".")
        }

        ## wrap graphs and aggregators in objects
        names(agTmp) <- parms
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        ba <- new("discreteAggregator", x=val)
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
                               width=min(220, lp*120), pdf=pdf)
        fid <- frameIDs[i]
        dr <- lapply(seq_len(lp), function(x) attr(summary[[x]], "raw")[[sampID[i]]])
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                summaryAggregator=ba,
                                                frameAggregators=agTmp,
                                                frameGraphs=fGraphs,
                                                details=list(raw=dr))
    }
            
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name="time line",
                     type="time line", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))
}    

    


## Detect disturbances in the flow of cells over time
qaProcess.timeflow <- function(set, outdir, cutoff=2, name="time flow",
                               sum.dimensions=NULL, det.dimensions=c(3,2),
                               pdf=TRUE, ...)
{
    ## some sanity checking
    if(!is(set, "flowSet"))
        stop("'set' needs to be of class 'flowSet'")
    ls <- length(set)
    if(!is.character(outdir) || length(outdir)!=1)
        stop("'outdir' must be a valid file path")
    if(!is.numeric(cutoff) || length(cutoff)!=1)
        stop("'cutoff' must be numeric scalar")
	if(!is.null(sum.dimensions))
          checkClass(sum.dimensions, "numeric")
	checkClass(det.dimensions, "numeric")

    sampID <- sampleNames(set)
    
    ## create summary plot and its associated qaGraph object
 
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    binSize <- min(max(1, floor(median(fsApply(set, nrow)/100))), 500)
    sfile <- file.path(idir, "summary.pdf")
	if(is.null(sum.dimensions))
		pdf(file=sfile)
	else
		pdf(file=sfile, width=sum.dimensions[1], height=sum.dimensions[2])
	if( !all(fsApply(set,nrow) ==1) ){
		summary <- timeLinePlot(set[fsApply(set,nrow) !=1], colnames(set)[[1]], binSize=binSize,
                            varCut=cutoff, type="frequency")
	}else{
		plot(1,1,axes=F,ann=F,cex=0)
		text(1,1,"MISSING",col="red")
	}
    dev.off()
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)

    ## create graphs and aggregators for each frame and wrap in object
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
	
    for(i in seq_len(ls))
    { 	agTmp <- aggregatorList()
		#fnames <- NULL
	    tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), ".pdf",
                                      sep=""))
        pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
        if(nrow(set[[i]]) >1){
			sum <- timeLinePlot(set[i], colnames(set)[[1]],
						 main=paste("score=", signif(summary[i], 4), sep=""),
						 cex.main=2, binSize=binSize, type="frequency",
						 varCut=cutoff)
			ba <- new("binaryAggregator", passed= sum<=0)
			
		}else{
			plot(1,1,axes=F,ann=F,cex=0)
			text(1,1,"MISSING",col="red")
			ba <- new("binaryAggregator", passed=TRUE)
			sum <- NA
			
		}
		
        dev.off()
		agTmp[[1]] <- ba
	    fGraphs <- qaGraphList(imageFiles=tfile, imageDir=idir,
                               width=min(220, 120), pdf=pdf)
		
		#fg <- qaGraph(fileName=tfile, imageDir=idir, width=220, pdf=pdf)
        fid <- frameIDs[i]
		#frameProcesses[[fid]] <-
        #    qaProcessFrame(fid,ba,fg, details=list(qaScore=sum))
	
		nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        tm <- new("discreteAggregator", x=val)
		 frameProcesses[[fid]] <- 
				qaProcessFrame(frameID=fid, summaryAggregator=tm,
                                  frameAggregators=agTmp,
							   frameGraphs=fGraphs,
                               details=list(qaScore=sum))

        cat(".")
    }
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name, type="time flow",
                     summaryGraph=sgraph, frameProcesses=frameProcesses))
}    


createImageDir <- function(outdir, gid)
{
    id <- file.path(outdir, "images", gid)
    if(!file.exists(id))
        dir.create(id, recursive=TRUE)
    return(win2UnixPath(id)) 
}


## Detect unusually low cell counts
qaProcess.cellnumber <- function(set, grouping=NULL, outdir, cFactor=2,
                                 absolute.value=NULL, two.sided=FALSE,
                                 name="cell number", sum.dimensions=NULL,
                                 pdf=TRUE, ...)
{
    ## some sanity checking
    checkClass(set, "flowSet")
    checkClass(outdir, "character", 1)
    checkClass(cFactor, "numeric", 1)
    checkClass(absolute.value, c("numeric", "NULL"), 1)
    if(!is.null(grouping))
        if(!is.character(grouping) || ! grouping %in% colnames(pData(set)))
            stop("'grouping' must be a character indicating one of the ",
                 "phenotypic variables in 'set'")
	if(!is.null(sum.dimensions))
	    checkClass(sum.dimensions, "numeric")
    checkClass(two.sided, "logical", 1)
    checkClass(name, "character", 1)
    checkClass(pdf, "logical", 1)

    ## deal with groups if there are any
    ls <- length(set)
    grps <- if(!is.null(grouping)) pData(set)[,grouping] else rep(1,ls)
    grpsi <- split(1:ls, grps, drop=TRUE)
    sset <- split(set, grps)
    
    ## create summary plot and its associated qaGraph object
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    cellNumbers <- as.numeric(fsApply(set, nrow))
	cellNumbers[cellNumbers ==1] <- NA
    sfile <- file.path(idir, "summary.pdf")
	if(is.null(sum.dimensions))
		pdf(file=sfile)
	else
		pdf(file=sfile, width=sum.dimensions[1], height=sum.dimensions[2])
    col <- "gray"
    par(mar=c(10.1, 4.1, 4.1, 2.1), las=2)
	if( !all(is.na(cellNumbers))){
		barplot(cellNumbers, col=col, border=NA, names.arg=sampleNames(set),
				cex.names=0.8, cex.axis=0.8)
		abline(h=mean(cellNumbers), lty=3, lwd=2)
	}else{
		plot(1,1,axes=F,ann=F,cex=0)
		text(1,1,"MISSING",col="red")
	}
    dev.off()
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350, pdf=pdf)

    ## create aggregators for each frame and wrap in object
    frameIDs <- sampleNames(set)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    for(i in seq_len(ls))
    {
        thisGrp <-
            if(is.null(grouping)) "1" else as.character(pData(set)[i,grouping])
		cellNum <- cellNumbers[grpsi[[thisGrp]]]
		cellNum <- cellNum[!is.na(cellNum)]
        var <- .mySd(cellNum)
        m <- mean(cellNum)
        if(is.null(absolute.value))
        {
            summary <- if(!two.sided) m-cellNumbers[i] else abs(m-cellNumbers[i])
            co <- var*cFactor
            ba <- new("numericAggregator", x=cellNumbers[i],
                      passed=if(!is.na(summary)) summary<co else TRUE )
        }else{   
            summary <- cellNumbers[i]
            co <- absolute.value
            ba <- new("numericAggregator", x=cellNumbers[i],
                      passed=if(!is.na(summary)) summary > co else TRUE)
        }
        fid <- frameIDs[i]
        frameProcesses[[fid]] <-
            qaProcessFrame(fid, ba, details=list(qaScore=summary, mean=m, sd=var,co=co,
                                    two.sided=two.sided, absolute.value=absolute.value,
									cFactor=cFactor))
        cat(".")
    }## end for(i in ...

    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name, type="cell number",
                     summaryGraph=sgraph, frameProcesses=frameProcesses))
}   




## Similar to qaProcess.marginevents but this will compare boundary events across
## panels.
qaProcess.BoundaryPlot <- function(flowList, dyes=NULL, outdir="QAReport",
                                   cutoff=3,sum.dimensions=NULL, det.dimensions=NULL,
                                   pdf=TRUE,name="Boundary",side= "both",...)
{
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    sfiles <- NULL
    alqLen<- length(flowList)
    Pats  <- lapply(flowList, sampleNames)
    patientID <- unique(unlist(Pats))
    
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    parLbl<-vector(mode="character",length=alqLen)
    legend <-vector(mode="character",length=alqLen)
    if(is.null(dyes)){
	dupes <- locateDuplicatedParameters(flowList)
    if(length(dupes)==0)
			stop("Duplicated parameters do not appear in the 
                  list of flowSets provided")
    }else{
	dupes <- as.character(dyes)
    }
    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
    panelFlag <-list()
	boundPerc <- list()

    for (cellType in dupes ){
        res<-data.frame()
        panelFlag[[cellType]] <-list()
		boundPerc[[cellType]] <- list()
		for( i in patientID){
            perc<-matrix(nrow=alqLen,ncol=1)
            panelFlag[[cellType]][[i]] <- TRUE
            for(j in seq_len(alqLen)){
				par <- locateParameter(flowList,cellType,j,i)
				if(length(par)!=0 && !is.na(par) && nrow(flowList[[j]]@frames[[i]])!=1 ){
					legend[j] <-"green"
					parLbl[j] <- paste(j," ",par)   #####
					bf <- boundaryFilter(par,side=side)
					tol <- bf@tolerance
					dat <- exprs(flowList[[j]]@frames[[i]])[,par]
					ranges <- range(flowList[[j]]@frames[[i]])[,par]
					if(length(dat) ==0){
						bound <- NULL
					}else{
						bound <- 	switch(bf@side, 
								both={(dat >  (ranges[1] + tol)) & (dat < (ranges[2] - tol))},
								upper = dat < (ranges[2] - tol), lower = dat >  (ranges[1] + tol)
								)
						bound[is.na(bound)] <- TRUE	
					} ## end of dat==0
					trueCount <- sum(bound==T)
					count <- length(dat)
					perc[j,]  <- (1- trueCount/count)*100									
				}else{
					legend[j] <-"white"
					parLbl[j] <- paste(j," ")
				} ## end of length(par)!=0        
  
            }  ## end of alqLen
		
            colnames(perc)<-cellType
            #legend=rep("green",length(perc))
            legend[perc>cutoff]<-"red"
            panelFlag[[cellType]][[i]] <-all(perc[!is.na(perc)]<=cutoff)
			boundPerc[[cellType]][[i]] <- perc
			passed <- rep(TRUE,alqLen)
			passed[perc > cutoff] <- FALSE
            newres<-data.frame(Patient=rep(i,alqLen),Aliquot=seq_len(alqLen),
                              passed=factor(c(passed),levels=c(TRUE,FALSE)),
                              data=perc,check.names=F)
            res<-rbind(res,newres)     
            formula <- paste("`","Aliquot","`","~","`",cellType,"`",sep="")      
            tempgrph[[cellType]][[i]]<-
                         barchart( eval(parse(text=formula)), 
                                data=newres,origin = 0,
                                col=myCol[unique(newres[,"Aliquot"])],
                                key=list(space="right",points=list(col=legend),
                                text=list(parLbl),col=myCol))                                                           
            cat(".")
		} ## end of patientID
     
    	sfile <- file.path(idir, paste("summary_", cellType, ".pdf", sep=""))
		if(is.null(sum.dimensions))
			pdf(file=sfile)
		else
			pdf(file=sfile, width=sum.dimensions[1],height=sum.dimensions[2])
        formula <- paste("`","Aliquot","`","~","`",cellType,"`","|","Patient",sep="")   
        print(barchart( eval(parse(text=formula)), 
                data=res,col=(c("green","red")),origin = 0,drop.unused.levels=F,
		main="Percentage of margin events",
                groups=res[,"passed"],
 		key=simpleKey(text=as.character(c("failed","passed")),
			space="right",points=FALSE,col=c("red","green")))
        )
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")	
    } ## end of dupes
    sfile <- paste(idir, "summary.pdf", sep="/")
    # system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 # lp, "x1 ", sfile, sep=""))

	system(paste("convert ", " -density 240x240 +append ",
		  paste(sfiles, collapse=" ")," ",sfile, sep=""))
	
	
    sgraph <- qaGraph(fileName=sfile, imageDir=idir,width=max(350,200*lp),pdf=pdf)
	frameProcesses <- list()
    cat("\nCreating frame plots...")

    for(i in seq_len(ls)) ##over patient
    {  
		fnames <- NULL
        agTmp <- aggregatorList()
        for(j in seq_len(lp)){  ##over nrow(dyes)
			tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
                                          gsub("\\..*$", "", j), ".pdf", sep=""))
			if(is.null(det.dimensions))
				pdf(file=tfile)
			else
				pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
			print(tempgrph[[j]][[patientID[i]]])
			dev.off()
			fnames <- c(fnames, tfile)
			agTmp[[j]] <- new("binaryAggregator", passed=panelFlag[[j]][[patientID[i]]])
			#   agTmp[[j]] <- new("discreteAggregator", passed=panelFlag[[j]][patientID[i]],x=1)
			##names(agTmp[[j]]) <-paste(dyes[j])
                        names(agTmp)[j] <-paste(dyes[j])
			cat(".")
		}##end of lp
        names(agTmp) <-dupes
		nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail[!is.na(nfail)])==1) factor(2) else factor(0)
     	if(sum(nfail[!is.na(nfail)])==0)
             val <- factor(1)
        ba <- new("discreteAggregator", x=val,passed= as.logical(sum(nfail==0)))
		fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
				  width=200, pdf=pdf)
		fid <- patientID[i]
		frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
						    summaryAggregator=ba,
						    frameAggregators=agTmp,
						    frameGraphs=fGraphs,details=list(bPerc = boundPerc,thresh=cutoff))
        cat(".")
    } ## end of ls
    cat("\n")
	output<-qaProcess(id=gid, name=name, type="BoundaryEvents",	summaryGraph=sgraph, 
						frameProcesses=frameProcesses)
	return(output)
}

qaProcess.2DStatsPlot <- function(
		flowList,
		dyes=c("FSC-A","SSC-A"),
		outdir="QAReport",
		thresh=0.25, # numeric between 0/1 indicating outlier boundary defualt 0.25
		func=mean,	
		sum.dimensions=NULL,
		det.dimensions=NULL,
		pdf=TRUE,
		name ="2DStats",
		...)
{	
	if(is(dyes,"character")){
		if(length(dyes)!=2)
			dyes<- matrix(dyes,byrow=TRUE,ncol=2)
		else
			dyes <- matrix(dyes,ncol=2)
	}
	cat("creating summary plots...")
	gid <- guid()
	idir <- createImageDir(outdir, gid)
	sfiles <- NULL
	alqLen<- length(flowList)
	Pats  <- lapply(flowList, sampleNames)
	patientID <- unique(unlist(Pats))
	ls <- length(patientID)
	lp <- nrow(dyes)
	myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
	tempgrph<-list()
	parLbl<-vector(mode="character",length=alqLen)
	panelFlag<-list()
	pcoutVals <- list()
	
	for(cellType in seq_len(nrow(dyes))){
        tempgrph[[cellType]] <- list()
        panelFlag[[cellType]]<-list()
		pcoutVals[[paste(dyes[cellType,1],"/",dyes[cellType,2],sep="")]] <- list()
        outRes<-data.frame()
        for( i in patientID){
            res<-data.frame()
            panelFlag[[cellType]][i]<-TRUE            
            for(j in seq_len(alqLen)){
                par1 <- locateParameter(flowList,dyes[cellType,1],j,i)
                par2 <- locateParameter(flowList,dyes[cellType,2],j,i)
                if(length(par1)!=0 && length(par2)!=0 && nrow(flowList[[j]]@frames[[i]])!=1 && !is.na(par1) && !is.na(par2)) {
					parLbl[j] <- paste(j," ",par1,"/",par2)
					eps<- c(.Machine$double.eps, - .Machine$double.eps)
					fFrame <- flowList[[j]]@frames[[i]]
					valRange1<-range(fFrame[,par1])
					ranges <- t(valRange1) +eps
					ef <- char2ExpressionFilter(
						   paste("`", par1, "`>", ranges[1]," & `",par1,"`<",
						   ranges[2], sep="",
						   collapse=""), filterId=as.character(cellType))
					ff <-filter(fFrame,ef)
					fFrame <- Subset(fFrame,ff)
					valRange2<-range(fFrame[,par2])
					ranges <- t(valRange2) +eps
					ef <- char2ExpressionFilter(
						   paste("`", par2, "`>", ranges[1]," & `",par2,"`<",
						   ranges[2], sep="",
						   collapse=""), filterId=as.character(cellType))
					ff <-filter(fFrame,ef)
					fFrame <- Subset(fFrame,ff)
					value1 <- func(exprs(fFrame[,par1]))
					value2 <- func(exprs(fFrame[,par2]))
					newres<-data.frame(Aliquot=j,x=value1,y=value2, Patient=i,check.names=FALSE)
					res<-rbind(res,newres)
                }else{
                    parLbl[j] <- paste(j," ") 
                }
			} ### end of alqLen
			if(nrow(res)<1)
			     res<- data.frame(Aliquot=NA,x=NA,y=NA, Patient=i,check.names=F)
			colnames(res) <-c("Aliquot",paste(dyes[cellType,1]),paste(dyes[cellType,2]),"Patient")
            tempIndx <- NA
			tp <- rep(NA,alqLen)
			pcoutVals[[paste(dyes[cellType,1],"/",dyes[cellType,2],sep="")]][[i]] <- NA
			if(nrow(res) >1){
				tempIndx <- which(pcout(x=as.matrix(res[,2:3]))$wfinal<= thresh)
				outLier <- rep(FALSE,nrow(res))
				if(length(tempIndx)>0)
					panelFlag[[cellType]][[i]]<-FALSE
				tp[res[,"Aliquot"]] <- pcout(x=as.matrix(res[,2:3]))$wfinal
				pcoutVals[[paste(dyes[cellType,1],"/",dyes[cellType,2],sep="")]][[i]] <- tp
			}else{
			    pcoutVals[[paste(dyes[cellType,1],"/",dyes[cellType,2],sep="")]][[i]] <-tp
			}
			legend <-rep("white",alqLen)
			legend[res[,"Aliquot"]] <-"green"
			if(!is.na(tempIndx) && length(tempIndx)>0){
				legend[res["Aliquot"][tempIndx,]] <- "red" 
			}
		    
			### generate a graph for each patient
			formula<-paste("`",dyes[cellType,1],"`"," ", "~"," ","`",dyes[cellType,2],"`","|","Patient",sep="")
			tempgrph[[cellType]][[i]] <- xyplot(eval(parse(text=formula)),data=res,
                            groups=Aliquot,
                            ylim = if(nrow(res) >1) extendrange(res[,dyes[cellType,1]],f=3) else res[,dyes[cellType,1]],
                            xlim = if(nrow(res) >1) extendrange(res[,dyes[cellType,2]],f=3) else res[,dyes[cellType,2]],
                            col=myCol[unique(res[,"Aliquot"])],
                            key=list(space="right",points=list(pch=19,col=legend),
							text=list(parLbl),col=myCol),pch=19,cex=2
							)
			### create the summary figure with outlier/non outliers labelled
			if(nrow(res) >0){
				res["Aliquot"] <- "Non-Outlier"
				if(length(tempIndx) >0 && !is.na(tempIndx))
					res["Aliquot"][tempIndx,] <- "Outlier"
			}
			outRes <- rbind(outRes,res)
		    
		}  ### end of patientID
		cat(".")
		sfile <- file.path(idir, paste("summary_", paste(dyes[cellType,1],"_",dyes[cellType,2],sep=""), ".pdf", sep=""))
		if(is.null(sum.dimensions))
			pdf(file=sfile)
		else
			pdf(file=sfile, width=det.dimensions[1],height=det.dimensions[2])
		print(xyplot(eval(parse(text=formula)),data=outRes, 
						auto.key=list(space="right"), groups=Aliquot,
						ylim= if(nrow(outRes) >1) extendrange(outRes[,dyes[cellType,1]],f=1.8) else outRes[,dyes[cellType,1]],
						xlim= if(nrow(outRes) >1) extendrange(outRes[,dyes[cellType,2]],f=1.8) else outRes[,dyes[cellType,2]],
						col=c("green","red"),
						key=simpleKey(text=as.character(c("OutLier","Non-outlier")),
							space="right",points=F,col=c("red","green")),pch=19
							))
		dev.off()
		sfiles <- c(sfiles, sfile)
	} ### end of cellType
	sfile <- paste(idir, "summary.pdf", sep="/")
	system(paste("convert ", " -density 240x240  +append ",
		  paste(sfiles, collapse=" ")," ",sfile, sep=""))
	# system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                                    # lp, "x1 ", sfile, sep=""))
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=350,pdf=pdf)
    frameProcesses <- list()
    cat("\nCreating frame plots...")
	for(i in seq_len(ls)) ##over patient
    {  
        fnames <- NULL
        agTmp <- aggregatorList()
        for(j in seq_len(lp)){  ##over nrow(dyes)
            tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
                                            gsub("\\..*$", "", j), ".pdf",
                                            sep=""))
            pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
            print(tempgrph[[j]][[patientID[i]]])
            dev.off()
			if(!is.null(tempgrph[[j]][[patientID[i]]])){
            fnames <- c(fnames, tfile)}
            agTmp[[j]] <- new("binaryAggregator", passed=panelFlag[[j]][[patientID[i]]])
            names(agTmp)[j] <-paste(dyes[j,1],"/",dyes[j,2],sep="")
            cat(".")
        }
        dyeNames<-apply(dyes,1,function(x){ paste("",x[1],"/",x[2],"",sep="")
                        })
        names(agTmp) <-dyeNames
        nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
        if(sum(nfail)==0)
            val <- factor(1)
        ba <- new("discreteAggregator", x=val, passed =as.logical(sum(nfail==0)))
        fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
                        width=200, pdf=pdf)
        fid <- patientID[i]
        frameProcesses[[fid]] <- qaProcessFrame(frameID=fid, summaryAggregator=ba,
									frameAggregators=agTmp, frameGraphs=fGraphs,
									details =list(thresh=thresh, stat=pcoutVals))
        cat(".")
    } 
    cat("\n")	
	qaProcess(id=gid, name=name, type="2DStat", summaryGraph=sgraph, 
					  frameProcesses=frameProcesses)    
}

qaProcess.DensityPlot <- function(
	flowList, dyes=NULL, outdir="QAReport",	alpha=0.05,	absolute.value=NULL,
	sum.dimensions=NULL, det.dimensions=NULL, pdf=TRUE,	name="Density",...)
{
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    sfiles <- NULL
    alqLen<- length(flowList)
    Pats  <- lapply(flowList, sampleNames)
    patientID <- unique(unlist(Pats))
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    if(is.null(dyes)){
		dupes <- locateDuplicatedParameters(flowList)
		if(length(dupes)==0)
			stop("Duplicated parameters do not appear in the 
                  list of flowSets provided")
    }else{
		dupes <- as.character(dyes)
    }

    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
    valRange <- data.frame()

    for (cellType in dupes ){
        formula<-paste("~","`",cellType,"`","|","Patient",sep="")
        tempgrph[[cellType]]<-list()
        tempDist[[cellType]]<-list()
        tempInput <-list()
        parLbl<-vector(mode="character",length=alqLen)
        ymax<-0
        for( i in patientID){
            res<-data.frame()
            tempStat<-matrix(ncol=512,nrow=alqLen)
            for(j in seq_len(alqLen)){
                par <- locateParameter(flowList,cellType,j,i)
                if(length(par)!=0 && !is.na(par) && nrow(flowList[[j]]@frames[[i]])!=1){
					parLbl[j] <- paste(j," ",par)
                    eps<- c(.Machine$double.eps, - .Machine$double.eps)
					valRange<-range(flowList[[j]]@frames[[i]][,par])
                    ranges <- t(valRange) +eps                                 
                    ef <- char2ExpressionFilter(
								paste("`", par, "`>", ranges[1]," & `",par,"`<",ranges[2], 
								sep="",collapse=""), filterId=cellType)
					ff <-filter(flowList[[j]]@frames[[i]][,par],ef)
                    value=exprs(Subset(flowList[[j]]@frames[[i]][,par],ff))
                    colnames(value) <- cellType
            		dens <- density(value,n=512,from=valRange[1,],to=valRange[2,])
					newres <- data.frame(Patient=rep(i,512),Aliquot=rep(j,512),
										 dx=dens$x,dy=dens$y,check.names=FALSE)
				    res<-rbind(res,newres)
					tempStat[j,] <- dens$y
                    tempInput[[j]] <-value
                    yrng<-range(tempStat[j,])[2]
                    if(yrng>ymax)
                        ymax<-yrng
				}else{
                    parLbl[j] <- paste(j," ") 
                    tempInput[[j]] <-NA
                }
			} ###end of alqLen	
		
			if( !all(is.na(unlist(tempInput))) && 
				length(which(unlist(lapply(tempInput,length)) > 1))>1 ){
			    ## divide by the maximum value in range to get inputs in the range of 0 to 1
				## before calculating KL distance
			    dst <- KLdist.matrix(lapply(tempInput,"/",valRange[2,]),symmetrize=TRUE)
				## normalize by total count 
                tempDist[[cellType]][[i]] <- sum(dst,na.rm=T)/length(which(!is.na(dst)==T))
				tempgrph[[cellType]][[i]] <- 
						xyplot(dy ~ dx | Patient,data=res,groups=Aliquot,
								col=myCol[unique(res[,"Aliquot"])], ,ylab="Density",
								xlab=as.character(cellType),
								key=simpleKey(text=parLbl,space="right", points=F,col=myCol),
								plot.points=F,pch = ".", cex = 2, type = c("l"))
			}else{
				tempDist[[cellType]][[i]] <- NA
				m<-data.frame(x=1,y=1,z=1)
                tempgrph[[cellType]][[i]] <- densityplot(~x,data=m,main= "MISSING",xlab="")
			}
            cat(".")

        }### end of patientID
  
        xrange<-c(0.9*valRange[1,],1.1*valRange[2,])
        sfile <- file.path(idir, paste("summary_", cellType,".pdf", sep=""))
		if(is.null(sum.dimensions))
			pdf(file=sfile)
		else
			pdf(file=sfile, width=sum.dimensions[1],height=sum.dimensions[2])
	
		print(densityplot(~x|patientID, 
              data = list(patientID = 
		      factor(names(tempgrph[[cellType]]),
              levels = names(tempgrph[[cellType]])),
		      x = seq_along(tempgrph[[cellType]])), 
              xlim =xrange,ylim=c(0,1.05*ymax),
              xlab=as.character(cellType),
              key=simpleKey(text=parLbl,space="right",
		      points=F,col=myCol),
              lwd=2,
              panel = function(x, y, plot.list) {
                  do.call(panel.xyplot, 
			  trellis.panelArgs(
				  tempgrph[[cellType]][[x]], 1))
                  }
              ))
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
   }
  
    sfile <- paste(idir, "summary.pdf", sep="/")
	system(paste("convert ", " -density 240x240 +append ",
		  paste(sfiles, collapse=" ")," ",sfile, sep=""))
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width =max(350,200*lp),
					  pdf=pdf)
    frameProcesses <- list()
    cat("\ncreating frame plots...")

    threshFlag<-list()
    for(i in seq_len(lp)){   #over dupes
         threshFlag[[i]]<-rep(TRUE,ls)
         tmpVal <- unlist(tempDist[[dupes[i]]])
         tmpIndx <- which(tmpVal < mean(tmpVal[!is.na(tmpVal)]) )
		 if(is.null(absolute.value)){
				 outNames <-  names(calout.detect(tmpVal[!is.na(tmpVal)],
				 alpha=alpha,method="GESD" )$val)
				 threshFlag[[i]][which(names(tmpVal) %in% outNames)] <-FALSE       
				 threshFlag[[i]][tmpIndx] <- TRUE
		 }else{ 
				 outNames <- names(which(tmpVal >absolute.value))
		         threshFlag[[i]][which(names(tmpVal) %in% outNames)] <-FALSE     
		 }         
    }

    for(i in seq_len(ls)){ #over patient
	fnames <- NULL
        agTmp <- aggregatorList()
	for(j in seq_len(lp)){   #over dupes
	    tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
                                          gsub("\\..*$", "", j), ".pdf", sep=""))
		if(is.null(det.dimensions))
			pdf(file=tfile)
		else
     	    pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
			
	    print(tempgrph[[dupes[j]]][[patientID[i]]])
	    dev.off()
	    fnames <- c(fnames, tfile)
            val <- unlist(tempDist[[dupes[j]]])
		agTmp[[j]] <- new("numericAggregator", passed=threshFlag[[j]][i],
					 x=if(is.na(tempDist[[dupes[j]]][[patientID[i]]])) as.numeric(NA) else
						as.numeric(formatC(tempDist[[dupes[j]]][[patientID[i]]],digits=4))) 				 
        cat(".")
	}
    
	names(agTmp) <- dupes
	nfail <- !sapply(agTmp, slot, "passed")
            val <- if(sum(nfail)==1) factor(2) else factor(0)
     	    if(sum(nfail)==0)
             val <- factor(1)

	ba <- new("discreteAggregator", x=val,passed =as.logical(sum(nfail)==0))
	fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir,
				  width=200, pdf=pdf)
	fid <- patientID[i]
	frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
						    summaryAggregator=ba,
						    frameAggregators=agTmp,
						    frameGraphs=fGraphs,details=list(absolute.value=absolute.value,
															alpha=alpha))
    }
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name,
                     type="Density", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))
}


qaProcess.ECDFPlot <- function(flowList,
	dyes=NULL,
	outdir="QAReport",
	alpha = 0.05,
	absolute.value=NULL,
	sum.dimensions=NULL,
	det.dimensions=NULL,
	pdf=TRUE,
	name="ECDF",...
){
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    sfiles <- NULL
    alqLen<- length(flowList)
    Pats  <- lapply(flowList, sampleNames)
    patientID <- unique(unlist(Pats))
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    if(is.null(dyes)){
		dupes <- locateDuplicatedParameters(flowList)
		if(length(dupes)==0)
			stop("Duplicated parameters do not appear in the 
                  list of flowSets provided")
    }else{
	dupes <- as.character(dyes)
    }

    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
    valRange <-data.frame()
    for (cellType in dupes ){
		formula<-paste("~","`",cellType,"`","|","Patient",sep="")
		tempgrph[[cellType]]<-list()
		tempDist[[cellType]]<-list()
		xmax<-0
		xmin<-100000
		parLbl<-vector(mode="character",length=alqLen)
		tempInput <- list() 
		for( i in patientID){
			res <-data.frame()
			pointCount <- 512
			p<-ppoints(pointCount)
			for(j in seq_len(alqLen)){
				par <- locateParameter(flowList,cellType,j,i)
				if(length(par)!=0 && !is.na(par) && nrow(flowList[[j]]@frames[[i]])!=1){
					parLbl[j] <- paste(j," ",par)
					eps<- c(.Machine$double.eps, - .Machine$double.eps)
					valRange<-range(flowList[[j]]@frames[[i]][,par])
					ranges <- t(valRange) +eps
					ef <- char2ExpressionFilter(paste("`", par, "`>", ranges[1]," & `",par,"`<",
											ranges[2], sep="",collapse=""), filterId=cellType)
					ff <-filter(flowList[[j]]@frames[[i]][,par],ef)
					value=exprs(Subset(flowList[[j]]@frames[[i]][,par],ff))
					colnames(value) <- cellType
					quant <-  quantile(x = value,probs=p)
					newres<-data.frame(Patient=rep(i, 512), Aliquot=rep(j,512),
							  dx = p, dy = quant,check.names=FALSE)
					res<-rbind(res,newres)
					valRange<-range(flowList[[j]]@frames[[i]][,par])
					tempInput[[j]] <-value
				}else{
					tempInput[[j]] <-NA
					parLbl[j] <- paste(j," ")
				} ## end of length(par)
					
			} ##end of alqLen

			if( !all(is.na(unlist(tempInput))) &&  length(which(unlist(lapply(tempInput,length)) > 1))>1 ){
				tempgrph[[cellType]][[i]]<- xyplot(dx ~ dy | Patient,data=res,groups=Aliquot,
								col=myCol[unique(res[,"Aliquot"])],ylab="Emperical CDF",xlab=as.character(cellType),
								key=simpleKey(text=parLbl,space="right", points=F,col=myCol),
								plot.points=F,pch = ".", cex = 2, type = c("l"))
				dst <- KLdist.matrix(lapply(tempInput,"/",valRange[2,]),symmetrize=TRUE)
				tempDist[[cellType]][[i]]<-sum(dst,na.rm=T)/
									  length(which(!is.na(dst)==T))
			}else{
				tempDist[[cellType]][[i]]<- NA
				m<-data.frame(x=1:2,y=1:2,z=c(0,0))
				tempgrph[[cellType]][[i]] <- ecdfplot(~z,data=m,main= "MISSING",xlab="")
			}
			cat(".")
		}  ## end of patientID
		xrange<-c(0.9*valRange[1,],1.1*valRange[2,])
		sfile <- file.path(idir, paste("summary_", cellType, ".pdf", sep=""))
		if(is.null(sum.dimensions))
			pdf(file=sfile)
		else
			pdf(file=sfile, width=sum.dimensions[1],height=sum.dimensions[2])

	   print(densityplot(~ x | patientID, data = list(patientID = 
						 factor(names(tempgrph[[cellType]]), levels = names(tempgrph[[cellType]])),
                         x = seq_along(tempgrph[[cellType]])),
                        xlim=xrange, ylim=c(0,1.2),
                        xlab=as.character(cellType),
                        ylab="ECDF",
                        key=simpleKey(text=parLbl,space="right",
                              points=F,col=myCol),
                        lwd=2,
                        panel = function(x, y, plot.list) {
                           do.call(panel.xyplot,
                                 trellis.panelArgs(
                                       tempgrph[[cellType]][[x]], 1))
                        }
                        ))      
		dev.off()	
		sfiles <- c(sfiles, sfile)
		cat(".")
	
    } ##end of dupes
    sfile <- paste(idir, "summary.pdf", sep="/")
	system(paste("convert ", " -density 240x240 +append ",
		  paste(sfiles, collapse=" ")," ",sfile, sep=""))
	# system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 # lp, "x1 ", sfile, sep=""))
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, 
				width= max(350,200*lp), pdf=pdf)

    frameProcesses <- list()
    cat("\ncreating frame plots...")
    threshFlag<-list()
    for(i in seq_len(lp)){   #over dupes
		threshFlag[[i]]<-rep(TRUE,ls)
        tmpVal <- unlist(tempDist[[dupes[i]]])
		tmpIndx <- which(tmpVal < mean(tmpVal[!is.na(tmpVal)]))
		if(is.null(absolute.value)){
				 outNames <-  names(calout.detect(tmpVal[!is.na(tmpVal)],
				 alpha=alpha,method="GESD" )$val)
				 threshFlag[[i]][which(names(tmpVal) %in% outNames)] <-FALSE       
				 threshFlag[[i]][tmpIndx] <- TRUE
		 }else{ 
				 outNames <- names(which(tmpVal >absolute.value))
		         threshFlag[[i]][which(names(tmpVal) %in% outNames)] <-FALSE     
		 }		
    }								      

    for(i in seq_len(ls)){ #over patient
            fnames <- NULL
            agTmp <- aggregatorList()
	    	for(j in seq_len(lp)){   #over dupes
                tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
                                         gsub("\\..*$", "", j), ".pdf",sep=""))
				if(is.null(det.dimensions))
					pdf(file=tfile)
				else
					pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])

                print(tempgrph[[dupes[j]]][[patientID[i]]])
                dev.off()
                fnames <- c(fnames, tfile)
                val <- unlist(tempDist[[dupes[j]]])
				agTmp[[j]] <- new("numericAggregator", passed=threshFlag[[j]][i],
					 x=if(is.na(tempDist[[dupes[j]]][[patientID[i]]])) as.numeric(NA) else
						as.numeric(formatC(tempDist[[dupes[j]]][[patientID[i]]],digits=4))) 				 
                cat(".")
            } ##end of dupes
            names(agTmp) <- dupes
			nfail <- !sapply(agTmp, slot, "passed")
            val <- if(sum(nfail)==1) factor(2) else factor(0)
     	    if(sum(nfail)==0)
				val <- factor(1)
            ba <- new("discreteAggregator", x=val,passed =as.logical(sum(nfail)==0))
            fGraphs <- qaGraphList(imageFiles=fnames,imageDir=idir,
                                     width=200, pdf=pdf)
            fid <- patientID[i]
            frameProcesses[[fid]] <- qaProcessFrame(frameID=fid,
                                                      summaryAggregator=ba,
                                                      frameAggregators=agTmp,
                                                      frameGraphs=fGraphs,
													  details=list(absolute.value=absolute.value,
																	alpha=alpha))
		    cat(".")
    }
    cat("\n")
    return(qaProcess(id=gid, name=name,
                     type="ECDF", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))
}


qaProcess.KLDistPlot <- function(
	flowList,
	dyes=NULL,
	outdir="QAReport",
	alpha=0.05,
	absolute.value=NULL,
	sum.dimensions=NULL,
	det.dimensions=NULL,
 	pdf=TRUE,
	name="KLDist", ...
){
    cat("creating summary plots...")
    gid <- guid()
    idir <- createImageDir(outdir, gid)
    sfiles <- NULL
    alqLen<- length(flowList)
    Pats  <- lapply(flowList, sampleNames)
    patientID <- unique(unlist(Pats))
    myCol<- colorRampPalette(brewer.pal(9, "Set1"))(alqLen)
    if(is.null(dyes)){
	dupes <- locateDuplicatedParameters(flowList)
    if(length(dupes)==0)
			stop("Duplicated parameters do not appear in the 
                  list of flowSets provided")
    }else{
	dupes <- as.character(dyes)
    }

    lp<-length(dupes)
    ls <- length(patientID)
    tempgrph<-list()
    tempDist<-list()
    sfiles <- NULL
    colorFun<-colorRampPalette(c("yellow","red"))
	for (cellType in dupes ){
		formula<-paste("~","`",cellType,"`","|","Patient",sep="")
		tempgrph[[cellType]]<-list()
		tempDist[[cellType]]<-list()
        parLbl<-vector(mode="character",length=alqLen)
        outRes<-data.frame()
		for( i in patientID){
            res<-data.frame()
            tempList<- list()
            for(j in seq_len(alqLen)){
                par <- locateParameter(flowList,cellType,j,i)
                if(length(par)!=0 && !is.na(par) && nrow(flowList[[j]]@frames[[i]])!=1 ){
					parLbl[j] <- paste(j," ",par)
					eps<- c(.Machine$double.eps, - .Machine$double.eps)
					valRange<-range(flowList[[j]]@frames[[i]][,par])
					ranges <- t(valRange) +eps
					ef <- char2ExpressionFilter(
					paste("`", par, "`>", ranges[1]," & `",par,"`<",
					ranges[2], sep="",collapse=""), filterId=cellType)
					ff <-filter(flowList[[j]]@frames[[i]][,par],ef)
					summary(ff)    	                
					value=exprs(Subset(flowList[[j]]@frames[[i]][,par],ff))
                           	colnames(value) <- cellType
                           	tempList[[j]]<-value
                }else{
					tempList[[j]]<-NA
                    parLbl[j] <- paste(j," ")
				}
            }## end of alqLen
			if( !all(is.na(unlist(tempList))) &&  length(which(unlist(lapply(tempList,length)) > 1))>1 ){
				dst <- KLdist.matrix(lapply(tempList,"/",valRange[2,]),symmetrize=TRUE) 
				tempDist[[cellType]][[i]]<-sum(dst,na.rm=T)/length(which(!is.na(dst)==T))
   	        	pm<-as.matrix(dst)
				diag(pm)<-NA
        		z<-as.matrix(as.vector(pm),ncol=1)
	        	x <- as.character(sapply(parLbl,function(x){
        		   rep(x,length(parLbl))
        			}))
        		y <- rep(parLbl,length(parLbl))
        		res <- data.frame(x=factor(x,levels=parLbl),y=factor(y,levels=parLbl),z=z)
				tempgrph[[cellType]][[i]]<-levelplot(z~x*y,data=res,xlab="Aliquot",ylab="Aliquot",
                                           main=cellType,scales = list(x = list(rot = 90)),
                                           col.regions=colorFun,colorkey=list(col=colorFun))    
    						Patient=rep(i,nrow(res))
					        tm<-cbind(res,Patient)
					        outRes<-rbind(tm,outRes)
			}else{
				tempDist[[cellType]][[i]] <- NA
				m<-data.frame(x=1,y=1,z=1)
				tempgrph[[cellType]][[i]] <- levelplot(z~x*y,data=m,main= "MISSING",xlab="",colorkey=list(col=colorFun))
			}
			cat(".")
		} ## end of patientID
        sfile <- file.path(idir, paste("summary_", cellType, ".pdf", sep=""))
		if(is.null(sum.dimensions))
			pdf(file=sfile)
		else
			pdf(file=sfile, width=det.dimensions[1],height=det.dimensions[2])
        print(grph<-levelplot(z~x*y|Patient,data=outRes,xlab="Aliquot",ylab="Aliquot",
                               main=cellType,scales = list(x = list(rot = 90)),
                               col.regions=colorFun,colorkey=list(col=colorFun)))
        dev.off()	
        sfiles <- c(sfiles, sfile)
        cat(".")
    } ## end of cellType

    sfile <- paste(idir, "summary.pdf", sep="/")
	system(paste("convert ", " -density 240x240 +append ",
		  paste(sfiles, collapse=" ")," ",sfile, sep=""))
    # system(paste("montage ", paste(sfiles, collapse=" "), " -geometry +0+0 -tile ",
                 # lp, "x1 ", sfile, sep=""))
    sgraph <- qaGraph(fileName=sfile, imageDir=idir, width=max(350,200*lp),pdf=pdf)
    frameProcesses <- list()
    cat("\ncreating frame plots...")
    threshFlag<-list()
    for(i in seq_len(lp)){   #over dupes
		threshFlag[[i]]<-rep(TRUE,ls)
		tmpVal <- unlist(tempDist[[dupes[i]]])       
		tmpIndx <- which(tmpVal < mean(tmpVal[!is.na(tmpVal)]))
		if(is.null(absolute.value)){
			outNames <-  names(calout.detect(tmpVal[!is.na(tmpVal)],
			alpha=alpha,method="GESD" )$val)
			threshFlag[[i]][which(names(tmpVal) %in% outNames)] <-FALSE       
			threshFlag[[i]][tmpIndx] <- TRUE
		}else{ 
			outNames <- names(which(tmpVal > absolute.value))
			threshFlag[[i]][which(names(tmpVal) %in% outNames)] <-FALSE     
		} 
    }

    for(i in seq_len(ls)){ #over patient
		fnames <- NULL
        agTmp <- aggregatorList()
		for(j in seq_len(lp)){   #over dupes	
			tfile <- file.path(idir, paste("frame_", sprintf("%0.4d", i), "_",
									 gsub("\\..*$", "", j), ".pdf",sep=""))
			if(is.null(det.dimensions))
				pdf(file=tfile)
			else
				pdf(file=tfile, width=det.dimensions[1], height=det.dimensions[2])
			print(tempgrph[[dupes[j]]][[patientID[i]]])
			dev.off()
			fnames <- c(fnames, tfile)
			val <- unlist(tempDist[[dupes[j]]])
			agTmp[[j]] <- new("numericAggregator", passed=threshFlag[[j]][i],
					 x=if(is.na(tempDist[[dupes[j]]][[patientID[i]]])) as.numeric(NA) else
						as.numeric(formatC(tempDist[[dupes[j]]][[patientID[i]]],digits=4))) 				 
			cat(".")
		}
    
		names(agTmp) <- dupes
		nfail <- !sapply(agTmp, slot, "passed")
        val <- if(sum(nfail)==1) factor(2) else factor(0)
     	if(sum(nfail)==0)
            val <- factor(1)
		ba <- new("discreteAggregator", x=val,passed =as.logical(sum(nfail)==0))
		fGraphs <- qaGraphList(imageFiles=fnames, imageDir=idir, width=200, pdf=pdf)
		fid <- patientID[i]
		frameProcesses[[fid]] <- qaProcessFrame(frameID=fid, summaryAggregator=ba,
						    frameAggregators=agTmp, frameGraphs=fGraphs,
							details =list(absolute.value=absolute.value, alpha=alpha)
							)
		}
    ## create qaProcess object
    cat("\n")
    return(qaProcess(id=gid, name=name,
                     type="KLD", summaryGraph=sgraph,
                     frameProcesses=frameProcesses))

}

