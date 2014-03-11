## Set up a HTML header
myOpenHtmlPage <- function(name, title = "", path="../") 
{
    name = gsub("\\.html$", "", name)
    con = file(paste(name, ".html", sep = ""), open = "wt")
    writeLines(paste("<html>\n<head>\n<title>", title, "</title>", sep=""),
               con)
    writeLines(paste("<link rel=\"stylesheet\" type=\"text/css\" ",
                     "href=\"", path, "qaReport.css\">", sep=""), con)
    writeLines(paste("<script src=\"", path, "qaReport.js\" ",
                     "type=\"text/javascript\"></script>", sep=""), con)
    writeLines(paste("</head>\n<body style=\"font-family: helvetica,arial",
                     ",sans-serif;\">", sep=""), con)
    return(con)
}

## include links to pdf versions of images if pdf=TRUE and if the vectorized
## version of the image exists.
pdfLink <- function(vimg, bimg, class, id, pdf=TRUE, qaGraph=NULL)
{
    res <- character(length(vimg))
    for(i in seq_along(vimg))
    {
        cl <- ifelse(missing(class), "'", sprintf(" %s'", class))
        res[[i]] <- if(pdf && !is.null(qaGraph) &&
                       file.exists(qaGraph[[i]]@fileNames["vectorized"]))
            sprintf(paste("<a href='%s' target='QAdetails'>\n<img class='%s",
                          "src='%s'id='img_%s'>\n</a>\n"),
                    vimg[i], cl, bimg[i], id[i]) else
        sprintf("<img class='nolink%s src='%s' id='img_%s'>\n",
                cl, bimg[i], id[i])
    }
    return(res)
}



## make sure the QAprocesses and the flowSets match
checkInputs <- function(sets, procs, gproc, grouping)
{    
    

   ## Force inputs to be lists
    single <- FALSE
    if(!is.list(sets))
    {
        single <- TRUE              ## single flowSet
        sets <- list(sets)
        procs <- list(procs)
        if(!is.null(gproc) && !is.list(gproc))
            gproc <- list(gproc)
    }

    
 

    if(!is.null(procs))
    {
       
        all.list <- is.list(procs) & length(procs) == length(sets)
        procs <- lapply(procs, function(x) if(is(x, "qaProcess")) list(x) else x)
        #correct.class <- all(sapply(sets, is, "flowSet")) &&
        #all(sapply(procs, function(x)  sapply(x,is, "qaProcess")))
        if(!all.list)
            stop("When 'set' is a list of flowSets, 'processes' needs to be ",
                 "a list of lists of qaProcess objects, where each list corresponds ",
                 "to the results for a particular flowSet.")
        for(i in seq_along(procs))
        {  
            sn <- sampleNames(sets[[i]])
            qn <- matrix(sapply(procs[[i]], function(x) names(x@frameProcesses)),
                         nrow=length(procs[[i]]), byrow=TRUE)
            if(!all(apply(qn, 2, function(x) length(unique(x))==1)))
                stop("Some of the qaProcess objects in processes[[", i,
                     "]] are not compatible.")
            if(length(intersect(sn, qn[1,])) != length(sn))
                stop("The sampleNames of set[[", i, "]] do not match the names ",
                     "of the qaProcess frames in processes[[", i, "]].")
        }
    }
   # else
   # {
   #     single <- TRUE
   #     sets <- if(!is.list(sets)) list(sets) else sets[1]
   # }

     ## Add the globalProcess to the list of processes and also add a dummy 
    ## flowSet to the list of flowSets if needed
 
	if(!is.null(gproc))
    {  
        dfs <- flowCore:::copyFlowSet(sets[[1]])
        opd <- pData(dfs)
        gpat <- if(is.null(grouping)) "sampleid" else paste("sampleid", grouping, sep="|")
        csel <- grep(gpat, tolower(colnames(opd)))
        pData(dfs) <- opd[,csel, drop=FALSE]
        opa <- parameters(dfs[[1]])
        popa <- pData(opa)
        popa[,-1] <- NA
        pData(opa) <- popa
        dfs[[1]]@parameters <- opa
        if(!single)
            gproc <- list(gproc) 
        sets <- if(is.null(procs[[1]])) list("comparison across panels"=dfs) else
        c("comparison across panels"=dfs, sets)
        procs <- if(is.null(procs[[1]])) (gproc) else c(gproc, procs)
    }
 
   return(list(sets=sets, procs=procs, single=single, gproc=!is.null(gproc[[1]])))
}



## Copy the images contained within the qaGraph objects of the process list to
## the output directory (which is only needed when they have been created some
## place else...)
copyGraphs <- function(procs, outdir)
{
    procs <- unlist(procs)
    for(p in procs)
    {   
        relBase <- file.path(outdir, "images",  p@id)
        if(!file.exists(relBase))
            dir.create(relBase, recursive=TRUE)
        sg <- p@summaryGraph@fileNames
        file.copy(sg, relBase)
        for(f in p@frameProcesses)
        {
            sfg <- f@summaryGraph@fileNames
			
			if(length(sfg) >0){
				file.copy(sfg, relBase)
			}
            fg <- sapply(f@frameGraphs, slot, "fileNames")
            if(length(fg) >0)
                file.copy(fg, relBase)
        }
    }
    invisible()
}




## create HTML report for (lists of) QA processes 
writeQAReport  <- function(set, processes=NULL, globalProcess=NULL, outdir="./qaReport",
                           grouping=NULL, pagebreaks=TRUE,
                           pdf=TRUE)
{    
    ## making sure the inputs are correct
    inputs <- checkInputs(set, processes, globalProcess, grouping)
    checkClass(outdir, "character", 1)
    if(file.exists(file.path(outdir, "index.html")))
        warning("Target directory already exists. Content may be ",
                "overwritten")
    if(!file.exists(file.path(outdir, "images")))
        dir.create(file.path(outdir, "images"), recursive=TRUE)
    checkClass(pagebreaks, "logical", 1)
    checkClass(pdf, "logical", 1)
   
    ## We only need panel tabs if 'set' is a list of 'flowSets'
    ## FIXME: Need to somehow factor in the globalProcess argument
    single <- inputs$single
    set <- inputs$sets
    processes <- inputs$procs


    ## For the overview page, we need to match samples across panels
    sID <- all(sapply(set, function(x) "sampleid" %in% 
                   tolower(colnames(pData(x)))))
    if(!single && !sID)
        warning("Some of the panels in 'set' don't have a global ",
                "sample identifier.\nUnable to create overview.")
    
    ## copy infrastructure
    sdir <- system.file("htmlTemplates", package = "flowQ")
    idir <- file.path(outdir, "images")
    file.copy(dir(sdir, full.names=TRUE), idir, overwrite=TRUE)

    ## We have to make sure that all images are copied to the output directory
	
    copyGraphs(processes, outdir)
    
    ## iterate over panels
   
    for(s in seq_along(set)){
        
        ## rearange set according to grouping if necessary
        grps <- NULL
        if(!is.null(grouping)){
            if(!is.character(grouping) ||
               !grouping %in% varLabels(phenoData(set[[s]])))
                stop("'grouping' must be a factor variable in the flowSet ",
                     "phenoData")
            sord <- order(pData(set[[s]])[, grouping])
            if(!all(sord==seq_len(length(set[[s]]))))
                set[[s]] <- set[[s]][sord]
            grps <- pData(set[[s]])[, grouping]
        }

        ## open a file connection  
        ifile <- if(!sID && s==1) "index" else if(!sID)
            paste("index", s-1, sep="") else  paste("index", s, sep="")
        con <- myOpenHtmlPage(file.path(outdir, ifile), "qatest", "images/")
        
        ## setup of table and table header row
                                        # process <- if(length(set)>1) processes[[s]] else processes
        
        process <- processes[[s]]
        writeLines("\n<table class=\"QA\">", con)
        pIDs <- sapply(process, slot, "id")
        pNames <- sapply(process, slot, "name")
        pTypes <- sapply(process, slot, "type")
        relBase <- file.path("images", pIDs)
        sumLinks <- file.path(relBase, sapply(process, function(x)
                                              names(x@summaryGraph)))
        sumVecLinks <- gsub("\\..*$", ".pdf", sumLinks)
        nrAggr <- sapply(process, function(x)
                         length(x@frameProcesses[[1]]@frameAggregators))+1
        th <- paste("\n<th colspan=\"", nrAggr, "\" ",
                    "id=\"", pIDs, "_sumHeader\">\n",
                    "<div id=\"", pIDs, "_button",
                    "\" onClick=\"toggleImage('", pIDs, "')\">\n",
                    pNames, "\n</div>\n</th>", sep="")
        esel <- sapply(process, function(x) length(x@summaryGraph@fileNames)>0)
        th[!esel] <- paste("\n<th colspan=\"",
                           nrAggr[!esel], "\" ",
                           "id=\"", pIDs[!esel], "_sumHeader\">\n",
                           "<div id=\"", pIDs[!esel],
                           "_button",
                           "\">\n", pNames[!esel], "\n</div>\n</th>", sep="",
                           collapse="\n")
        th <- paste(th, collapse="\n")
        pd <- pData(set[[s]][[1]]@parameters)[,c("name", "desc", "minRange",
                                                 "maxRange")]
        writeLines(paste("\n<tr class=\"QAHeader\">\n\n<th>\n",
                         "<div id=\"parameters_button",
                         "\" onClick=\"toggleImage('parms')\">\n",
                         "flow set details\n</div>\n",  
                         "</th>\n", th, "\n\n</tr>", sep=""), con)
        graphs <- lapply(process, function(x) x@summaryGraph)
        td <- paste("\n<td colspan=\"", nrAggr, "\"",
                    " id=\"", pIDs, "_sumBack\">\n",
                    pdfLink(vimg=sumVecLinks, bimg=sumLinks, id=pIDs, pdf=pdf,
                            qaGraph=graphs),
                    "</td>", sep="", collapse="\n")
        writeLines(paste("\n<tr class=\"QASummary\">\n\n<th>",
                         "\n<span id=\"img_parms\" style=\"display:none;\">",
                         sep=""), con)
        writeLines(pd, con)
        writeLines(paste("</span>\n</th>\n", td, "\n\n</tr>\n", sep=""), con)
        
        frameIDs <- sampleNames(set[[s]])
        classes <- paste("QAFrameHeader",
                         c(" even", " odd")[(seq_along(frameIDs)%%2)+1], sep="")
        names(classes) <- frameIDs
        fpp <- if(pagebreaks) 14 else 10e20
        lf <- length(frameIDs)
        nrPages <- (lf %/% fpp) +1
        counter <- 1
        page <- 1
        lastGrp <- grps[1]
        for(f in frameIDs){
            showRow <- ifelse(counter>fpp, "none", "table-row")
            pd <- pData(set[[s]])[f,,drop=FALSE]
            phi <- paste("<span id=\"img_pd_", counter,
                         "\" style=\"display:none;\"",
                         ">", sep="", collapse="\n")
            ## new table row and column header for aggregators
            if(is.null(grouping)){## no grouping is specified
                writeLines(paste("\n<tr class=\"", classes[f], "\" id=\"frow1_",
                                 counter, "\" style=\"display:",
                                 showRow, ";\">\n\n<th id=\"",
                                 f, "_sumHeader\">\n<div class=\"QARowButton\"",
                                 " id=\"", f, "_button\" onClick=\"toggle",
                                 "Image('pd_", counter, "')\">\n<div class=\"",
                                 "QAFrameHeaderNr\">", counter, "</div>\n<span ",
                                 "class=\"QAFrameHeaderID\">", f, "</span>",
                                 "\n", phi, sep=""), con)
                writeLines(pd, con)
            }else{##grouping
                caption <-  paste("<div class=\"QAFrameHeaderNr\">",
                                  counter, "</div><span class=\"QAFrameHeader",
                                  "ID\">", f, "</span>",
                                  "<span class=\"QAFrameHeaderGrp\">",
                                  grps[counter], "</span>", sep="")
                if(counter==1 || grps[counter]!=lastGrp){
                    writeLines(paste("\n<tr class=\"", classes[f],
                                     "\" id=\"frow1_",
                                     counter, "\" style=\"display:", showRow,
                                     ";\">\n\n<th id=\"",
                                     f, "_sumHeader\">\n<div class=\"QAGrp",
                                     "Header\">", grps[counter],
                                     "</div>\n<div ",
                                     "class=\"QARowButton\" id=\"", f,
                                     "_button\" onClick=\"toggleImage(",
                                     "'pd_", counter, "')\">", caption,
                                     "\n</div>\n", phi, sep=""), con)
                    writeLines(pd, con)                      
                }else{  
                    writeLines(paste("\n<tr class=\"", classes[f],
                                     "\" id=\"frow1_",
                                     counter, "\" style=\"display:",
                                     showRow, ";\">\n\n<th id=\"",
                                     f, "_sumHeader\">\n<div class=\"QARow",
                                     "Button\" id=\"", f, "_button\" ",
                                     "onClick=\"toggle",
                                     "Image('pd_", counter, "')\">", caption,
                                     "\n</div>\n", phi, sep=""), con)
                    writeLines(pd, con)
                }
            }
            writeLines(paste("</span>\n</div>\n</th>\n"), con)
            
            ## aggregators first, each process is one column
            for(p in seq_along(process)){
                
                thisProcess <- process[[p]]@frameProcesses[[f]]
                pid <- thisProcess@id
                sid <- thisProcess@summaryGraph@id
                writeLines(paste("\n<td class=\"QASumAggr\" id=\"", pid,
                                 "_sumHeader\" align=\"center\">", sep=""), con)
                nrDetAgr <- length(thisProcess@frameAggregators)
                offset <- ""
                ## add trigger to access details
                if((counter==1 || (counter-1) %% fpp == 0) && nrDetAgr>0){
                    writeLines(paste("<div class=\"QADetTrigger\" id=\"",
                                     pIDs[p], "_detTriggerIn_", page, "\" ",
                                     "onClick=\"",
                                     "toggleDetails(", nrDetAgr, ", ",
                                     length(frameIDs), ", ", p, 
                                     ", '", pIDs[p], "', ", nrPages, ")\">",
                                     "\n+\n</div>", sep=""), con)
                    writeLines(paste("<div class=\"QADetTrigger\" id=\"",
                                     pIDs[p], "_detTriggerOut_", page, "\" ",
                                     "onClick=\"",
                                     "toggleDetails(", nrDetAgr, ", ",
                                     length(frameIDs), ", ", p, 
                                     ", '",  pIDs[p], "', ", nrPages, ")\"",
                                     " style=\"display:none;\">\n&#150;\n</div>",
                                     sep=""), con)
                    offset <- paste(" style=\"position:relative; left:-7px;",
                                    "margin-left:15px; margin-right:0px;",
                                    "z-index:0;\"")
                    page <- page+1
                }
                ## add summary aggregator and link to image if necessary
                if(length(sid)>0)
                    writeLines(paste("<div class=\"QAFrameButton\" ",
                                     "id=\"", pid, "_button\" onClick=\"",
                                     "toggleImage('", sid, "')\"", offset,
                                     ">", sep=""), con)
                else
                    writeLines(paste("<div class=\"QAFrameButtonNoSel\"",
                                     offset,
                                     ">", sep=""), con)
                writeLines(thisProcess@summaryAggregator, con)
                writeLines("</div>\n</td>", con)
                ## add detailed aggregators and links to images if necessary
                for(d in seq_along(thisProcess@frameAggregators)){
                    fname <- names(thisProcess@frameAggregators)[d]
                    fname  <- ifelse(is.null(fname), "", fname)
                    did <- thisProcess@frameGraphs[[d]]@id
                    id <- paste(p, d, counter, sep="_")
                    writeLines(paste("\n<td class=\"QADetAggr\" id=\"row_", id,
                                     "_1\" align=\"center\">", sep=""), con)
                    if(length(did)>0)
                        writeLines(paste("<div class=\"",
                                         "QADetButton\" ", "id=\"button_",id,
                                         "\" onClick=\"toggleImage('", id,
                                         "')\">", fname, sep=""), con)
                    else
                        writeLines("<div>", con)
                    writeLines(thisProcess@frameAggregators[[d]], con)
                    writeLines("</div>\n</td>", con)
                }## end for d
            }## end for p
            writeLines("</tr>", con)
            
            ## new table row and column header for images
            writeLines(paste("\n<tr class=\"", classes[f], "\" id=\"frow2_",
                             counter, "\" style=\"display:", showRow,
                             ";\">\n\n<td ",
                             "class=\"QASumGraph\">\n</td>", sep=""), con)
            ## now the images, each process is one column
            for(p in seq_along(process)){
                thisProcess <- process[[p]]@frameProcesses[[f]]
                pid <- thisProcess@id
                sid <- thisProcess@summaryGraph@id
                writeLines(paste("\n<td class=\"QASumGraph\" id=\"",
                                 pid, "_sumBack\" align=\"center\">",
                                 sep=""), con)
                if(length(sid)>0){
                    sGraph <- file.path(relBase[p],
                                        names(thisProcess@summaryGraph))
                    sVecGraph <- gsub("\\..*$", ".pdf", sGraph)
                    writeLines(pdfLink(vimg=sVecGraph, bimg=sGraph, class="QASumGraph",
                                       id=sid, pdf=pdf,
                                       qaGraph=list(thisProcess@summaryGraph)),
                               con)
                }
                writeLines("</td>", con)
                for(d in seq_along(thisProcess@frameAggregators)){
                    did <- thisProcess@frameGraphs[[d]]@id
                    id <- paste(p, d, counter, sep="_")
                    writeLines(paste("\n<td class=\"QADetGraph\" id=\"row_", id,
                                     "_2\" align=\"center\">", sep=""), con)
                    if(length(did)>0){
                        fGraph <- file.path(relBase[p],
                                            names(thisProcess@frameGraphs[[d]]))
                        fVecGraph <- gsub("\\..*$", ".pdf", fGraph)
                        writeLines(pdfLink(vimg=fVecGraph, bimg=fGraph, class="QADetGraph",
                                           id=id, pdf=pdf,
                                           qaGraph=list(thisProcess@frameGraphs[[d]])),
                                   con)
                    }
                    writeLines("</td>", con)
                }## end for d
            }## end for p
            lastGrp <- grps[counter]
            counter <- counter+1
            writeLines("\n</tr>\n", con)
        }## end for f
        writeLines("\n</table>", con)
        
        ## the page navigation
        writeLines(paste("<div class=\"QAPagesTile\"><table width=\"100%\"",
                         " style=\"padding-right:35px;\"><tr><td align=",
                         "\"left\">", sep=""), con)
        if(lf>fpp){
            from <- c(seq(1, lf, fpp))
            nt <- length(from)
            to <- c(seq(fpp, lf, fpp), lf)[1:nt]
            writeLines(paste("<span class=\"QAPages\" id=\"pages_", 1:nt,
                             "\" onClick=\"togglePages(", from, ", ", to, ", ",
                             nt, ", ", lf, ")\">Frames ", from, "-", to,
                             "</span>", sep=""), con)
        }
        
        np <- length(set)
        iFiles <- if(!sID && np>1) c("", 1:(np-1)) else as.character(1:np)
   
        if(np>1){

            writeLines("\n</td><td align=\"right\">", con)
            if(!inputs$gproc) 
               pnams <- paste("Panel ", 1:np) 
                     else
               pnams <-   c("Multipanel", paste("Panel ", 1:(np-1)))
            if(!is.null(names(set)[s]) && nchar(names(set)[s])>0)
                pnams[s] <- paste(pnams[s], " <i><small>(", names(set)[s],
                               ")</i></small>", sep="")
            panels <- paste("<span class=\"QAPanels\" id=\"panels_", 1:np,
                            "\"n><a class=\"QAPanels\" href=\"index",
                            iFiles, ".html\">",
                            pnams, "</a></span>", sep="")
            if(!is.null(names(set)))
                panels[s] <- gsub("QAPanels","QAPanelsAct", panels[s])
            if(sID)
                panels <- c(paste("<span class=\"QAPanels\" id=\"panels_0",
                                  "\"n><a class=\"QAPanels\" href=\"index",
                                  ".html\">Summary</a></span>", sep=""),
                            panels)
            writeLines(panels, con) 
        }else{
        	pnams <- paste(" ") 

        }
        writeLines("</td></tr></table></div>", con)  
        closeHtmlPage(con)
    }##end s

    ## We create an overview page if we have multiple panels
    nps <- sapply(processes, length)
    oneOnly <- all(nps==1)
    if(sID && !oneOnly)
    {
        #if(single)
        #    processes <- list(processes)
        ifile <- "index"
        con <- myOpenHtmlPage(file.path(outdir, ifile), "qatest", "images/")
        on.exit(closeHtmlPage(con))
        summary <- failedProcesses(processes, set, pnams)
        writeLines(summary, con)
    }
    else
    {
        file.copy(file.path(outdir, "index1.html"), file.path(outdir, "index.html"))
    }
    return(file.path(outdir, "index.html"))
}



## Create QA output based on a flowSet and a list of QA functions.
## This is a very basic convenience function, for more complex experiments
## including panels use writeQAReport directly
qaReport <- function(set, qaFunctions, outdir="./qaReport", argLists,
                     grouping=NULL, ...)
{
    processes <- list()
    for(i in seq_along(qaFunctions)){
        cat(paste("quality process ", i, ":\n", sep=""))
        if(missing(argLists))
            processes[[i]] <- do.call(qaFunctions[i],
                                      list(set=set, outdir=outdir,
                                           grouping=grouping))
        else{
            argLists[[i]]$set <- set
            argLists[[i]]$outdir <- outdir
            argLists[[i]]$grouping <- grouping
            processes[[i]] <- do.call(qaFunctions[i], argLists[[i]])
        }
        save(processes, file=file.path(outdir, "processes.rda"))
    }
    writeQAReport(set, processes, outdir=outdir, grouping=grouping, ...)
}





## count numbers of failed qaProcesses in a list of lists of such objects.
## Each item in the outer list is a list of qaProcess objects for a single
## panel. The output is an object for which HTML output can be generated
## via a writeLines method.
## The phenoData of the list of flowSets that gets passed as the second
## argument has to contain a column sampleIDs which provides the mapping
## of samples over panels.
failedProcesses <- function(processes, set, pnams)
{
    ## some sanity checking first
    sampleIDs <- lapply(set, function(x) {
        pd <- pData(x)
        pdid <- match("sampleid", tolower(colnames(pd)))
        pd[,pdid]})
    if(!all(listLen(lapply(sampleIDs, unique))==listLen(sampleIDs)))
        stop("'SampleIDs' must be unique in each panel")
    sids <- unlist(sampleIDs)
    comSids <- unique(sids)
    fids <- unlist(lapply(set, function(x) rownames(pData(x))))
    ##if(any(duplicated(fids)))
    ##    stop("The 'FrameIDs' in the whole experiment are not unique")
    #allChannels <- c(unique(unlist(lapply(set, colnames))), "global")
    allChannels <- c(unique(unlist(lapply(unlist(processes), cnams))), "global")
    res <- ranges <- mapping <- vector(length(set), mode="list")
    ## iterate over panels
    for(i in seq_along(processes)){
        fmat <- matrix(0, ncol=length(allChannels), nrow=length(comSids),
                       dimnames=list(comSids, allChannels))
        clist <- slist <- NULL
        nrSum <- 0
        ## iterate over qaProcess objects for one panel
        for(j in seq_along(processes[[i]])){
            nrSamp <- 0
            mlist <- NULL
            ## iterate over samples in the qaProcess object
            for(pro in processes[[i]][[j]]@frameProcesses){
                ## match frameID to sampleID
                samp <- sids[match(pro@frameID, fids)]
                mlist <- rbind(mlist, c(samp, pro@frameID, nrSamp+1))
                ## check for the multiple channels and iterate over those 
                    channels <- names(pro@frameAggregators)
                    clist <- c(clist, channels)
                    
                    for(chan in seq_along(pro@frameAggregators)){
                        fmat[samp, channels[chan]] <-
                            fmat[samp, channels[chan]] + 
                                as.numeric(!pro@frameAggregators[[chan]]@passed)
                    }
                slist <- c(slist, samp)
                nrSamp <- nrSamp+1
            }
        }
        colnames(mlist) <- c("sample", "frame", "number")
        res[[i]] <- fmat
        mapping[[i]] <- mlist
        ## We want to sum up over qaProcesses for one panel
        rtemp <- numeric(length(allChannels))
        names(rtemp) <- allChannels
        if(!is.null(clist)){
            tc <- table(clist)/nrSamp
            rtemp[names(tc)] <- tc
        }
        rtemp["global"] <- nrSum/nrSamp
        ranges[[i]] <- rtemp
    }
    ## we also want an overall summary
    sum <- res[[1]]
    if(length(res)>1)
        for(i in 2:length(res))
            sum <- sum+res[[i]]
    names(res) <- names(ranges) <- names(mapping) <- names(set)
    os <- cbind(rowSums(sum), sapply(res, rowSums))
    colnames(os) <- c("global", 1:(ncol(os)-1))
    return(new("qaProcessSummary", panels=res, summary=sum, ranges=ranges,
               mapping=mapping, pnams=pnams, overallSum=os))
}

### writes the QA report into a tab delimited file
writeQATextReport  <- function(set, processes=NULL, globalProcess=NULL, fileName="textReport.txt")
{    
    ## making sure the inputs are correct
    inputs <- checkInputs(set, processes,globalProcess, grouping=NULL)
    if(file.exists(fileName))
        warning("Target File already exists. Content may be ",
                "overwritten")
    single <- inputs$single
    set <- inputs$sets
    processes <- inputs$procs
    con = file(fileName, open = "wt")
	for(s in seq_along(set)){           
			## rearange set according to grouping if necessary
			grps <- NULL
			process <- processes[[s]]			
			tblList <- lapply(process,txtFormatQAObject)
			tbl <- tblList[[1]]
			
			for( x in tblList[-1]){
				
				tbl <- cbind(tbl,x)
			}		
			cat('\t',file=con)
			suppressWarnings(write.table(tbl,file=con,col.names=T,append=T,sep="\t",quote=F))	
			cat("\n",file=con)		
	}
	close(con)
}


cnams <- function(proc)
  unique(unlist(lapply(proc@frameProcesses, function(x) names(x@frameAggregators))))
