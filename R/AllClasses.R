## ===========================================================================
## virtual outlier test
## ---------------------------------------------------------------------------
setClass("outlier", 
         representation("VIRTUAL",test="character",
                        parameters="ANY"))


## ===========================================================================
## outlierResult
## ---------------------------------------------------------------------------
## A container for the results of outliers test 
## ---------------------------------------------------------------------------
setClass("outlierResult",
         representation(frameId="character", filterDetails="list"),
         contains="outlier",
         prototype=list(frameId=character(0), filterDetails=list()))




## ===========================================================================
## Virtual qaAggregator
## ---------------------------------------------------------------------------
## A class describing an aggregated QA value for a single flowFrame. Derived
## subclasses describe the various subtypes of aggregators. Slot 'frameID'
## stores the reference to the flowFrame (i.e., the sample name in the
## respective flowSet) and slot 'passed' contains a logical
## indicating whether the QA requirements have been met. Dedicated write
## methods of these subclasses produce the appropriate HTML output.
## ---------------------------------------------------------------------------
setClass("qaAggregator",
         representation("VIRTUAL", passed="logical"),
         prototype=list(passed=TRUE))



## ===========================================================================
## binaryAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the most simple binary type,
## i.e. indicating whether a QA requirement has been passed or not 
## ---------------------------------------------------------------------------
setClass("binaryAggregator",
         contains="qaAggregator")

binaryAggregator <- function(passed=TRUE)
    new("binaryAggregator", passed=passed)



## ===========================================================================
## discreteAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the most discrete type,
## i.e. it can have three states: passed(1), warning(2) and failed(0) 
## ---------------------------------------------------------------------------
setClass("discreteAggregator",
         representation(x="factor"),
         contains="qaAggregator",
         prototype=list(x=factor(1)))

discreteAggregator <- function(x)
{
    if(!is.factor(x)){
        if(! x %in% 0:2)
            stop("'x' must be in 0,1 or 2", call.=FALSE)
        x <- factor(x, levels=0:2, ordered=TRUE)
    }else if(!all(levels(x) %in% as.character(0:2)))
        stop("x need to be factor with levels 0,1 and 2", call.=FALSE)
    if(x == 1)
        new("discreteAggregator", x=x, passed=TRUE)
    else
        new("discreteAggregator", x=x, passed=FALSE)
}



## ===========================================================================
## factorAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the factor type, i.e.
## indicating the states of the QA results from a selection of
## different outcomes
## ---------------------------------------------------------------------------
setClass("factorAggregator",
         representation(x="factor"),
         contains="qaAggregator")

factorAggregator <- function(x, passed=TRUE)
{
    if(!missing(x) && !is.factor(x))
        x <- factor(x)
    new("factorAggregator", x=x, passed=passed)
}




## ===========================================================================
## stringAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the string type, i.e. a
## character vector that was created by the QA process with a textual
## description of the result
## ---------------------------------------------------------------------------
setClass("stringAggregator",
         representation(x="character"),
         contains="qaAggregator")

stringAggregator <- function(x, passed=TRUE) new("stringAggregator", x=x,
                                passed=passed)



## ===========================================================================
## numericAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the numeric type, i.e. a
## numeric skalar that was created by the QA process
## ---------------------------------------------------------------------------
setClass("numericAggregator",
         representation(x="numeric"),
         contains="qaAggregator")

numericAggregator <- function(x, passed=TRUE) new("numericAggregator", x=x,
                                 passed=passed)



## ===========================================================================
## rangeAggregator
## ---------------------------------------------------------------------------
## A class describing aggregated QA value of the range type, i.e. a
## numeric value within a defined range of values (e.g. a percentage)
## ---------------------------------------------------------------------------
setClass("rangeAggregator",
         representation(min="numeric", max="numeric"),
         contains="numericAggregator")
         
rangeAggregator <- function(x, min, max, passed=TRUE)
{
    if(x<min || x>max)
        stop("'x' must be in the range of 'min' and 'max'")
    new("rangeAggregator", min=min, max=max, x=x, passed=passed)
}


## ===========================================================================
## aggregatorList
## ---------------------------------------------------------------------------
## A list of aggregators for a whole flow set. All elements of the list must
## be aggregators. The class should be initiated via it's
## constructor
## ---------------------------------------------------------------------------
setClass("aggregatorList",
         contains="list")

setMethod("initialize", "aggregatorList",
          function(.Object, ...) {
              if(length(list(...))>0){
                  if(is.list(..1))
                      input <- ..1
                  else
                      input <- list(...)
                  if(!all(sapply(input, is, "qaAggregator")))
                      stop("All items of an aggregator list must ",
                           "inherit from class 'qaAggregator'")
                  .Object@.Data=input
              }
              return(.Object)
          })

aggregatorList <- function(...)
    new("aggregatorList", ...)



## ===========================================================================
## qaGraph
## ---------------------------------------------------------------------------
## A class that describes graphical output of a QA operation for a single
## flowFrame. This contains information about the type, the dimensions
## and the name of the image file. During initiation of the class we make
## sure that both bitmap and vector versions of the image are present and
## that the files are copied to the correct location. Object of the class
## should be created using the constructor
## ---------------------------------------------------------------------------
setClass("qaGraph",
         representation(fileNames="character",
                        dimensions="matrix",
                        types="character",
                        id="character"))

## ImageMagick is very particular about how to specify a path. This function
## tries to accomodate all possible inputs (realtive and absolute)
constructPath <- function(path){
    curDir <- getwd()
    if(path == "")
        return(curDir)
    curDirSplt <- strsplit(curDir, "[/\\]")[[1]]
    pathSplt <-  strsplit(path, "[/\\]")[[1]]
    if(.Platform$OS.type=="windows"){
		if(curDirSplt[[1]] !=  pathSplt[[1]]){
	        #curDir
			paste( c(curDir,pathSplt),collapse=.Platform$file.sep, sep="")
        }else{
			common <- suppressWarnings(min(which(curDirSplt[-1] != pathSplt[-1])))
			if(is.infinite(common)){
				curDir
			}else if(common == length(pathSplt)){
	    # else
	    # {
		# if(common==1)
		# {		path
        #    	   paste(c(curDir, pathSplt[-1]), 
 		#        collapse=.Platform$file.sep, sep="")
		# }
		
		
		    path
			}else{
				file.path(paste(curDirSplt[1:(common)], 
	               collapse=.Platform$file.sep, sep=""), 
		       paste(pathSplt[(common+1):length(pathSplt)],
		       collapse=.Platform$file.sep, sep=""))  
			}
	    } 
	
	      
    }else{	
        common <- suppressWarnings(min(which(curDirSplt != pathSplt)))
		if(is.infinite(common)){
			curDir
		}else{
       	    if(common==1)
			{
                file.path(curDir, path)
			} else if((common+1) == length(pathSplt)){
				path
			}else{
                file.path(paste(curDirSplt[1:(common-1)], 
	                  collapse=.Platform$file.sep, sep=""), 
		       	  paste(pathSplt[common:length(pathSplt)],
		       	  collapse=.Platform$file.sep, sep=""))
			}
		}
    }
}


win2UnixPath <- function(path)
{
    gsub("\\", "/", path.expand(path), fixed=TRUE)
}

setMethod("initialize", "qaGraph",
          function(.Object, fileName, imageDir, width=NULL, empty=FALSE, pdf=TRUE){
              if(!empty){
	          fileName <- win2UnixPath(fileName)
	          imageDir <- win2UnixPath(imageDir)
                  ## check arguments
                  if(!is.null(width) && (length(width)!=1 ||
                                         !is.numeric(width)))
                      stop("'width' must be numeric scalar")
                  if(!file.exists(fileName))
                      stop("Unable to find file '", fileName, "'")
          
                  ## get file information
                  if(!file.exists(imageDir))
                      dir.create(imageDir, recursive=TRUE)
		  fn <- basename(fileName)
		  file.copy(fileName, imageDir)
		  inf <- basename(system(paste("identify", shQuote(fileName)),
                                  intern=TRUE))
                  iInf <- gsub(paste(".*", fn, " ", sep=""), "", inf)
                  imageInfo <- c(fileName, strsplit(iInf, " ")[[1]][1:2])
                  names(imageInfo) <- c("file", "type", "dimensions")
                  bname <- basename(gsub("\\..*$", "", fileName))
                  dims <- as.numeric(strsplit(imageInfo["dimensions"],"x")[[1]])
                  newDims <- dims
                  if(!is.null(width)){
                      scaleFac <- width/dims[1]
                      if(scaleFac!=1)
                          newDims <-  dims*scaleFac
                  }
                  ft <- c("vectorized", "bitmap")
                  cf <- file.path(imageDir, fn)
                  
                  ## convert image to vectorized or bitmap version
          if(tolower(imageInfo["type"])=="pdf")
		  {
                      ## original image is vectorized
                      convType <- "jpg"
                      newFileName <- file.path(imageDir, paste(bname, convType,
		      		     	        sep="."))
                      if(!file.exists(newFileName))
                        system(paste("convert"
                                ,shQuote(fileName)
                                ,"-channel RGBA -separate -resize ", paste(newDims,collapse="x")
                                , "-combine "
                                , shQuote(newFileName)
                            )
                        )
                            
                      type <- c(ifelse(pdf, "pdf", NA), "jpg")
					#  files <-c(paste( constructPath(cf),fileName,sep="/"),
					#           paste(constructPath(cf),newFileName,sep="/"))
                     files <- c(constructPath(fileName), constructPath(newFileName))
					  if(file.exists(fileName) && pdf==FALSE){
						file.remove(fileName)
						files[1] <- NA
					  }
							   
           }else{
                      ## original image is bitmap
                      convType <- "pdf"
                      newFileName <- file.path(imageDir, paste(bname, convType,
                                                               sep="."))
                      if(!file.exists(newFileName) && pdf)
                          system(paste("convert", fileName, newFileName))
                      if(!file.exists(cf))
                          system(paste("convert -resize", paste(newDims,
                                                                collapse="x"),
                                       shQuote(fileName), shQuote(cf)))
                      type <- c(ifelse(pdf, "pdf", NA), tolower(imageInfo["type"]))
                      #files <- c(newFileName, constructPath(cf))
					  files <- c(constructPath(fileName), constructPath(newFileName))
					  # files <-c(paste( constructPath(cf),fileName,sep="/"),
					           # paste(constructPath(cf),newFileName,sep="/"))
                  }
                  ## fill qaGraph object
                  .Object@dimensions <-  matrix(c(dims, newDims), ncol=2,
                                                byrow=TRUE,
                                                dimnames=list(ft, c("width",
                                                "height")))
                  names(type) <- names(files) <- ft
                  .Object@types <- type
                  .Object@fileNames <- files
                  .Object@id=guid()
              }
              return(.Object)
          })


qaGraph <- function(...)
    new("qaGraph", ...)

 

## ===========================================================================
## qaGraphList
## ---------------------------------------------------------------------------
## A list of qaGraphs. Elements of the list are individual qaGraphs. This
## mainly exists for the sake of a constructor to facilitate object creation
## ---------------------------------------------------------------------------
setClass("qaGraphList",
         contains="list")

setMethod("initialize", "qaGraphList",
          function(.Object, imageFiles, imageDir, width=NULL, pdf=TRUE) {
              if(!all(file.exists(imageFiles)))
                  stop("'imageFiles' must be character vector of ",
                       "paths to image files")
              
              input <- lapply(imageFiles, qaGraph, imageDir=imageDir,
                              width=width, pdf=pdf)
              .Object@.Data=input
              return(.Object)
          })

qaGraphList <- function(...)
    new("qaGraphList", ...)



## ===========================================================================
## qaProcessFrame
## ---------------------------------------------------------------------------
## A class that bundles all information about the QA output for a single
## flowFrame. Slots of this class are:
##   id: a unique identifier for the whole block
##   frameID: the id of the respective flowFrame
##   summaryAggregator: a binaryAggregator indicating status
##   summaryGraph: an optional image linked to the summaryAggregator
##   frameAggregators: a list of aggregators for this frame
##   frameGraphs: a list of optional images linked to the respective
##                frameaggregators
## ---------------------------------------------------------------------------
setClass("qaProcessFrame",
         representation(id="character",
                        frameID="character",
                        summaryAggregator="qaAggregator",
                        summaryGraph="qaGraph",
                        frameAggregators="aggregatorList",
                        frameGraphs="qaGraphList",
                        details="list"))



setMethod("initialize", "qaProcessFrame",
          function(.Object, frameID, summaryAggregator, summaryGraph,
                   frameAggregators, frameGraphs, details){
              .Object@id <- guid()
              .Object@frameID <- frameID
              .Object@summaryAggregator <- summaryAggregator
              if(missing(summaryGraph))
                  summaryGraph <- new("qaGraph", empty=TRUE)
              .Object@summaryGraph <- summaryGraph
              if(xor(missing(frameAggregators), missing(frameGraphs)))
                 stop("Both 'frameAggregators' and 'frameGraphs' must be ",
                       "specified")
              else if(!missing(frameAggregators)){
                  .Object@frameAggregators <- frameAggregators
                  .Object@frameGraphs <- frameGraphs
              }
              if(!missing(details))
                  .Object@details <- details
              return(.Object)
          })

qaProcessFrame <- function(...)
    new("qaProcessFrame", ...)



## ===========================================================================
## qaProcess
## ---------------------------------------------------------------------------
## A class that describes a QA process and the associated graphical output
## that is generated during this process. Objects of this class encapsulate
## all the necessary information to generate the HTML output including all
## links that is bundled in a QA report as generated by function
## 'qaReport'
## ---------------------------------------------------------------------------
setClass("qaProcess",
         representation(id="character",
                        name="character",
                        type="character",
                        frameIDs="character",
                        summaryGraph="qaGraph",
                        frameProcesses="list"))
                        


setMethod("initialize", "qaProcess",
          function(.Object, id, name, type, summaryGraph,
                   frameProcesses) {
              names(frameProcesses) <- sapply(frameProcesses, slot, "frameID")
              fids <- names(frameProcesses)
              if(any(duplicated(fids)))
                  stop("IDs of the frame processes are not unique",
                       call.=FALSE)
              if(missing(name))
                  name <- "anonymous QA Process"
              if(missing(summaryGraph))
                  summaryGraph <- qaGraph(empty=TRUE)
		      .Object@id <- id
              .Object@name <- name
              .Object@type <- type
              .Object@frameIDs <- fids
              .Object@summaryGraph <- summaryGraph
              .Object@frameProcesses <- frameProcesses
			  validProcess(.Object)
              return(.Object)
          })

qaProcess <- function(id, name, type, summaryGraph, frameProcesses)
{
    if(missing(name))
        name <- "anonymous QA Process"
    if(missing(summaryGraph))
        summaryGraph <- qaGraph(empty=TRUE)
    new("qaProcess", id=id, name=name, type=type,
        summaryGraph=summaryGraph, frameProcesses=frameProcesses)
}



## ===========================================================================
## qaProcessSummary
## ---------------------------------------------------------------------------
## A numeric matrix of summary values over all qaProcesses, split
## by fluorochromes and samples. This mainly exists to allow method
## dispatch
## ---------------------------------------------------------------------------
setClass("qaProcessSummary",
         representation(panels="list",
                        summary="matrix",
                        ranges="list",
                        mapping="list",
                        pnams="character",
			overallSum="matrix")
         )
