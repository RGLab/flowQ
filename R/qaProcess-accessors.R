## ===========================================================================
## qaProcessFrame
## ---------------------------------------------------------------------------

## display details about process
setMethod("show", signature("qaProcessFrame"),
          function(object)
              cat("Quality process output for frame '", object@frameID, "' ",
                  ifelse(object@summaryAggregator@passed, "", "not "),
                  "passing the requirements\n", sep="")
          )



## ===========================================================================
## qaProcess
## ---------------------------------------------------------------------------

## display details about process
setMethod("show", signature("qaProcess"),
          function(object)
              cat("Quality process '", object@name, "' of type '",
                  object@type, "'\n", sep="")
          )


## validity checking of a qaProcess object. Among the object's integrity
## this checks for the presence of all image files.
validProcess <- function(object)
{
    ## Are all mandatory slots present?
    if(is.null(object@id) || length(object@id)==0)
        stop("Slot 'id' can not be empty")
    if(is.null(object@frameIDs) || length(object@frameIDs)==0)
        stop("Slot 'frameIDs' can not be empty")
    if(is.null(object@name) || length(object@name)==0)
        stop("Slot 'name' can not be empty")
    if(is.null(object@frameProcesses) || length(object@frameProcesses)==0)
        stop("Slot 'frameProcesses' can not be empty")
    
   
    ## Are all image files existing?
    imageFiles <- unlist(c(sapply(object@frameProcesses, function(x)
                                  c(x@summaryGraph@fileNames["bitmap"],
                                    sapply(x@frameGraphs, function(y)
                                           y@fileNames["bitmap"]))),
                           object@summaryGraph@fileNames["bitmap"]))
    imageFiles <- imageFiles[!is.na(imageFiles)]
    missing <- !sapply(imageFiles, file.exists)
    if(all(missing))
        stop("No image files available. Check file paths")
    if(any(missing))
        stop(paste("Unable to find image file", imageFiles[missing], "\n"))
    ## Do the frameIDs match
    mismatch <- ! names(object@frameProcesses) %in% object@frameIDs
    if(any(mismatch))
       stop(paste("IDs for frame", which(mismatch), "do not match\n"))

    ## Are frameProcesses valid?
    if(!all((sapply(object@frameProcesses, function(x)
                   length(x@frameAggregators))) ==
            (sapply(object@frameProcesses, function(x)
                    length(x@frameGraphs)))))
        stop("Each qaGraph needs an associated aggregator")
    return(TRUE)
}
