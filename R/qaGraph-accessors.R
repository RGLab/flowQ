## ===========================================================================
## qaGraph
## ---------------------------------------------------------------------------

## display details about image files
setMethod("show", signature("qaGraph"),
          function(object)
              cat("QA process image information\n")
          )

## return file names
setMethod("names", signature("qaGraph"),
          function(x){
          if(length(x@fileNames)==0)
              return(NULL)
          basename(x@fileNames["bitmap"])
          })



## ===========================================================================
## qaGraphList
## ---------------------------------------------------------------------------

## display details about list
setMethod("show", signature("qaGraphList"),
          function(object)
              cat("List of", length(object), "images\n")
          )
