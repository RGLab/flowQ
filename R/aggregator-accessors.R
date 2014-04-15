## ===========================================================================
## binaryAggregator
## ---------------------------------------------------------------------------

## write method to create HTML output
setMethod("writeLines", signature("binaryAggregator", "file", "missing","missing"),
          function(text, con){
              if(!is.na(text@passed)){
                  if(text@passed)
                      writeLines(paste("<img class=\"QABinAggr\"",
                                       " src=\"images/bulbGreen.png\">"), con)
                  else
                      writeLines(paste("<img class=\"QABinAggr\"",
                                       " src=\"images/bulbRed.png\">"), con)
              }
          })

## display details about aggregator
setMethod("show", signature("binaryAggregator"),
          function(object){
              cat("Binary quality score ", ifelse(object@passed, "", "not "),
                  "passing the requirements\n", sep="") 
          })



## ===========================================================================
## discreteAggregator
## ---------------------------------------------------------------------------

## write method to create HTML output
setMethod("writeLines", signature("discreteAggregator", "file", "missing","missing"),
          function(text, con){
              switch(as.character(text@x),
                     "1"=writeLines(paste("<img class=\"QABinAggr\"",
                     " src=\"images/bulbGreen.png\">"), con),
                     "0"=writeLines(paste("<img class=\"QABinAggr\"",
                     " src=\"images/bulbRed.png\">"), con),
                     "2"=writeLines(paste("<img class=\"QABinAggr\"",
                     " src=\"images/bulbYellow.png\">"), con),
                     stop("Unknown state"))
          })

## display details about aggregator
setMethod("show", signature("discreteAggregator"),
          function(object){
              states <- c("failed", "passed", "warn")
                cat("Discrete quality score ",
                    ifelse(object@passed, "", "not "), "passing the ",
                    "requirements with state ",
                    states[as.integer(object@x)], "\n",
                    sep="") 
          })



## ===========================================================================
## factorAggregator
## ---------------------------------------------------------------------------

## write method to create HTML output
setMethod("writeLines", signature("factorAggregator", "file", "missing","missing"),
          function(text, con){
              col <- ifelse(text@passed, "green", "red")
              lx <- levels(text@x)
              fcol <- rep("lightgray", length(lx))
              fcol[match(text@x, lx)] <- col
              writeLines("<b>", con)
              writeLines(paste("<span class=\"QAFactAggr\" ",
                               "style=\"color:", fcol, ";\">", lx,
                               "</span>", sep=""), con)
              writeLines("</b><br>", con)
          })


## display details about aggregator
setMethod("show", signature("factorAggregator"),
          function(object){
              cat("Factorized quality score ", ifelse(object@passed, "",
                                                      "not "),
                  "passing the requirements of value=", as.character(object@x),
                  "\n", sep="") 
          })



## ===========================================================================
## stringAggregator
## ---------------------------------------------------------------------------

## write method to create HTML output
setMethod("writeLines", signature("stringAggregator", "file", "missing","missing"),
          function(text, con){
              col <- ifelse(text@passed, "green", "red")
              writeLines(paste("<b><span class=\"QAStringAggr\" ",
                               "style=\"color:", col, ";\">", text@x,
                               "</span></b>", sep=""), con)
          })

## display details about aggregator
setMethod("show", signature("stringAggregator"),
          function(object){
              cat("Textual quality score ", ifelse(object@passed, "", "not "),
                  "passing the requirements of value=", object@x,
                  "\n", sep="") 
          })


## ===========================================================================
## numericAggregator
## ---------------------------------------------------------------------------

## write method to create HTML output
setMethod("writeLines", signature("numericAggregator", "file", "missing","missing"),
          function(text, con){
              col <- ifelse(text@passed, "green", "red")
              writeLines(paste("<div class=\"QANumAggr\" ",
                               "style=\"color:", col, ";\">", signif(text@x,2),
                               "</div>", sep=""), con)
          })

## display details about aggregator
setMethod("show", signature("numericAggregator"),
          function(object){
              cat("Numeric quality score ", ifelse(object@passed, "", "not "),
                  "passing the requirements of value=", object@x,
                  "\n", sep="") 
          })


## ===========================================================================
## rangeAggregator
## ---------------------------------------------------------------------------

## write method to create HTML output
setMethod("writeLines", signature("rangeAggregator", "file", "missing","missing"),
          function(text, con){
              x <- text
              perc <- (x@x-x@min)/diff(c(x@min, x@max))*100
              class <- ifelse(x@passed, "passed", "failed")
              writeLines(paste("<div class=\"QARange\">\n",
                               "<div class=\"", class, "\" style=\"width:", perc,
                               "%\"></div>\n",
                               "<div class=\"label\">", signif(perc,1), "%</div>\n",
                               "</div>", sep=""), con)
          })

## display details about aggregator
setMethod("show", signature("rangeAggregator"),
          function(object){
              cat("Range quality score ", ifelse(object@passed, "", "not "),
                  "passing the requirements of value=", object@x,
                  "\n", sep="") 
          })



## ===========================================================================
## aggregatorList
## ---------------------------------------------------------------------------

## display details about list
setMethod("show", signature("aggregatorList"),
          function(object)
              cat("List of", length(object), "aggregators\n")
          )

