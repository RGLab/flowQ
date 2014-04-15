\name{factorAggregator-class}
\docType{class}
\alias{factorAggregator-class}
\alias{factorAggregator}
\alias{show,factorAggregator-method}
\alias{writeLines,factorAggregator,file,missing,missing-method}

\title{Class "factorAggregator" }

\description{ Abstraction of a factor type of aggregator with possible
  states coded by the factor levels }

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("factorAggregator", ...)} or using the constructor
\code{factorAggregator(x, passed)}, where \code{x} is a factor, or an
object which can be coerced to a factor, and \code{passed} is a logical
scalar.  }

\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"factor"} coding the outcome state }
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Extends}{
Class \code{"\linkS4class{qaAggregator}"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "factorAggregator")}: print
  object details  }
    \item{writeLines}{\code{signature(text = "factorAggregator", con = "file", 
        sep = "missing", useBytes= "missing")}:   write to HTML file connection}
	 }
       }

  
\author{ Florian Hahne }


\seealso{
  \code{\link{qaProcess.marginevents}}, \code{\link{qaReport}},
  \code{\link{qaProcess}}, \code{\link{qaProcess.timeline}},
  \code{\linkS4class{discreteAggregator}},
  \code{\linkS4class{binaryAggregator}},
  \code{\linkS4class{numericAggregator}},
  \code{\linkS4class{stringAggregator}},
  \code{\linkS4class{rangeAggregator}},
}
\examples{
showClass("factorAggregator")
}
\keyword{classes}

 
