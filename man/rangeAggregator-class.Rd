\name{rangeAggregator-class}
\docType{class}
\alias{rangeAggregator-class}
\alias{rangeAggregator}
\alias{show,rangeAggregator-method}
\alias{writeLines,rangeAggregator,file,missing,missing-method}

\title{Class "rangeAggregator"  }

\description{ Abstraction of a range type of aggregator where possible
  states are within certain ranges (e.g. percentages) }

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("rangeAggregator", ...)} or using the constructor
\code{rangeAggregator(x, min, max, passed)}, where \code{x}, \code{min}
and \code{max} are numeric scalars, with \code{x} in the range of
\code{[min, max]}, and \code{passed} is a logical scalar.  }

\section{Slots}{
	 \describe{
    \item{\code{min}:}{Object of class \code{"numeric"}, the range minimum }
    \item{\code{max}:}{Object of class \code{"numeric"}, the range maximum }
    \item{\code{x}:}{Object of class \code{"numeric"}, the value within
      the range }
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements }
  }
}
\section{Extends}{
Class \code{numericAggregator}, directly.
Class \code{qaAggregator}, by class "numericAggregator",
	 distance 2. 
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "rangeAggregator")}: print
  object details }
    \item{writeLines}{\code{signature(text = "rangeAggregator", con = "file",
         sep = "missing",useBytes ="missing")}:  write to HTML file connection }
	 }
       }

\author{ Florian Hahne }


\seealso{
  \code{\link[flowQ]{qaProcess.marginevents}},
  \code{\link[flowQ]{qaReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}},
  \code{\link[flowQ]{qaProcess.timeline}},
  \code{\link[flowQ:discreteAggregator-class]{discreteAggregator}},
  \code{\link[flowQ:factorAggregator-class]{factorAggregator}},
  \code{\link[flowQ:numericAggregator-class]{numericAggregator}},
  \code{\link[flowQ:stringAggregator-class]{stringAggregator}},
  \code{\link[flowQ:binaryAggregator-class]{binaryAggregator}}
}
\examples{
showClass("rangeAggregator")
}
\keyword{classes}


 
