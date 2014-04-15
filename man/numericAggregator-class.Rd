\name{numericAggregator-class}
\docType{class}
\alias{numericAggregator-class}
\alias{numericAggregator}
\alias{show,numericAggregator-method}
\alias{writeLines,numericAggregator,file,missing,missing-method}

\title{Class "numericAggregator" }

\description{ Abstraction of a numeric type of aggregator for which possible
  states are coded by a numeric value}

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("numericAggregator", ...)} or using the constructor
\code{numericAggregator(x, passed)}, where \code{x} is a numeric scalar,
and \code{passed} is a logical scalar.  }

\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"numeric"} coding the outcome
    state.}
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements. }
  }
}
\section{Extends}{
Class \code{\link[flowQ:qaAggregator-class]{qaAggregator}}, directly.
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "numericAggregator")}: print
  object details}
    \item{writeLines}{\code{signature(text = "numericAggregator", con = "file",
        sep = "missing", useBytes= "missing")}: write to HTML file connection }
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
  \code{\link[flowQ:binaryAggregator-class]{binaryAggregator}},
  \code{\link[flowQ:stringAggregator-class]{stringAggregator}},
  \code{\link[flowQ:rangeAggregator-class]{rangeAggregator}}
}
\examples{
showClass("numericAggregator")
}
\keyword{classes}

    
