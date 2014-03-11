\name{discreteAggregator-class}
\docType{class}
\alias{discreteAggregator-class}
\alias{discreteAggregator}
\alias{show,discreteAggregator-method}
\alias{writeLines,discreteAggregator,file,missing,missing-method}

\title{Class "discreteAggregator" }

\description{ Abstraction of a discrete type of aggregator with possible
states "passed", "not passed" and "warning"}

\section{Objects from the Class}{ Objects can be created by calls of the
form \code{new("discreteAggregator", ...)} or using the constructor
\code{discreteAggregator(x)}, where \code{x} is an integer value in
\code{[0,1,2]} or a factor with levels \code{0}, \code{1} and \code{2}.
}

\section{Slots}{
	 \describe{
    \item{\code{x}:}{Object of class \code{"factor"} One in \code{0}
	 (not passes),
	 \code{1} (passed) or \code{2} (warning)}
    \item{\code{passed}:}{Object of class \code{"logical"} indicating
      whether the process has passed the QA requirements (not eveluated) }
  }
}
\section{Extends}{
Class \code{\link[flowQ:qaAggregator-class]{qaAggregator}}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "discreteAggregator")}: print
  object details }
\item{writeLines}{\code{signature(text = "discreteAggregator",
    con = "file", sep = "missing", useBytes="missing")}:  write to HTML
	 file connection } 
	 }
 }

 
\author{ Florian Hahne }


\seealso{
  \code{\link[flowQ]{qaProcess.marginevents}},
  \code{\link[flowQ]{qaReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}},
  \code{\link[flowQ]{qaProcess.timeline}},
  \code{\link[flowQ:binaryAggregator-class]{binaryAggregator}},
  \code{\link[flowQ:factorAggregator-class]{factorAggregator}},
  \code{\link[flowQ:numericAggregator-class]{numericAggregator}},
  \code{\link[flowQ:stringAggregator-class]{stringAggregator}},
  \code{\link[flowQ:rangeAggregator-class]{rangeAggregator}}
}
\examples{
showClass("discreteAggregator")
}
\keyword{classes}
 
