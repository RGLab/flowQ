\name{qaProcess-class}
\docType{class}
\alias{qaProcess-class}
\alias{qaProcess}
\alias{initialize,qaProcess-method}
\alias{show,qaProcess-method}

\title{ Abstraction of the results of a QA process }

\description{ QA processes create graphical output which can be bundled
  in a single HTML document. This class stores all information that is
  needed by \code{\link[flowQ]{writeQAReport}} to produce such HTML reports.}

\section{Objects from the Class}{ Objects should be created using the
  constructor functions. See \code{\link[flowQ]{qaProcess.timeline}} and
  \code{\link[flowQ]{qaProcess.marginevents}} for details. When writing new QA
  process functions, the constructors
  \code{\link[flowQ:qaProcessFrame-class]{qaProcessFrame}} and 
  \code{qaProcess} should be used. The latter expects the mandatory
  arguments \code{id}, \code{type} and \code{frameProcesses} and also
  accepts the optional arguments \code{name} and \code{summaryGraph}.
  See the vignette of this package for details.}

\section{Slots}{
	 \describe{
    \item{\code{id}:}{Object of class \code{"character"}, the objects
      unique identifier }
    \item{\code{name}:}{Object of class \code{"character"}, the name of
      the process }
    \item{\code{type}:}{Object of class \code{"character"}, the type of
	 process } 
    \item{\code{frameIDs}:}{Object of class \code{"character"}, the
      identifiers of the \code{\link[flowCore:flowSet-class]{flowSets}}
	 to which the subprocesses are linkes }
    \item{\code{summaryGraph}:}{Object of class \code{"qaGraph"}, a
      graphical summary of the processe's outcome }
    \item{\code{frameProcesses}:}{Object of class \code{"list"}, more
      detailed information for each
	 \code{\link[flowCore:flowFrame-class]{flowFrame}} } 
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "qaProcess")}: constructor }
  }
  
}

\author{ Florian Hahne}

\seealso{
  \code{\link[flowQ:qaGraphList-class]{qaGraphList}},
  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ:qaProcessFrame-class]{qaProcessFrame}}
}

\examples{
showClass("qaProcess")
}
\keyword{classes}
