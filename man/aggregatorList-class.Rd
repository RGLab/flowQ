\name{aggregatorList-class}
\docType{class}
\alias{aggregatorList-class}
\alias{aggregatorList}
\alias{initialize,aggregatorList-method}
\alias{show,aggregatorList-method}

\title{Class "aggregatorList" }
\description{ A list of qaAggregators }
\section{Objects from the Class}{
  Objects should be created using the constructor:

  \code{aggregatorList(\dots)}, where \code{\dots} are objects
  inheriting form class
  \code{\link[flowQ:qaAggregator-class]{qaAggregator}} or a list of such
  objects. 
}

\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"list"}, the list data }
  }
}
\section{Extends}{
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "aggregatorList")}:
      constructor }
    
    \item{show}{\code{signature(object = "aggregatorList")}: print
      object details }
  }
}
       
\details{This class directly extends class \code{"list"} and is intended
  to exclusively hold objects of class
  \code{qaAggregator}, where each list item represents the
  outcome of a QA subprocess for a single
  \code{\link[flowCore:flowFrame-class]{flowFrame}}. It mainly exists to
  allow for method dispatch and should never be populated manualy;
  instead, use the constructor \code{qaAggregatorList} which checks for
  valid objects.} 

\author{ Florian Hahne }


\seealso{
  \code{\link[flowQ:qaGraph-class]{qaGraph}},
  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}},
  \code{\link[flowQ:qaAggregator-class]{qaAggregator}}
}
\examples{
showClass("aggregatorList")
}
\keyword{classes}
  
