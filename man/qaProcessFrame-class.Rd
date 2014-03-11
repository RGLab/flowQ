\name{qaProcessFrame-class}
\docType{class}
\alias{qaProcessFrame-class}
\alias{qaProcessFrame}
\alias{initialize,qaProcessFrame-method}
\alias{show,qaProcessFrame-method}

\title{Class "qaProcessFrame" }
\description{ Abstraction of subitems within a qQA process}
\section{Objects from the Class}{
  Objects should be created using the constructor:

  \code{qaProcessFrame(frameID, summaryAggregator, summaryGraph,
    frameAggregators, frameGraphs, details)} where \code{frameID} is the
    ID of the flowFrame the process is linked to,
    \code{summaryAggregator} is an object inheriting from class
    \code{\link[flowQ:qaAggregator-class]{qaAggregator}} which
    summarizes the outcome, \code{summaryGraph} is an object of class
    \code{\link[flowQ:qaGraph-class]{qaGraph}} which is the overview
    graph of the process for the whole frame, \code{details} is a list
    containing any additional information regarding the QA process,
    \code{frameAggregators} is an object of class
    \code{\link[flowQ:aggregatorList-class]{aggregatorList}} and
    \code{frameGraphs} is an 
    object of class
    \code{\link[flowQ:qaGraphList-class]{qaGraphList}}. The latter two
    are the collections of aggregators and graphs for each subprocess. Only
    \code{frameID} and \code{summaryAggregator} are mandatory
    arguments.}

\section{Slots}{
	 \describe{
    \item{\code{id}:}{Object of class \code{"character"}, a unique ID of
    the object}
    \item{\code{frameID}:}{Object of class \code{"character"}, ID of the
    flowFrame the process is linked to }
    \item{\code{summaryAggregator}:}{Object of class
    \code{"qaAggregator"}, an aggregator summarizing the output of the process }
    \item{\code{summaryGraph}:}{Object of class \code{"qaGraph"}, a
    graphical summary of the process }
    \item{\code{frameAggregators}:}{Object of class
    \code{"aggregatorList"}  a list of aggregators for the subprocesses}
    \item{\code{frameGraphs}:}{Object of class \code{"qaGraphList"} a
      list of graphical summaries for the subprocesses }
    \item{\code{details}:}{A list for any additional information }
  }
}

\section{Methods}{
  \describe{
    
    \item{initialize}{\code{signature(.Object = "qaProcessFrame")}:
      constructor }
    
  }
}
       
\details{ This class bundles graphs and aggregators for a single
  \code{\link[flowCore:flowFrame-class]{flowFrame}}. This allows to
  create processes 
  with subcomponents, where each item in the \code{frameAggregators} and
  \code{frameGraphs} lists corresponds to one subprocess, which can be
  used, for instance, to create individual plots for each flow
  channel. For QA processes without subcomponents, these slots would
  simply not be populated.}

\author{ Florian Hahne }


\seealso{
  \code{\link[flowQ:qaGraphList-class]{qaGraphList}},
  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}}
}
\examples{
showClass("qaProcessFrame")
}
\keyword{classes}
