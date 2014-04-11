\name{qaProcessSummary}
\docType{class}
\alias{qaProcessSummary-class}
\alias{qaProcessSummary}
\alias{writeLines,qaProcessSummary,file,missing,missing-method}
\alias{show,qaProcessSummary-method}

\title{Class "qaProcessSummary" }
\description{ An internal class to represent QA summaries}
\section{Objects from the Class}{ The class is internal and not ment for
  interactive use. 

}
\section{Slots}{
  
  \describe{
    
    \item{\code{panels}:}{Object of class \code{"list"}}
    
    \item{\code{ranges}:}{Object of class \code{"matrix"}}
    
    \item{\code{summary}:}{Object of class \code{"list"}}
    
    \item{\code{mapping}:}{Object of class \code{"list"}}

    \item{\code{pnams}:}{Object of class \code{"character"}}

    \item{\code{overallSum}:}{Object of class \code{"matrix"}} 
	      
  }
}

\section{Methods}{
  \describe{
    \item{writeLines}{\code{signature(text = "qaGraphList", con =
  "file")}: HTML output }
    \item{show}{\code{signature(object = "qaGraphList")}: print object details }
  }
}

    
\author{ Florian Hahne }


\keyword{classes}
  