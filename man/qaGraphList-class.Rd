\name{qaGraphList-class}
\docType{class}
\alias{qaGraphList-class}
\alias{qaGraphList}
\alias{initialize,qaGraphList-method}
\alias{show,qaGraphList-method}

\title{Class "qaGraphList" }
\description{ A list of \code{\link[flowQ:qaGraph-class]{qaGraph}} objects}
\section{Objects from the Class}{ Objects should be created using the
  constructor:

  \code{qaGraphList(imageFiles, imageDir, width, pdf=TRUE)}, where
  \code{imageFiles} are paths to image files, imageDir is the
  destination path for the images and \code{width} is the final width to
  which the bitmap images are converted. \code{pdf} can be used to
  toggle the generation of vectorized versions of the images.
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
    \item{initialize}{\code{signature(.Object = "qaGraphList")}: constructor }
    \item{show}{\code{signature(object = "qaGraphList")}: print object details }
  }
}

\details{This class directly extends class \code{"list"} and is intended
  to exclusively hold objects of class
  \code{\link[flowQ:qaGraph-class]{qaGraph}}, 
  where each list item represents the grapical output of a QA subprocess
  for a single \code{\link[flowCore:flowFrame-class]{flowFrame}}. It
  mainly exists to allow for method dispatch and should never be
  populated manualy; instead, use the constructor \code{qaGraphList}
  which makes sure, that all image files are converted into the
  appropriate types and sizes and copied to the expected file location.}
    
\author{ Florian Hahne }


\seealso{
  \code{\link[flowQ:qaGraph-class]{qaGraph}},
  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}}
}
\examples{
showClass("qaGraphList")
}
\keyword{classes}
  
