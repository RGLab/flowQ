\name{qaGraph-class}
\docType{class}
\alias{qaGraph-class}
\alias{qaGraph}
\alias{initialize,qaGraph-method}
\alias{names,qaGraph-method}
\alias{show,qaGraph-method}

\title{Class "qaGraph"}
\description{ Abstraction of the graphical output created in the cause
  of a QA process }
\section{Objects from the Class}{
  Objects should be created using the constructor:

  \code{qaGraph(fileName, imageDir, width, empty=FALSE, pdf=TRUE)},
  where \code{fileName} is a path to an image file, \code{imageDir} is
  the destination path for the images and the optional argument
  \code{width} is the final width to which the bitmap images are
  converted. The \code{pdf} argument controls wether vectorized versions
  of the image are to be produced. For the special case of an empty
  object without any images one can use option \code{empty=TRUE}, in
  which case the constructor ignores all other arguments.

  During object instantiation the image file will be converted, resized
  and copied if necessary.
}

\section{Slots}{
	 \describe{
    \item{\code{fileNames}:}{Object of class \code{"character"} The paths
    to the image files, both the vectorized and unverctorized versions}
    \item{\code{dimensions}:}{Object of class \code{"matrix"} The
      dimensions of the image files, both for the vectorized and
      unvectorized version}
    \item{\code{types}:}{Object of class \code{"character"} The file
      extensions for both versions }
    \item{\code{id}:}{Object of class \code{"character"} A unique
      identifier for the images }
  }
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "qaGraph")}: constructor }
    \item{names}{\code{signature(x = "qaGraph")}: returns the file name
      of the bitmap version of the image}
    \item{show}{\code{signature(object = "qaGraph")}: print object details }
	 }
}

\author{ Florian Hahne }


\seealso{ \code{\link[flowQ:qaGraphList-class]{qaGraphList}},
  \code{\link[flowQ]{writeQAReport}},
  \code{\link[flowQ:qaProcess-class]{qaProcess}}
}
\examples{
showClass("qaGraph")
}
\keyword{classes}
