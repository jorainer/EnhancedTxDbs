\name{loadEnhancedTxDb}

\alias{loadEnhancedTxDb}

\title{
  Loading \code{\linkS4class{EnhancedTxDb}} objects
}

\description{
  Takes a file name and loads a previously saved
  \code{\linkS4class{EnhancedTxDb}} object.
}

\usage{
loadEnhancedTxDb(fn, ...)
}

\arguments{
  \item{fn}{
    File name without extension used for saving the object
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
}

\details{
  Taking a file name and loading a previously saved
  \code{\linkS4class{EnhancedTxDb}} object. Only the name of the file,
  without file extensions should be provided.
}

\author{
  D. Bindreither
}

\seealso{
  \link{TranscriptDb}
}

\examples{
\dontrun{
  path <- system.file("extdata", package="EnhancedTxDbs")
  fn <- paste(path, "ens75chr19", sep="/")
  EnhancedTxDb <- loadEnhancedTxDb(fn)
}
}
