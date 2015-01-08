\name{saveEnhancedTxDb}

\alias{saveEnhancedTxDb}
\alias{saveEnhancedTxDb,EnhancedTxDb,character-method}

\title{
  Serializing \code{\linkS4class{EnhancedTxDb}} objects
}

\description{
  Takes a \code{\linkS4class{EnhancedTxDb}} object and saves it to the
  given file name.
}

\usage{
saveEnhancedTxDb(eTxDb, fn, ...)
\S4method{saveEnhancedTxDb}{EnhancedTxDb,character}(eTxDb, fn, ...)
}

\arguments{  
  \item{eTxDb}{
    \code{\linkS4class{EnhancedTxDb}} object which should be serialized
  }
  \item{fn}{
    File name without extension used for saving the object
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
}

\details{
  Taking a \code{\linkS4class{EnhancedTxDb}} and a file name to which
  the object should be saved. Only the name of the file, without file
  extensions should be provided.
}

\author{
  D. Bindreither
}

\seealso{
  \link{TranscriptDb}
}

\examples{
\dontrun{
  ## TranscriptDb-object
  path <- system.file("extdata", package="EnhancedTxDbs")
  fn <- paste(path, "ens67chr19", sep="/")
  EnhancedTxDb <- loadEnhancedTxDb(fn)
  saveEnhancedTxDb(EnhancedTxDb, "myEnhancedTxDb")
}
}