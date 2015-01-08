\name{EnhancedTxDb-class}
\docType{class}

% Class:
\alias{class:EnhancedTxDb}
\alias{EnhancedTxDb-class}

% Constructors:
\alias{EnhancedTxDb}
\alias{EnhancedTxDb,TranscriptDb-method}

% Accessors:
\alias{transcriptsBy}
\alias{transcriptsBy,EnhancedTxDb-method}
\alias{transcripts}
\alias{transcripts,EnhancedTxDb-method}
\alias{exonsBy}
\alias{exonsBy,EnhancedTxDb-method}
\alias{exons}
\alias{exons,EnhancedTxDb-method}
\alias{cdsBy}
\alias{cdsBy,EnhancedTxDb-method}
\alias{cds}
\alias{cds,EnhancedTxDb-method}

% "show" method:
\alias{show,EnhancedTxDb}

\title{EnhancedTxDb objects}

\description{
  The EnhancedTxDb class is a container for a
  \code{\linkS4class{TranscriptDb}} object and the associated annotation
  information.
}

\details{
  A EnhancedTxDb object i
  
}

\section{Constructors}{
  \describe{
    \item{}{
      \code{EnhancedTxDb(txdb, gnAnnot, txAnnont, exAnnot,...)}:
      Constructor to create a EnhancedTxDb object based on splice
      graph edge counts, the associated splice graph structure and the
      \describe{
        \item{\code{txdb}}{\code{\linkS4class{TranscriptDb}}}
        \item{\code{gnAnnot}}{\code{data.frame}}
        \item{\code{txAnnont}}{\code{data.frame}}
	\item{\code{exAnnont}}{\code{data.frame}}	
      }
    }
  }
}

\section{Accessors}{
  In the code below, \code{x} refers to a EnhancedTxDb object.
  \describe{
    \item{}{
      \code{transcriptsBy(x)}:
      Returns transcripts grouped by supplied elements
    }
    \item{}{
      \code{transcripts(x)}:
      Returns all transcripts
    }
    \item{}{
      \code{exonsBy(x)}:
      Returns all exons
    }
    \item{}{
      \code{exons(x)}:
      Returns all exons
    }
  }
}

\author{
  D. Bindreither
}

\seealso{
  \itemize{
    \item \link{TranscriptDb-class}
  }
}

\examples{
\dontrun{
etxdb <- EnhancedTxDb()
}
}
