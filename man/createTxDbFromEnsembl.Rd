\name{createTxDbFromEnsembl}

\alias{createTxDbFromEnsembl}

\title{
  Creates an EnhancedTxDb object from the ENSEMBL data base.
}

\description{
  Takes a data base connection to create a EnhancedTxDb based on
  ENSEMBL. In contrast to the normal \code{\linkS4class{TranscriptDb}}
  objects EnhancedTxDb objects allow to store also annotation information.
  At the moment the only used annotation information is the ENSEMBL biotype.
}

\usage{
createTxDbFromEnsembl(version, annotInfo=TRUE,
ensemblApi="/home/bioinfo/ensembl/73/",
bioperl="~/ENSEMBL/bioperl-1.2.3",
species="human", user="anonuser", host="madb.i-med.ac.at", 
pass="", chr=character(), ...)
}

\arguments{  
  \item{version}{
    Ensembl version from which the object should be created. E.g.: 67
  }
  \item{annotInfo}{
    Logical. If \code{TRUE} an EnhancedTxDb object is
    returned. Otherwise a \code{\linkS4class{TranscriptDb}} object is
    created.
  }
  \item{ensemblApi}{
    Data base from which the information should be retrieved. Default is
    \code{core}.
  }
  \item{bioperl}{
    Location of the BioPerl library.
  }
  
  \item{species}{
    Species for which the \code{\linkS4class{TranscriptDb}} or
    EnhancedTxDb object should be created. Default is \code{human}.
  }
  \item{user}{
    Username for the data base. Default is \code{anonuser}.
  }
  \item{host}{
    Sever where the data base is running.
  }
  \item{pass}{
    Password of the chosen data base user. Default no password is set.
  }
  \item{chr}{
    Utilized chromosomes. E.g.: chr=\code{c("X", "19")}
  }
  \item{...}{
    Arguments to be passed to or from methods.
  }
}

\details{
  Making a \code{\linkS4class{TranscriptDb}} or a EnhancedTxDb object
  from ENSEMBL containing annotation information from ENSEMBL.
  The function utilizes a perl script to retrieve the required
  information from ENSEMBL. To make the function work it is essential to
  have the ENSEMBL Perl API installed on your machine and that the
  provided data base connection is correct.
}

\value{ 
  Returns either a \code{\linkS4class{EnhancedTxDb}} object containing
  the coordinates and annotation of the provided ENSEMBL version or a
  \code{\linkS4class{TranscriptDb}} providing only the genomic
  coordinates of the transcripts, exons, cds and genes.
}

\author{
  D. Bindreither
}

\seealso{
  \link{TranscriptDb}, \link{GRangesList}
}

\examples{
\dontrun{
  ## TranscriptDb-object
  txdb <- createTxDbFromEnsembl(67, annotInfo=FALSE, chr=19)
  txdb

  ## EnhancedTxDb-object
  etxdb <- createTxDbFromEnsembl(67, chr=19)
  etxdb
}
}