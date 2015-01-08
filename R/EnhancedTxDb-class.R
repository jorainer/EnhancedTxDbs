### ==================================================================
### EnhancedTxDb object
### ------------------------------------------------------------------
### Stores TranscriptDb objects together with annotation information from ENSEMBL

setClass("EnhancedTxDb",
         representation(txdb="TxDb",
                        txAnnot="data.frame",
                        gnAnnot="data.frame",
                        exAnnot="data.frame"),
         prototype(txdb=NULL,
                   txAnnot=NULL,
                   gnAnnot=NULL,
                   exAnnot=NULL)
         )

## Validity method
.valid.EnhancedTxDb <- function(x)
{
     if(! all(names(exonsBy(x@txdb, "tx")) %in% rownames(x@txAnnot))){
       msg <- paste("Transcript IDs in TranscriptDb object have to match",
                    "transcript IDs in the transcript table\n")
     }
}

setValidity2("EnhancedTxDb", .valid.EnhancedTxDb)

## Show method
setMethod("show", "EnhancedTxDb",
function(object){
  cat("EnhancedTxDb object\n")
  metadata <- metadata(object@txdb)
  cat(paste("EnhancedTxDb:", installed.packages()["EnhancedTxDbs","Version"], "\n"))
  cat(apply(metadata(object@txdb)[c(-1,-2),],1, function(x)
            paste(x[1], ": ", x[2], "\n", collapse="", sep="")))
})

## Constructor
setGeneric("EnhancedTxDb",
    function(txdb, txAnnot, gnAnnot, exAnnot, ...) standardGeneric("EnhancedTxDb"))

setMethod(EnhancedTxDb, "TranscriptDb",
    function(txdb, txAnnot, gnAnnot, exAnnot, ..., verbose=FALSE)
{
    new("EnhancedTxDb", txdb=txdb, txAnnot=txAnnot, gnAnnot=gnAnnot, exAnnot=exAnnot, ...)
})

setMethod(EnhancedTxDb, "TxDb",
    function(txdb, txAnnot, gnAnnot, exAnnot, ..., verbose=FALSE)
{
    new("EnhancedTxDb", txdb=txdb, txAnnot=txAnnot, gnAnnot=gnAnnot, exAnnot=exAnnot, ...)
})

## Accessor for transcripts
setMethod("transcriptsBy", "EnhancedTxDb",
function(x, ...)
{
  by <- match.arg(by)
  txs <- transcriptsBy(x@txdb, ...)

  gn <- rep(names(txs), elementLengths(txs))
  txs <- unlist(txs, use.names=FALSE)
    df <- DataFrame(values(txs), x@txAnnot[values(txs)[["tx_id"]], ])
  colnames(df) <- c(colnames(values(txs)),
                      colnames(x@txAnnot[values(txs)[["tx_id"]], , drop=FALSE]))
    values(txs) <- df
  split(txs, gn)
})

setMethod("transcripts", "EnhancedTxDb",
function(x, ...)
{
  txs <- transcripts(x@txdb, ...)
  df <- DataFrame(values(txs), x@txAnnot[values(txs)[["tx_id"]], ])
  colnames(df) <- c(colnames(values(txs)),
                    colnames(x@txAnnot[values(txs)[["tx_id"]], , drop=FALSE]))
  values(txs) <- df
  txs
})

## Accessors for CDS
setMethod("cds", "EnhancedTxDb",
function(x, ...)
{
  cds(x@txdb, ...)
})

setMethod("cdsBy", "EnhancedTxDb",
function(x, ...)
{
  cdsBy(x@txdb, ...)
})


## Accesors for exons
setMethod("exons", "EnhancedTxDb",
function(x, ...)
{
  exons(x@txdb, ...)
})

setMethod("exonsBy", "EnhancedTxDb",
function(x, ...)
{
  exonsBy(x@txdb, ...)
})

setReplaceMethod("isActiveSeq","EnhancedTxDb",
function(x, value){
  GenomicFeatures::isActiveSeq(x@txdb) <- value
  x
})

## Set active seqisActiveSeq<-"
setMethod("isActiveSeq", "EnhancedTxDb", function(x){isActiveSeq(x@txdb)})
setMethod("seqinfo", "EnhancedTxDb", function(x){seqinfo(x@txdb)})
