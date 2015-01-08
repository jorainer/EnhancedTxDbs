## Methods for saving and loading EnhancedTxDb objects

## saving
setMethod("saveEnhancedTxDb", c("EnhancedTxDb", "character"),
          function(eTxDb, fn, ...)
{
  txdb <- eTxDb@txdb
  txAnnot <- eTxDb@txAnnot
  gnAnnot <- eTxDb@gnAnnot
  exAnnot <- eTxDb@exAnnot
  save(txAnnot, gnAnnot, exAnnot, file=paste(fn, "Rda", sep="."))
  saveDb(txdb, file=paste(fn, "sqlite", sep="."))
})


## loading
loadEnhancedTxDb <- function(fn, ...)
{
  txAnnot <- gnAnnot <- exAnnot <- NULL
  load(paste(fn, "Rda", sep="."))
  txdb <- loadDb(paste(fn, "sqlite", sep="."))
  EnhancedTxDb(txdb=txdb, txAnnot=txAnnot,
               gnAnnot=gnAnnot, exAnnot=exAnnot)
}
