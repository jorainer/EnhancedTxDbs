## ExTab should require all needed information for
## creating a TranscriptDb object out of it
## For creating such a object four tables are essential
## - transcripts
## - splicings
## - genes
## - chrominfo

## This function comutes the exon ranks for all genes and transcripts
computeExonRanks <- function(exTab){
  exTab <- exTab[order(exTab$chrom, exTab$gene_id, exTab$tx_id,
                       exTab$exon_start*exTab$strand,
                       exTab$exon_end*exTab$strand), ]

  exR.tmp <- seq(1:nrow(exTab))

  tx_shift <- c(exTab$tx_name[-1],
                exTab$tx_name[nrow(exTab)])

  ## mark first exon within each transcript
  bm <- exTab$tx_name != tx_shift
  luTab <- setNames(c(0, exR.tmp[bm]), c(exTab$tx_name[1], tx_shift[bm]))
  luTab <- luTab[exTab$tx_name]
  exTab$exon_rank <- exR.tmp - luTab
  exTab
}

computeExonCds <- function(exTab){
  cat("Computing exon CDS from transcript CDS\n")

  ## This function computes the cds for the individual exons
  fwdStrand <- exTab$strand == 1

  ## Compute exon cds starts
  start <-
    ifelse(exTab$cds_start < exTab$exon_start,  exTab$exon_start,
           ifelse(exTab$cds_start > exTab$exon_end,
                  rep(NA, nrow(exTab)), exTab$cds_start))

  ## Compute exon cds ends
  end <-
    ifelse(exTab$cds_end > exTab$exon_end,  exTab$exon_end,
           ifelse(exTab$cds_end < exTab$exon_start,
                  rep(NA, nrow(exTab)), exTab$cds_end))

  start <- ifelse(start > end, rep(NA, nrow(exTab)), start)
  end <- ifelse(end < start, rep(NA, nrow(exTab)), end)

  exTab$cds_end <- end
  exTab$cds_start <- start
  exTab
}


.createTxDbFromEnsembl <-
  function(version, annotInfo, ensemblApi, bioperl,
           species, user, host, pass, chr, ...)
{
  fn <- system.file("script", "get-tables-required-for-TxDb.pl",
                    package="EnhancedTxDbs")

  ## create the tables needed for creating a TxDb object
  Sys.setenv(ENS=ensemblApi)

  ## add ensembl API to the perl search path
  con <- file("./bash.txt", open="w")
  writeLines(text=paste0("PERL5LIB=${PERL5LIB}:", normalizePath(bioperl), "\n PERL5LIB=${PERL5LIB}:",
             normalizePath(ensemblApi), "API/ensembl/modules\nPERL5LIB=${PERL5LIB}:",
             normalizePath(ensemblApi),
             "API/ensembl-compara/modules\nPERL5LIB=${PERL5LIB}:", normalizePath(ensemblApi),
             "API/ensembl-variation/modules\nPERL5LIB=${PERL5LIB}:", normalizePath(ensemblApi),
             "API/ensembl-functgenomics/modules\nexport PERL5LIB\n"), con=con)
  close(con)

  system("source ./bash.txt")
  unlink("./bash.txt")

  cmd <- paste("perl ", fn, " -s ",
               species," -e ", version, " -w ", chr,
               " -d core -u ", user, " -h ",
               host, " -p ", pass, sep="")
  system(cmd)

  Sys.unsetenv("ENS")

  ## Read chromosome information from perl out put
  chrTab <- read.table("chr.txt", sep="\t", header=TRUE,
                       stringsAsFactors=FALSE)
  system("rm chr.txt")
  rownames(chrTab) <- chrTab$chrom
  chrTab$is_circular <- ifelse(chrTab$is_circular == 0, FALSE, TRUE)

  ## Read the exon table created by the perl script
  exTab <- read.table("exons.txt", sep="\t", header=TRUE,
                      stringsAsFactors=FALSE,
                      na.strings = c("NA", "NULL"))

  exTab <- exTab[exTab$chrom %in% chrTab$chrom, ]

  system("rm exons.txt")

  ## Ordering of the exon table according chromosome, gene, transcript
  ## and genomic coordinates
  exTab <- exTab[order(exTab$chrom, exTab$gene_id, exTab$tx_name,
                       exTab$exon_start*exTab$strand,
                       exTab$exon_end*exTab$strand), ]


  ## create transcript ids
  exTab$tx_id <- as.integer( factor(paste(exTab$tx_name, exTab$chrom, sep=":")) )

  ## create the exon id
  exTab$exon_id <- as.integer(factor(exTab$exon_name))

  exTab <- computeExonRanks(exTab)
  exTab <- computeExonCds(exTab)

  exTab$strand <- ifelse(exTab$strand == -1, "-", "+")


  ## Creating the transcript table
  txInfo <- c("tx_id", "tx_name", "chrom",
              "strand", "tx_start",
              "tx_end", "transcript_biotype")

  transcripts <- exTab[!duplicated(exTab$tx_id), txInfo]


  idx <- grep(colnames(transcripts), pattern="^chrom$|^strand$")
  colnames(transcripts)[idx] <- paste("tx_", colnames(transcripts)[idx], sep="")


  ## Creating the splicings table
  splInfo <- c("tx_id", "exon_rank", "exon_id", "exon_name", "chrom",
               "strand", "exon_start", "exon_end", "cds_start", "cds_end")


  splicings <- unique(exTab[, splInfo])


  idx <- grep(colnames(splicings), pattern="^chrom$|^strand$")
  colnames(splicings)[idx] <- paste("exon_", colnames(splicings)[idx], sep="")


  genes <- unique(exTab[, c("tx_id", "gene_id", "gene_name",
                            "gene_biotype", "chrom")])

  genes <- genes[, colnames(genes) != "chrom"]

  ## chromosomes must be supplied as character to makeTranscriptDb
  transcripts$tx_chrom <- as.character(transcripts$tx_chrom)
  splicings$exon_chrom <- as.character(splicings$exon_chrom)
  chrTab$chrom <- as.character(chrTab$chrom)


  ## exon id must be unique per chromosome otherwise problems with
  ## makeTranscriptDb function.
  splicings$exon_id = as.integer(factor(paste(splicings$exon_id,
    splicings$exon_chrom, sep=":")))

  ## Generate the TranscriptDb object
  txdb <- makeTranscriptDb(transcripts[, !colnames(transcripts) %in%
                                       c("transcript_biotype")],
                           splicings,
                           genes=genes[, ! colnames(genes) %in%
                             c("gene_name", "gene_biotype")],
                           chrominfo=chrTab)

  if(annotInfo) {

    ## Remove transcript ids from the table
    genes <- unique(genes[, -1, drop=FALSE])


    ## get a transcripts by genes object
    txsByGns <- transcriptsBy(txdb, "gene")
    exsByTxs <- exonsBy(txdb, "tx")
    exsByGns <- exonsBy(txdb, "gene")


    ## add biotype information
    gn <- rep(names(txsByGns), elementLengths(txsByGns))
    cat(paste(class(txsByGns), "\n"))

    txsByGns.flat <- GenomicRanges::unlist(txsByGns, use.names=FALSE)


    txn <- paste(as.character(seqnames(txsByGns.flat)),
                 values(txsByGns.flat)[["tx_name"]], sep=":")
    values(txsByGns.flat)["transcript_biotype"] <-
      transcripts[txn, "transcript_biotype"]
    txsByGns <- split(txsByGns.flat, gn)


    rownames(transcripts) <- transcripts$tx_id
    values(exsByTxs) <- transcripts[names(exsByTxs),]


    txdb <- EnhancedTxDb(txdb=txdb,
                         txAnnot=transcripts[, c("transcript_biotype"), drop=FALSE],
                         gnAnnot=genes, exAnnot=exTab)

  }

  txdb

}

createTxDbFromEnsembl <-
  function(version, annotInfo=TRUE,
           ensemblApi="/home/bioinfo/ensembl/73/", bioperl="~/ENSEMBL/bioperl-1.2.3",
           species = "human", user = "anonuser", host = "madb.i-med.ac.at", pass = "", chr=character(), ...)
{
  chr <- ifelse(length(chr) == 0, "NA", paste("'", paste(chr, collapse = " "), "'", sep=""))
  .createTxDbFromEnsembl(version=version, annotInfo=annotInfo, ensemblApi=ensemblApi, bioperl=bioperl,
                         species=species, user=user, host=host, pass=pass, chr=chr, ...)
}







