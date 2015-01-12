getChrominfo <-
    function(version, ensemblApi, species, user="anonymous", host="ensembldb.ensembl.org", pass="", outfile, ...){
        fn <- system.file("script", "get-chrominfo.pl",
                          package="EnhancedTxDbs")
        if( missing( species ) ){
            species <- "human"
        }
        if( missing( outfile ) ){
            outfile <- paste0( "chrominfo_Ensembl_", version, "_", species, ".txt" )
        }
        ## create the tables needed for creating a TxDb object
        Sys.setenv(ENS=ensemblApi)

        ## add ensembl API to the perl search path
        ##con <- file("./bash.txt", open="w")
        ##writeLines(text=paste0("PERL5LIB=${PERL5LIB}:", normalizePath(bioperl), "\n PERL5LIB=${PERL5LIB}:",
        ##writeLines(text=paste0( "PERL5LIB=${PERL5LIB}:",
        ##             normalizePath(ensemblApi), "API/ensembl/modules\nPERL5LIB=${PERL5LIB}:",
        ##             normalizePath(ensemblApi),
        ##             "API/ensembl-compara/modules\nPERL5LIB=${PERL5LIB}:", normalizePath(ensemblApi),
        ##             "API/ensembl-variation/modules\nPERL5LIB=${PERL5LIB}:", normalizePath(ensemblApi),
        ##             "API/ensembl-functgenomics/modules\nexport PERL5LIB\n"), con=con)
        ##close(con)

        ##system("source ./bash.txt")
        ##unlink("./bash.txt")

        cmd <- paste0("perl ", fn, " -s ", species," -e ", version,
                     " -U ", user, " -H ", host, " -P ", pass, " -o ", outfile )
        system(cmd)
        Sys.unsetenv("ENS")

        Dat <- read.table( outfile, sep="\t", as.is=TRUE, header=TRUE )
        return( Dat )
}
