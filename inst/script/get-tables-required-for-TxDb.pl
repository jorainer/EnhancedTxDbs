## This script retrives the exon and transcript information required
## for creating a TranscriptDb object later in R

my $script_version="0.0.2";

## uses environment variable ENS pointing to the
## ENSEMBL API on the computer
use lib $ENV {ENS};
use strict;
use Getopt::Std;

## unification function for arrays
use List::MoreUtils qw/ uniq /;

## initialize option variables
my %option=();
my $ensembl_version="none";
my $ensembl_database="core";
my $species = "human";
my $chrPreSel="NA";


## connecting to the ENSEMBL data base
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::ApiVersion;
#my $user = "anonuser";
#my $host = "madb.i-med.ac.at";
my $user = "anonymous";
my $host = "ensembldb.ensembl.org";
my $pass = "";
my $registry = 'Bio::EnsEMBL::Registry';

getopts("e:d:H:hP:U:s:w:",\%option);

if( defined( $option{ h } ) ){
  ## print help and exit.
  print( "\nget-tables-required-for-TxDb version ".$script_version.".\n" );
  print( "Retrieves gene/transcript/exon annotations from Ensembl and stores it as tabulator delimited text files that can be used to generate a Bioconductor TxDb database.\n\n" );
  print( "usage: perl get-tables-required-for-TxDb -e:d:P:hH:U:s:w\n" );
  print( "-e (required): the Ensembl version (e.g. -e 75). The function will internally check if the submitted version matches the Ensembl API version and database queried.\n" );
  print( "-d (optional): the Ensembl database, defaults to core.\n" );
  print( "-H (optional): the hostname of the Ensembl database; defaults to ensembldb.ensembl.org.\n" );
  print( "-h (optional): print this help message.\n" );
  print( "-P (optional): the password to access the Ensembl database.\n" );
  print( "-U (optional): the username to access the Ensembl database.\n" );
  print( "-s (optional): the species; defaults to human.\n" );
  print( "-w (optional): a list of chromosomes; by default all genes will be fetched.\n" );
  exit 0;
}

if( defined($option{ d } ) ){
	$ensembl_database=$option{ d };
}
if( defined($option{ s } ) ){
	$species=$option{ s };
}

if( defined( $option{ U } ) ){
	$user=$option{ U };
}
if( defined( $option{ H } ) ){
	$host=$option{ H };
}
if( defined( $option{ P } ) ){
	$pass=$option{ P };
}
if( defined( $option{ w } ) ){
	$chrPreSel=$option{ w };
}
if( defined( $option{ e } ) ){
	$ensembl_version=$option{ e };
}
else{
	die("the ensembl version has to be specified with the -e parameter (e.g. -e 75)");
}

my $api_version="".software_version()."";
if( $ensembl_version ne $api_version ){
    die "The submitted Ensembl version (".$ensembl_version.") does not match the version of the Ensembl API (".$api_version."). Please configure the environment variable ENS to point to the correct API.";
}


$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -verbose => "1" );

my $gene_adaptor = $registry->get_adaptor( $species, $ensembl_database, "gene" );
my $slice_adaptor = $registry->get_adaptor( $species, $ensembl_database, "slice" );


## provide some information about the script
print "\nCreating tables required for TxDb objects!\n \
Querying Ensembl database $ensembl_database, species: $species";

if( defined( $option{a} ) ){
print "\nadding additional transcript definitions from the table $option{a}";
}
print "\n\n";

my $slice;
## get all gene ids defined in the database...
my @gene_ids = ();
if ($chrPreSel eq "NA") {
  @gene_ids = @{$gene_adaptor->list_stable_ids()};
} else {
  my @retrChrs = split(' ', $chrPreSel);
  foreach my $retrChr (@retrChrs ){
    $slice = $slice_adaptor->fetch_by_region( 'chromosome', $retrChr);
    @gene_ids = (@gene_ids, @{ $gene_adaptor->fetch_all_by_Slice($slice) });
  }

}

my $counta=0;



## chromosome list
my @chr_ids = ();

## write exon information into table
my $outfile="exons.txt";
open( OUT , ">$outfile");

print OUT "gene_id\tgene_name\tgene_biotype\tchrom\ttx_name\ttx_start";
print OUT "\ttx_end\tstrand\tcds_start\tcds_end\texon_name\texon_start";
print OUT "\texon_end\ttranscript_biotype\tcoord_system\n";


foreach my $gene_id ( @gene_ids ){
	$counta++;
	if( ($counta % 2000) == 0 ){
		print "processed $counta genes\n";
	}
	my $orig_gene;

	## check if gene id or gene
	if ($chrPreSel eq "NA") {
	  $orig_gene = $gene_adaptor->fetch_by_stable_id( $gene_id );
	} else {
	  $orig_gene = $gene_id;
	  $gene_id = $gene_id->stable_id;
	}
	if( defined $orig_gene ){
	  my $do_transform=1;
	  my $gene  = $orig_gene->transform("chromosome");
	  if( !defined $gene ){
	    #next;
	    ## gene is not on known defined chromosomes!
	    $gene = $orig_gene;
	    $do_transform=0;
	  }
	  my $coord_system = $gene->coord_system_name;
	  my $chrom = $gene->slice->seq_region_name;

	  ## store chromosomes per genes

	  push(@chr_ids, $chrom);
	  my $gene_external_name= $gene->external_name;
	  my $gene_biotype = $gene->biotype;
	  my @transcripts = @{ $gene->get_all_Transcripts };
	  ## ok looping through the transcripts
	  foreach my $transcript ( @transcripts ){
	    if( $do_transform==1 ){
	      ## just to be shure that we have the transcript in chromosomal coordinations.
	      $transcript = $transcript->transform("chromosome");
	    }
	    my $tx_start = $transcript->start;
	    my $tx_end = $transcript->end;
	    my $strand = $transcript->strand;

	    ## caution!!! will get undef if transcript is non-coding!
	    my $transcript_cds_start = $transcript->coding_region_start;
	    if( !defined( $transcript_cds_start ) ){
	      $transcript_cds_start = "NULL";
	    }
	    my $transcript_cds_end = $transcript->coding_region_end;
	    if( !defined( $transcript_cds_end ) ){
	      $transcript_cds_end = "NULL";
	    }
	    my $tx_name = $transcript->stable_id;
	    my $transcript_biotype = $transcript->biotype;
	    my @exons = @{ $transcript->get_all_Exons() };
	    foreach my $exon (@exons){
	      if( $do_transform==1 ){
		$exon->transform("chromosome");
	      }
	      my $exon_start = $exon->start;
	      my $exon_end = $exon->end;
	      my $exon_id = $exon->stable_id;
	      ## insert all this into the database...
	      my $query_string = "$gene_id\t$gene_external_name\t$gene_biotype";
	      $query_string=$query_string."\t$chrom\t$tx_name\t$tx_start\t$tx_end";
	      $query_string=$query_string."\t$strand\t$transcript_cds_start";
	      $query_string=$query_string."\t$transcript_cds_end\t$exon_id\t";
	      $query_string=$query_string."$exon_start\t$exon_end\t";
	      $query_string=$query_string."$transcript_biotype\t$coord_system";
	      print OUT "$query_string\n";
	    }
	  }
	}
}

close(OUT);

## unify chromosme list
my @chr_ids = uniq @chr_ids;

## retrive chromosome information for all chromosomes
## in contained in the exon table


open( OUT , ">chr.txt");
print OUT "chrom\tlength\tis_circular\n";

foreach my $chr_id ( @chr_ids ){
  ## get a slice on the entire chromosome
  my $chr_slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr_id );
  my $name = $chr_slice->seq_region_name;
  my $length = $chr_slice->length;
  my $is_circular = $chr_slice->is_circular;
  print OUT "$name\t$length\t$is_circular\n";
}

close(OUT)
