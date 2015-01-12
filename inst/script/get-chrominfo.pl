my $script_version="0.0.1";

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
my $user = "anonymous";
my $host = "ensembldb.ensembl.org";
my $pass = "";
my $registry = 'Bio::EnsEMBL::Registry';
my $api_version="".software_version()."";

getopts("e:d:H:hP:U:s:o:",\%option);

if( defined( $option{ h } ) ){
  ## print help and exit.
  print( "\nget-chrominfo.pl version ".$script_version.".\n" );
  print( "Retrieves chromosome information required to build a TranscriptDB object by the GenomicFeatures package.\n\n" );
  print( "usage: perl get-chrominfo.pl -e:d:H:h:P:U:s:\n" );
  print( "-e (required): the Ensembl version (e.g. -e 75). The function will internally check if the submitted version matches the Ensembl API version and database queried.\n" );
  print( "-d (optional): the Ensembl database, defaults to core.\n" );
  print( "-H (optional): the hostname of the Ensembl database; defaults to ensembldb.ensembl.org.\n" );
  print( "-h (optional): print this help message.\n" );
  print( "-P (optional): the password to access the Ensembl database.\n" );
  print( "-U (optional): the username to access the Ensembl database.\n" );
  print( "-s (optional): the species; defaults to human.\n" );
  print( "-o (optional): the output file name; defaults to chrominfo-Ensembl-$api_version-$species.txt" );
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

if( $ensembl_version ne $api_version ){
    die "The submitted Ensembl version (".$ensembl_version.") does not match the version of the Ensembl API (".$api_version."). Please configure the environment variable ENS to point to the correct API.";
}

my $outfile = "chrominfo-Ensembl-$api_version-$species.txt";
if( defined( $option{o} ) ){
  $outfile=$option{o};
}

$registry->load_registry_from_db(-host => $host, -user => $user,
				 -pass => $pass, -verbose => "1" );

my $gene_adaptor = $registry->get_adaptor( $species, $ensembl_database, "gene" );
my $slice_adaptor = $registry->get_adaptor( $species, $ensembl_database, "slice" );

my @gene_ids = @{$gene_adaptor->list_stable_ids()};
## chromosome list
my @chr_ids = ();
my $counta=0;
foreach my $gene_id ( @gene_ids ){
  $counta++;
  if( ($counta % 2000) == 0 ){
    print "processed $counta genes\n";
  }
  my $tmp_gene = $gene_adaptor->fetch_by_stable_id( $gene_id );
  my $gene  = $tmp_gene->transform("chromosome");
  if( !defined $gene ){
    #next;
    ## gene is not on known defined chromosomes!
    $gene = $tmp_gene;
  }
  my $coord_system = $gene->coord_system_name;
  my $chrom = $gene->slice->seq_region_name;
  push(@chr_ids, $chrom);
}
## unify chromosme list
my @chr_ids = uniq @chr_ids;

open( OUT , ">".$outfile );
print OUT "chrom\tlength\tis_circular\n";

foreach my $chr_id ( @chr_ids ){
  ## get a slice on the entire chromosome
  my $chr_slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr_id );
  if( defined $chr_slice ){
    my $name = $chr_slice->seq_region_name;
    my $length = $chr_slice->length;
    my $is_circular = $chr_slice->is_circular;
    print OUT "$name\t$length\t$is_circular\n";
  }
}


## that would be such a nice and fast option, if the gtf would not contain all kind of genes
## encoded on all type of patched chromosomes.
#foreach my $chr ( @{ $slice_adaptor->fetch_all( "chromosome", undef, undef, 1 ) } ){
#  print OUT $chr->seq_region_name."\t".$chr->length."\t".$chr->is_circular."\n";
#}
close( OUT );


