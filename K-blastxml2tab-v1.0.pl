#!/usr/bin/perl -w

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

$XML::SAX::ParserPackage = 'XML::SAX::PurePerl';
 
my ($filein, $fileout);

GetOptions(
       'i|input:s'        => \$filein,
       'o|output:s'        => \$fileout,
#       'e|evalue:s'        => \$evalue,
);  

open (FILEOUT1, ">", $fileout);

my $in = Bio::SearchIO->new (-format => 'blastxml', 
                             -file   => $filein);
                           
while( my $result = $in->next_result ) {

    print $result, "\n";

  while( my $hit = $result->next_hit ) {
    
    while( my $hsp = $hit->next_hsp ) {

#	print $result->query_description, "\n";
    
    print FILEOUT1 $result->query_description, "\t", #query name

    	  		  $hit->name, "\t", #subject name
    	  		  $hsp->percent_identity, "\t", #percent identity
    	  		  $hsp->hsp_length, "\t", #alignment length 
    	          ($hsp->hsp_length - $hsp->gaps - $hsp->num_identical), "\t", #mismatches
    	  		  $hsp->gaps, "\t", #gaps
    	  		  $hsp->start('query'), "\t", #query start
    	  		  $hsp->end('query'), "\t", #query end
    	  		  $hsp->start('hit'), "\t", #subject start    	
    	  		  $hsp-> end('hit'), "\t", #subject end
    	  		  $hsp->evalue, "\t", #e-value
    	  		  $hsp->bits, "\n", #bit score

	}
	
 }
 	
}

close FILEOUT1;
