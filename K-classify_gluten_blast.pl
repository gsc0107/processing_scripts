#!/usr/bin/perl                                                                                                                                                                                                                                                                 
# Author: Ksenia Krasileva
# Contact: ksenia.krasileva@gmail.com                                                                                                                                                                                                  
                                           
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $usage="                                                                                                                                                                                                                                                                     
Usage: perl script.pl <options>                                                                                                                                                                                                                                                 
Required options:                                                                                                                                                                                                                                                               
       1|blast-nr    tabular blastx against nr ('6 std qcovs stitle') 
       2|blast-arms  tabular blastn against genome ('6 std qcovs')     
       f|fasta       input sequences in fasta format
       o|outfile     prefix for all output files 
";

my ($blast_nr_file, $blast_arms_file, $fasta_file, $out_file, $error_file, $out_blast, $out_fasta);

#Defaults                                                                                                                                                                                                                                                                       
#my $pident_cutoff = '0';
#my $qcov_cutoff = '0';

GetOptions(
       '1|blast-nr:s'        => \$blast_nr_file,
       '2|blast-arms:s'        => \$blast_arms_file,
       'f|fasta:s'        => \$fasta_file,
       'o|outfile:s'        => \$out_file,
       'e|errorfile:s'        => \$error_file,
       'b|outblast:s'        => \$out_blast,
       'of|outfasta:s'        => \$out_fasta,
#       'i|pid-cutoff:s'        => \$pident_cutoff,
#       'q|qcov-cutoff:s'        => \$qcov_cutoff,
    );

die $usage unless defined($blast_nr_file);

# Strings describing gluten genes in nr
# "alpha-gliadin" or "alpha/beta-gliadin" or "alpha gliadin" or "alfa gliadin"
# "gamma-gliadin"
# "omega gliadin"
# "LMW-glutenin" or "low molecular weight glutenin" or "lmw-gs" or "low-molecular-weight"
# "high molecular weight glutenin" or "HMW glutenin" or "high-molecular-weight"

# blast format (6 custom)
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle


my %blast_nr = load_blast($blast_nr_file);
#my %blast_arms = load_blast($blast_arms_file);

my @blast_format= ('query', 'subject', 'pident', 'length', 'mismatch', 'gapopen', 'qstart','qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'stitle');

# blast data structure
# $blast{$query}{$n}{'query'} = $query;

# classify sequences by gluten class by evaluating the top hit

my %gluten_class;

my ($qseqid, $stitle);

foreach my $key (sort keys %blast_nr){

    $qseqid = $blast_nr{$key}{'1'}{'query'};
    $stitle = $blast_nr{$key}{'1'}{'stitle'};
    
    if ($stitle=~ /glutenin|(lmw-gs)/i){

	if ($stitle=~ /(LMW)|(low)/i){
	
		$gluten_class{'lmw'}{$qseqid} = 1;
	
	}

	elsif ($stitle=~ /(HMW)|(high)/i){
	
		$gluten_class{'hmw'}{$qseqid} = 1;
	
	}

	else {
	
	    $gluten_class{'unknown'}{$qseqid} = 1; 
    
	}
    }

    elsif ($stitle=~ /gliadin/i){

	if ($stitle=~ /alpha/i){
	
		$gluten_class{'alpha'}{$qseqid} = 1;
	
	}

	elsif ($stitle=~ /gamma/i){
	
		$gluten_class{'gamma'}{$qseqid} = 1;
	
	}

	elsif ($stitle=~ /omega/i){
	
		$gluten_class{'omega'}{$qseqid} = 1;
	
	}

	else{

	    $gluten_class{'unknown'}{$qseqid} = 1;

	}
    }

    else{

	$gluten_class{'other'}{$qseqid} = 1;

    }

}
# create output files

#open OUTFILE_tab "gluten-all.tab";

    foreach my $g_class (sort keys %gluten_class){
	
	foreach my $id ( sort keys %{ $gluten_class{$g_class} } ){
	    
	    print $id, "\t", $g_class, "\t", $blast_nr{$id}{'1'}{'stitle'}, "\n";

			}
    }

#open OUTFILE_LMW;
#open OUTFILE_HMW;
#open OUTFILE_alpha;
#open OUTFILE_gamma;
#open OUTFILE_omega;
#open OUTFILE_other;


sub load_blast {

    my $filein = shift(@_);

    open (my $FILE1, "<", $filein) or die "cannot open $filein";

    my %blast;
    my $n=1;

    while ( my $line = <$FILE1> ) {

	chomp $line;
	my ($query, $subject, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $qcovs, $stitle) = split ("\t", $line);

	if ($blast{$query}){
	    $n++;
	}
	else{
	    $n=1;
	}

	$blast{$query}{$n}{'query'} = $query;
	$blast{$query}{$n}{'subject'} = $subject;
	$blast{$query}{$n}{'pident'} = $pident;
	$blast{$query}{$n}{'length'} = $length;
	$blast{$query}{$n}{'mismatch'} = $mismatch;
	$blast{$query}{$n}{'gapopen'} = $gapopen;
	$blast{$query}{$n}{'qstart'} = $qstart;
	$blast{$query}{$n}{'qend'} = $qend;
	$blast{$query}{$n}{'sstart'} = $sstart;
	$blast{$query}{$n}{'send'} = $send;
	$blast{$query}{$n}{'evalue'} = $evalue;
	$blast{$query}{$n}{'bitscore'} = $bitscore;
	$blast{$query}{$n}{'qcovs'} = $qcovs;
	$blast{$query}{$n}{'stitle'} = $stitle;

    }

close $FILE1;

return %blast;

}
