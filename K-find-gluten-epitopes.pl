#!/bin/perl
 
use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::SeqIO;

my ($file_pep, $file_nt, $out_file, $motif_file);

GetOptions(
       'p|inseq_pep:s'        => \$file_pep,
       'n|inseq_nt:s'        => \$file_nt,
       'o|outfile:s'        => \$out_file,
       'm|motif:s'        => \$motif_file,
#       'a|annotation:s'   => \$annotation_file,
    );

my $usage="                                                                                                                                                                                                                                                                     
Usage: perl script.pl <options>                                                                                                                                                                                                                                                 
Required options:                                                                                                                                                                                                                                                               
       '-p|--inseq_pep'  input fasta protein file
       '-n|--inseq_nt'   input fasta nucleotide file
       '-o|--outfile'    prefix for all output files
       '-m|--motif'      input tabular file with motives
";

die $usage unless ( (defined $file_pep) and (defined $motif_file) );

# Load all fasta sequences

my %sequences;

my $inseq_pep    = Bio::SeqIO->new(
									-file => $file_pep,
									-format => 'Fasta',
									);

while (my $seq = $inseq_pep->next_seq){

	my ($header)=split(" ", $seq->id);
	$sequences{$header}{'pep'}=$seq->seq();

}

if (defined $file_nt){

my $inseq_nt    = Bio::SeqIO->new(
									-file => $file_nt,
									-format => 'Fasta',
									);

while (my $seq = $inseq_nt->next_seq){

	my ($header)=split(" ", $seq->id);
	$sequences{$header}{'nt'}=$seq->seq();

}

}
# Load gluten family annotation


# Load immunogenic epitopes

# Epitope_Name	Epitope_Old_Name	Epitope_sequence_wheat	Modified_sequence_human	Critical_aa	Disease	Reference
# Epitope_sequence_wheat contains redundant entries

open FILE1, "<", $motif_file or die $usage;

my %epitopes_db;
my $header=<FILE1>;

while (my $line=<FILE1>){

	chomp $line;
	my @info=split("\t", $line);
	
	my $wheat_seq=$info[2];
	
	$epitopes_db{$wheat_seq}{'Name'}=$info[0]; 
	$epitopes_db{$wheat_seq}{'Wheat_sequence'}= $info[2];

#	$epitopes_db{$wheat_seq}{'Old_Name'} = $info[1];
#	push $epitopes_db{$wheat_seq}{'Critical_aa'}
	$epitopes_db{$wheat_seq}{'Disease'}= $info[5];
#	push $epitopes_db{$wheat_seq}{'Reference'}

}
#search for epitopes in each sequence

my %epitopes_present;

foreach my $seq_id (keys %sequences){

	foreach my $epitope (keys %epitopes_db){
		
		if ($epitopes_db{$epitope}{'Disease'} eq 'Celiac'){
		
			if ($sequences{$seq_id}{'pep'}=~ m/$epitope/g ){
			
				$epitopes_present{$seq_id}{$epitope}{'Positions'} = match_all_positions($epitope, $sequences{$seq_id}{'pep'});
				$epitopes_present{$seq_id}{$epitope}{'Count'} = count_matches($epitope, $sequences{$seq_id}{'pep'}); 
			}

		}
	}		
}

close FILE1;

#print header

print "Contig\tEpitopes\tTotal_epitopes\tProtein\tNucleotide_ORF\n";

foreach my $seq_id (keys %sequences){

 if ( defined $epitopes_present{$seq_id}){
    
    my @output;
    my $total_epitopes=0;
     
    foreach my $k (keys %{$epitopes_present{$seq_id}}){
    	
	    	my $output_epitope= $k . "[" . $epitopes_present{$seq_id}{$k}{'Count'} ."]";
	   	 	push @output, $output_epitope;
	    	$total_epitopes+=$epitopes_present{$seq_id}{$k}{'Count'};
    
    }
 	
 	print $seq_id, "\t", join (",", @output), "\t", $total_epitopes, "\t", 
 	$sequences{$seq_id}{'pep'}, "\t", $sequences{$seq_id}{'nt'}, "\n";
 
 }
 
 else{
 
 	print $seq_id, "\t", "No immunogenic epitopes found", "\t", "0", "\t", 
 	$sequences{$seq_id}{'pep'}, "\t", $sequences{$seq_id}{'nt'}, "\n";
 
 }


}

sub match_all_positions {
    my ($regex, $string) = @_;
    my @ret;
    while ($string =~ /$regex/g) {
        push @ret, [ $-[0], $+[0] ];
    }
    return @ret;
}

sub count_matches{

 my ($regex, $string) = @_;
 my @count = $string =~ /$regex/g;
 
 return scalar @count; 

}
