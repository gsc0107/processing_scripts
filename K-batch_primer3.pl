#!usr/bin/perl

use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Getopt::Long;

my $usage=
"perl script.pl <options> 
-i|--input  input sequences in fasta format
-o|--output output file
";

my $filein;
my $fileout="primer3";

GetOptions(

'i|input:s'  => \$filein,
'o|output:s' => \$fileout,

);

die $usage unless defined ($filein);


my $inseq    = Bio::SeqIO->new(
							-file => $filein,
						    -format => 'Fasta',
							);

open (FILEOUT1, ">", "$fileout.primers.txt");

while (my $seq = $inseq->next_seq){

	my $primer3_output= run_primer3($seq->id,$seq->seq());
	print FILEOUT1 $primer3_output, "\n";
	
	my $primer_left_0 = qx(echo "$primer3_output" | grep PRIMER_LEFT_0_SEQUENCE | cut -f 2 -d "=");
	chomp $primer_left_0;
	my $primer_left_0_tm = qx(echo "$primer3_output" | grep PRIMER_LEFT_0_TM | cut -f 2 -d "=");
	chomp $primer_left_0_tm;
	print $seq->id, "_fw", "\t", $primer_left_0, "\t", $primer_left_0_tm, "\n";
	
	my $primer_right_0 = qx(echo "$primer3_output" | grep PRIMER_RIGHT_0_SEQUENCE | cut -f 2 -d "=");
	chomp $primer_right_0;
	my $primer_right_0_tm = qx(echo "$primer3_output" | grep PRIMER_RIGHT_0_TM | cut -f 2 -d "=");
	chomp $primer_right_0_tm;
	print $seq->id, "_rev", "\t", $primer_right_0, "\t", $primer_right_0_tm, "\n";


}

sub run_primer3{

my ($id, $seq)=@_;

my $length=length($seq) . "-" . length($seq);

my $input="SEQUENCE_ID=$id
SEQUENCE_TEMPLATE=$seq
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=12
PRIMER_MAX_SIZE=35
PRIMER_MAX_NS_ACCEPTED=0
PRIMER_MIN_GC=10
PRIMER_PRODUCT_SIZE_RANGE=$length
P3_FILE_FLAG=1
PRIMER_EXPLAIN_FLAG=1
PRIMER_MIN_TM=55
PRIMER_MAX_TM=68
=";

my $output = qx(echo "$input" | primer3_core);

return $output;

}