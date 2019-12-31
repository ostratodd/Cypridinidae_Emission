#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

#This script converts a fasta sequence file into a comma delimited text file, e.g. to read into R
#USAGE: fasta2csv.pl <fastafile>
#It requires BioPerl's SeqIO package. Simply writes to screen. Use redirect > to write to a file.
# 
<<<<<<< HEAD
=======
#special use is that if sequence name contains __Lmax### use the ### for the last column for mutant analysis of color
>>>>>>> 94603285f2da46352c5f9e871aa8e8b05c35abc0

my $datafile = $ARGV[0];

my $seq_in = Bio::SeqIO->new(-file => $datafile); 
<<<<<<< HEAD

while ( my $seq = $seq_in->next_seq() ) {
=======
my $lmax;
my @lmaxes;

while ( my $seq = $seq_in->next_seq() ) {
	if($seq->id =~ m/Lmax/) {
		@lmaxes = split /Lmax/, $seq->id;
		$lmax = $lmaxes[1];
	}
>>>>>>> 94603285f2da46352c5f9e871aa8e8b05c35abc0
	print $seq->id.",";
	my $theseq = $seq->seq;
	my @residues = split //, $theseq;
	my $csv = join ",", @residues;
<<<<<<< HEAD
	print "$csv\n";
=======
	print "$csv";
	if($lmax){
		print ",".$lmax;
	}
	print "\n";
>>>>>>> 94603285f2da46352c5f9e871aa8e8b05c35abc0
}
