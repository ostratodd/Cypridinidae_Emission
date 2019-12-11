#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

#This script reads an aligned fasta file. After specifying a reference sequence from the file, 
#it creates a table of aligned sites' correspondence to the reference sequence
#USAGE: fasta2csv.pl <fastafile> "sequence name as reference"
#It requires BioPerl's SeqIO package. Simply writes to screen. Use redirect > to write to a file.
# 

my $datafile = $ARGV[0];
my $refseq = $ARGV[1];
my $seq_in;
my $alignN = 0;
my $refN = 0;
my @residues;

$seq_in = Bio::SeqIO->new(-file => $datafile); 
while ( my $seq = $seq_in->next_seq() ) {
	if($seq->id eq $refseq) {
		print "Aligned,".$seq->id."\n";
		my $theseq = $seq->seq;
		@residues = split //, $theseq;
		for(my $i=0; $i < @residues; $i++){
			if($residues[$i] eq '-') {
				$alignN++;				
			}else{
				$alignN++;				
				$refN++;	
			}
			print $alignN.",";
			print $refN."\n";
		}
	}
}
