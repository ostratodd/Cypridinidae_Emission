#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;


#This script reads an aligned fasta file with a sequence, and a text file with mutations. It ouputs
#the full sequence for each mutant, after adding mutations

#USAGE: fasta2csv.pl <fastafile> <Mutations File>
#It requires BioPerl's SeqIO package. Simply writes to screen. Use redirect > to write to a file.
# 

my $datafile = $ARGV[0];
my @residues;

my $seq_in = Bio::SeqIO->new(-file => $datafile); 
my $seq = $seq_in->next_seq();
my $wild = $seq->seq;
my $wtlmax = $ARGV[2];

print ">Wild_type__Lmax".$wtlmax."\n";
print $wild."\n";

#read mutation file
my @mutations;
my @mutant;
my $mutfile = $ARGV[1];
open(my $mfh, '<:encoding(UTF-8)', $mutfile)
  or die "Could not open file '$mutfile' $!";

my $count =1;
my $phenotype;
while (my $row = <$mfh>) {
	chomp $row;
	@mutations = split("\t", $row);
	my $newsequence = $wild;
	for(my $i =0; $i < @mutations; $i++){
		if($i == @mutations-1 ){
			$phenotype = $mutations[$i];
		}else{
			@mutant = split(/\,/ , $mutations[$i]);
			my $old = $mutant[0];
			my $site = $mutant[1];
			my $new = $mutant[2];
			my $tomutate = substr($wild, $site - 1,1);
			if($tomutate eq $old){
				$newsequence = substr($newsequence, 0, $site - 1).$new.substr($newsequence, $site);
			}else{
				print "\n\n****ERROR: $tomutate in wild type $old in mutation file do not match for site $site\n";
				exit;
			}
		}
	}
	print ">Mutation_$count"."__Lmax$phenotype\n";
	print $newsequence."\n";
	$count++;
}
