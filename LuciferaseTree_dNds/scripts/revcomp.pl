#!/usr/bin/perl
use strict;

my $filename = $ARGV[0];
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
while (my $row = <$fh>) {
  chomp $row;
  my $dna=$row;
  my $rcdna= & reverse_complement_IUPAC($dna);
  print "$rcdna\n";
}

sub reverse_complement_IUPAC {
        my $dna = shift;

	if ($dna =~ /\>/) {
		return($dna);
    	}


        # reverse the DNA sequence
        my $revcomp = reverse($dna);

        # complement the reversed DNA sequence
        $revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
        return $revcomp;
}
