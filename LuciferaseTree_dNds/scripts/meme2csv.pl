#!/usr/bin/perl
use strict;

my $filename = $ARGV[0];
open(my $fh, '<:encoding(UTF-8)', $filename)
  or die "Could not open file '$filename' $!";
 
while (my $row = <$fh>) {
  chomp $row;
  #Check for multiple pipes in a line, which is used in text table of meme results
  my $number = () = $row =~ /\|/gi;
  $row =~ s/ //g;
  if($number == 9){
	my @columns = split /\|/, $row;
	shift(@columns); #row leads with a pipe, so remove that
	my $csv = join ",", @columns;
	my $dashes = () = $csv =~ /\-/, $csv;
	if($dashes < 1){	#table has row of dashes for cosmetics, not necessary here
		print $csv."\n";
	}
  }
}
