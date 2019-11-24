#!/bin/perl

use strict;
use Bio::SeqIO;

# get the file name, somehow 
my $file = shift; 
my $seqio_object = Bio::SeqIO->new(-file => $file); 
my $seq_object = $seqio_object->next_seq;
