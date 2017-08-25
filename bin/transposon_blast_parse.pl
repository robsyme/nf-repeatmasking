#!/usr/bin/perl -w
 
# transposon_blast_parse.pl
# This script will parse the results of the unknown elements of RepeatModeler blasted against the transposon database into two files, one with identified elements from blast hits and the other of still unknown elements. 

# Megan Bowman
# 16 March 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my $usage = "\nUsage: $0 --blastx <blastx output file> --modelerunknown <full path of the modeler unknown file>\n";

my ($blastx, $m_unknown);

GetOptions('blastx=s' => \$blastx,
	  'modelerunknown=s' => \$m_unknown);

if ((!defined $blastx) || !defined ($m_unknown)) { 
  die $usage;
}

if ((!-e $blastx) || (!-e $m_unknown)) {
  die "File doesn't exist!\n";
}

my $out_blast = Bio::SearchIO->new ( -format =>'blast', -file => "$blastx");

open OUT1, ">identified_elements.txt" or die "Can't open identified_elements file for writing\n";
open OUT2, ">unknown_elements.txt" or die "Can't open unknown_elements file for writing\n";


my (%id_hash, %un_hash, %seq_hash, %query_hash, $seqobj); 

while (my $result = $out_blast ->next_result()) { 
  my $query_name = $result->query_name();
  my $query_desc = $result->query_description();
  my $num_hits1 = $result->num_hits();
  if ($num_hits1 == 0) {
      $un_hash{$query_name} = 1;
      next;
  }
  if ($num_hits1 > 0) {
    while (my $hit = $result->next_hit()) {
      my $sig = $hit->significance();
      if ($sig <= 1e-010) {	
	my $hit_name = $hit->name();
	$query_hash{$query_name} = $hit_name;
	last;
      }
      if ($sig > 1e-010) {	
	$un_hash{$query_name} = 1;
	last;
      }
    }
  }
}

my $in = Bio::SeqIO->new ( -format => 'fasta', -file => $m_unknown); 

my ($key1, $key2, $key, $key3, $key4);

while ($seqobj = $in ->next_seq()) {
  my $id = $seqobj->display_id();
  my $seq = $seqobj->seq();
  $seq_hash{$id} = $seq;
}

foreach $key (keys %query_hash) {
  foreach $key2 (keys %seq_hash) {
    if ($key eq $key2) {
      print OUT1 ">$query_hash{$key}\n";
      print OUT1 "$seq_hash{$key}\n";
    }
  }
}

foreach $key3 (keys %un_hash) {
  foreach $key4 (keys %seq_hash) {
    if ($key3 eq $key4) {
      print OUT2 ">$key3\n";
      print OUT2 "$seq_hash{$key3}\n";
    }
  }
}








	
