#!/usr/bin/perl -w
 
# outinner_blastx_parse.pl
# This script will parse the blastx results of the seqfile.outinner sequence against the transposon database and output those without hits to be used for downsteam steps.

# Megan Bowman
# 26 October 2015

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my $usage = "\nUsage: $0 --blastx <blastx output file> --outinner <full path to seqfile.outinner>\n";

my ($blastx, $outinner);

GetOptions('blastx=s' => \$blastx,
	  'outinner=s' => \$outinner);

if ((!defined $blastx) || !defined ($outinner)) { 
  die $usage;
}

if ((!-e $blastx) || (!-e $outinner)) {
  die "File doesn't exist!\n";
}

my $out_blast = Bio::SearchIO->new ( -format =>'blast', -file => "$blastx");

open OUT1, ">passed_outinner_sequence.fasta" or die "Can't open passed_outinner_sequence file for writing\n";

my (%id_hash, %seq_hash, %query_hash, %new_id_hash, $seqobj); 

while (my $result = $out_blast ->next_result()) { 
  my $query_name = $result->query_name();
  my $query_desc = $result->query_description();
  my $num_hits1 = $result->num_hits();
  if ($num_hits1 == 0) {
      $query_hash{$query_desc} = 1;
      next;
  }
  if ($num_hits1 > 0) {
      next;
  }
}
  

my $in = Bio::SeqIO->new ( -format => 'fasta', -file => $outinner); 

my ($key1, $key2, $key, $key3, $key4);

while ($seqobj = $in ->next_seq()) {
  my $id = $seqobj->display_id();
  my $desc = $seqobj->desc();
  my $seq = $seqobj->seq();
  my $new_desc = $desc;
  $new_desc =~ s/ /_/g;
  my $new_id = "$id" . "_" . "$new_desc";
  $new_id_hash{$desc} = $new_id;
  $seq_hash{$desc} = $seq;
}

foreach $key (keys %query_hash) {
  foreach $key2 (keys %seq_hash) {
    if ($key eq $key2) {
	print OUT1 ">$new_id_hash{$key2}\n";
	print OUT1 "$seq_hash{$key2}\n";
    }
  }
}










	
