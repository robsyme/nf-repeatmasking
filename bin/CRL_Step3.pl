#!/usr/bin/perl -w

# CRL_Step3.pl
# This script take the fasta files containing 50 bp upstream and downstream of each element, and aligns them with MUSCLE. Those 50bp alignments not passing the filters specified by the user will be removed from the pipeline. The output is a fasta file containing filtered elements for the next step of the pipeline. 

# Megan Bowman
# 21 January 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my ($directory, $id_1, $id_2, @seq1_array, @seq2_array, $step2, $pidentity, $seq_c, $seqobj1);


my $usage = "\nUsage: $0 --directory <full path of the directory with fasta files> --step2 <path to CRLStep2_Passed Elements file> --pidentity <percent identity> --seq_c <number of identical nucleotides>\n";

GetOptions('directory=s' => \$directory,
	   'step2=s'=> \$step2,
	   'pidentity=i' => \$pidentity,
	   'seq_c=i' => \$seq_c);

if (!defined ($directory) || !defined ($step2)) {
  die $usage;
}

if (!-e $directory) {
  die "$directory does not exist!\n"
}

if (!-e $step2) {
  die "$step2 does not exist!\n"
}

if (!defined ($pidentity) || !defined ($seq_c)) {
  $pidentity = 60;
  $seq_c = 25;
}

#---------------------------MUSCLE Alignments of 50bp up/downstream sequences---------------------------------------------------#
  
opendir DIR, $directory || die "\nUnable to open $directory for reading.\n\n";
 
my @files = readdir(DIR);

closedir DIR;

foreach my $file (@files) {
  if ($file eq ".." || $file eq "." ) { 
    next;
  }
  if ($file =~ /^Repeat/) { 
    my $open_me1 = $directory . "/" . $file;
    system ("muscle -in $open_me1 -out $file.afa");  
  }
}
#-------------------------------------Filtering of MUSCLE Alignments------------------------------------------------------------#


opendir DIR, $directory || die "\nUnable to open $directory for reading.\n\n";

my @files2 = readdir(DIR);

closedir DIR;

my $check = 1;

my (@down_array, %up_hash, $match_1, $match_2);

foreach my $file (@files2) {
  if ($file eq ".." || $file eq "." ) { 
    next;
  }  
  if ($file =~ /.fasta$/ ) { 
    next; 
  }
  if ($file =~ /.txt$/ ) { 
    next; 
  }
  if ($file =~ /.afa$/) {    
    my $open_me2 = $directory . "/" . $file;
  
    open OUT1, ">muscle_filtered.fasta" or die "Could not open the OUT file.\n";
    
    my $in = Bio::SeqIO->new ( -format => 'fasta', -file => "$open_me2"); 
    
    my $seq1;
  
    my $array_count = 0;
    
    my @id_array;
  
    while (my $seqobj1 = $in->next_seq()) {
      if ($array_count == 1) {
      $id_2 = $seqobj1->display_id();
      $seq1 = $seqobj1->seq();
	@seq2_array = split "" , $seq1;
	next;
    }
    $id_1 = $seqobj1->display_id(); 
    $seq1 = $seqobj1->seq();
    @seq1_array = split "" , $seq1;
    $array_count = 1;
    }
    
    my $seq_count = 0;
    my $gap_count = 0;
    
    for (my $i = 0; $i < $#seq1_array; $i++) {
      if ($seq1_array[$i] =~ /$seq2_array[$i]/) {
	$seq_count++;
	if ($seq1_array[$i] =~ /-/) {
	  $gap_count++;
	}
      }
    }
    
    my $percent_identity = (($seq_count-$gap_count)/$#seq1_array) * 100;
    
    if ($percent_identity < $pidentity && $seq_count < $seq_c) {
      if ($id_1 =~ /upstream/) {
	$up_hash{$id_1} = $id_1;
      }
      if ($id_1 =~ /downstream/) {
	push (@down_array, $id_1);
      }
    }
  }
}
push (@down_array, $id_1);

my $key1;

foreach $key1 (keys %up_hash) {
  if ($key1 =~ /^(.+_\d+_\d+)/) {
    my $match_1 = $1;
    for (my $i2 = 0; $i2 < $#down_array; $i2++) {
      if ($down_array[$i2] =~/^(.+_\d+_\d+)/) {
	$match_2 = $1;
    	if ($match_1 eq $match_2) { 
	  print OUT1 "$match_1\n";
	}
      } 
    }
  }
} 

close OUT1;

system ("rm *.afa");
#system ("rm Repeat_*");

#------------------------------------Creation of output file with filtered elements--------------------------------------------------#


open IN, "muscle_filtered.fasta" or die "Could not open final element list\n";
open OUT2, ">CRL_Step3_Passed_Elements.fasta" or die "Could not open the OUT2 file.\n";

my ($step2in, $seqobj2, $step2id, $seq_id, $silly_index, $artificial_key, $ltr_start, @id_array, %seq_hash, %seq_id_hash, $step2seq);

$step2in = Bio::SeqIO->new ( -format => 'fasta', -file => $step2); 

while (my $seqobj2 = $step2in->next_seq()) {
  $step2id = $seqobj2->display_id();
  $step2seq = $seqobj2->seq();
  if ($step2id =~ /^(.+)_\(/) {
    $seq_id = $1;
  }
  if ($step2id =~ /\(dbseq-nr_(\d+)\)_\[/) {
    $silly_index = $1;
  }
  if ($step2id =~ /\[(\d+),/) {
    $ltr_start = $1;
  }
  $artificial_key = $seq_id . "_" . $silly_index . "_" . $ltr_start;
  $seq_id_hash{$artificial_key} = $step2id;
  $seq_hash{$artificial_key} = $step2seq;
}

my $line;

while ($line = <IN>) {
  chomp $line;
  push (@id_array, $line);
}

push (@id_array, $line);

foreach my $key2 (keys %seq_id_hash) {
  for (my $i3 = 0; $i3 < $#id_array; $i3++) {
    if ($key2 eq $id_array[$i3]) {
      print OUT2 ">$seq_id_hash{$key2}\n";
      print OUT2 "$seq_hash{$key2}\n";
    }
  }
}

system ("rm muscle_filtered.fasta");


