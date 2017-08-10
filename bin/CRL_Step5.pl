#!/usr/bin/perl -w
 
# CRL_Step5.pl
# This is the second step and creates the final set of exemplar elements. 


# Megan Bowman
# 31 January 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my $usage = "\nUsage: $0 --LTR_blast <path of the LTR blast result File> --inner_blast <path of the inner blast result File> --step3 <path to CRL_Step3_passed elements file> --final <name of final output file> --pcoverage <percent coverage> --pidentity <percent identity>\n";

my ($LTR_blast, $inner_blast, $final, $pcoverage, $pidentity, $step3);

GetOptions('LTR_blast=s' => \$LTR_blast,
	   'inner_blast=s' => \$inner_blast,
	   'final=s' => \$final,
	   'step3=s' => \$step3,
	   'pcoverage=i' => \$pcoverage,
	   'pidentity=i' => \$pidentity);

if (!defined ($inner_blast) || !defined ($LTR_blast) || !defined ($final) || !defined ($step3)) {
    die $usage;
}

if (!-e $LTR_blast) {
    die "$LTR_blast does not exist!\n"
}

if (!-e $inner_blast) {
    die "$inner_blast does not exist!\n"
}

if (-e $final) {
  die "$final already exists!\n"
}

if (!-e $step3) {
    die "$step3 does not exist!\n"
}

if (!defined ($pcoverage) || !defined ($pidentity)) {
  $pcoverage = 90;
  $pidentity = 80;
}


my $step3in = Bio::SeqIO->new ( -format => 'fasta', -file => "$step3"); 
my $seqout = Bio::SeqIO->new (-format => 'fasta', -file => ">$final");
my $out_blast = Bio::SearchIO->new ( -format =>'blast', -file => "$inner_blast");

my (%inner_blast_hash, @iden1, @hsp_indiv1, $total_hsp1, $hsp_ind_len1, $query_length1);


#-----------------------------Create hashes containing the BLAST results--------------------------------------#

while (my $result = $out_blast ->next_result()) { 
  my @hit_array1;
  my $query_name1 = $result->query_name();
  my $query_length1 = $result->query_length();
  my $num_hits1 = $result->num_hits();
  if ($num_hits1 == 0) {
    last;
  }
  while (my $hit = $result->next_hit()) {
    my $hit_name = $hit->name();
    @iden1 = ();
    @hsp_indiv1 = ();
    $total_hsp1 = 0;
    while (my $hsp = $hit->next_hsp) {
      my $frac_identical = $hsp->frac_identical();
      push(@iden1, $frac_identical);
      $hsp_ind_len1 = $hsp->length('query');
      $total_hsp1 += $hsp_ind_len1;
      push(@hsp_indiv1, $hsp_ind_len1);
    }
    my $percent_coverage = ($total_hsp1/$query_length1) * 100;	
    my $percent_identity = 0;
    for my $i1 (0..$#hsp_indiv1) {
      $percent_identity += (($hsp_indiv1[$i1]/$total_hsp1) * $iden1[$i1]) * 100;  
    }
    if ($percent_identity >= $pidentity && $percent_coverage >= $pcoverage) {
      push (@hit_array1, $hit_name);
    }
  }
  $inner_blast_hash{$query_name1} = \@hit_array1;
}

my $LTR_blast_result = Bio::SearchIO->new ( -format =>'blast', -file => "$LTR_blast");

my (%LTR_blast_hash, @iden2, @hsp_indiv2, $total_hsp2, $hsp_ind_len2, $query_length2, @biggest_key);

while (my $result = $LTR_blast_result ->next_result()) { 
  my @hit_array2;
  my $query_name2 = $result->query_name();
  my $query_length2 = $result->query_length();
  my $num_hits2 = $result->num_hits();
  if ($num_hits2 == 0) {
    last;
  }
  while (my $hit = $result->next_hit()) {
    my $hit_name = $hit->name();
    @iden2 = ();
    @hsp_indiv2 = ();
    $total_hsp2 = 0;
    while (my $hsp = $hit->next_hsp) {
      my $frac_identical = $hsp->frac_identical();
      push(@iden2, $frac_identical);
      my $hsp_ind_len2 = $hsp->length('query');
      $total_hsp2 += $hsp_ind_len2;
      push(@hsp_indiv2, $hsp_ind_len2);
    }
    my $percent_coverage = ($total_hsp2/$query_length2) * 100;
    my $percent_identity = 0;  
    for my $i2 (0..$#hsp_indiv2) {
      $percent_identity += (($hsp_indiv2[$i2]/$total_hsp2) * $iden2[$i2]) * 100; 
    }
    if ($percent_identity >= $pidentity && $percent_coverage >= $pcoverage) {
      push (@hit_array2, $hit_name);
    }
  }
  $LTR_blast_hash{$query_name2} = \@hit_array2;	
}

#------------------------Identify exemplar sequences from BLAST results------------------------------------#

	 
my %exemplar_hash;
my %duplicated_elements;

while (scalar (keys %inner_blast_hash)) {  
  my @biggest_key = sort mysort (keys (%inner_blast_hash));
  my $biggest_key = $biggest_key[0];
  $exemplar_hash{$biggest_key} = 1;
  $duplicated_elements{$biggest_key} = 1;
  foreach my $key (@{$inner_blast_hash{$biggest_key}}) { 
    $duplicated_elements{$key} = 1;
    delete $inner_blast_hash{$key};
    delete $LTR_blast_hash{$key};
  }
  delete $inner_blast_hash{$biggest_key};
  delete $LTR_blast_hash{$biggest_key};
  foreach my $key (keys %inner_blast_hash) {
    foreach my $key2 (@{$inner_blast_hash{$key}}) {
      my @array;
      if (!exists $duplicated_elements{$key2}) { 
	push (@array, $key2);
      }
      $inner_blast_hash{$key} = \@array;
    }
    foreach my $key2 (@{$LTR_blast_hash{$key}}) {
      my @array;
      if (!exists $duplicated_elements{$key2}) { 
	push (@array, $key2);
      }
      $LTR_blast_hash{$key} = \@array;
    }
  }
}

while (scalar (keys %LTR_blast_hash)) {  
  my ($biggest_key) = sort { scalar (@{$LTR_blast_hash{$a}}) <=> scalar (@{$LTR_blast_hash{$b}}) } keys (%LTR_blast_hash);
  $exemplar_hash{$biggest_key} = 1;
  $duplicated_elements{$biggest_key} = 1;
  foreach my $key (@{$LTR_blast_hash{$biggest_key}}) { 
    $duplicated_elements{$key} = 1;
    delete $LTR_blast_hash{$key};
  }
  delete $LTR_blast_hash{$biggest_key};
  foreach my $key (keys %LTR_blast_hash) {
    foreach my $key2 (@{$LTR_blast_hash{$key}}) {
      my @array;
      if (!exists $duplicated_elements{$key2}) { 
	push (@array, $key2);
      }
      $LTR_blast_hash{$key} = \@array;
    }
  }
} 

#------------------------Create final fasta file of exemplar sequences------------------------------------#


my ($good_id, $step3_ltr_index, $step3_silly_index, %step3_seq_hash, $step3_ltr_start, $step3_silly_key, $step3_nr_index);

while (my $seqobj4 = $step3in->next_seq()) { 
  my $step3_seq = $seqobj4->seq();
  my $step3_id = $seqobj4->display_id();
  my $step3_desc = $seqobj4->desc();
  if ($step3_id =~ /^(.+)_\(/) {
    $good_id = $1;
  }
  if ($step3_id =~ /\(dbseq-nr_(\d+)\)_\[/) {
    $step3_nr_index = $1;
  }
  if ($step3_id =~ /\[(\d+),/) {
    $step3_ltr_start = $1;
  }
  $step3_silly_key = $step3_nr_index . "_" . $step3_ltr_start;
  $step3_seq_hash{$step3_silly_key} = $seqobj4;
}	

foreach my $key9 (sort keys %exemplar_hash) {
  foreach my $key10 (sort keys %step3_seq_hash) {
    if ($key10 eq $key9) { 
     $seqout->write_seq($step3_seq_hash{$key10});
   }
  }
}

sub mysort {
  scalar (@{$inner_blast_hash{$a}}) <=> scalar (@{$inner_blast_hash{$b}}); 
}
