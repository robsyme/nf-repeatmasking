#!/usr/bin/perl -w
 
# CRL_Step2.pl
# This script will take the seqfile.out file, remove elements containing 50 or more Ns, identify specific genomic regions containing LTR repeats in a whole genome sequence, find 50 bp upstream and downstream of each element within the genome sequence, and write it to a new fasta file. 


# Megan Bowman
# 18 December 2013

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;

my $usage = "\nUsage: $0 --step1 <path to step1 passed elements> --repeatfile <path to seqfile.out file> --resultfile <path to seqfile.result file> --sequencefile <path to whole genome sequence> --removed_repeats <name of removed repeat file>\n";

my($removed_repeats, $resultfile, $sequencefile, $repeatfile, $step1);

GetOptions('resultfile=s' => \$resultfile,
	   'repeatfile=s' => \$repeatfile,
	   'sequencefile=s' => \$sequencefile,
	   'removed_repeats=s' => \$removed_repeats,
	   'step1=s' => \$step1);



if (!defined ($resultfile) || !defined ($sequencefile) || !defined ($repeatfile) || !defined ($removed_repeats)) {
    die $usage;
}

if (!-e $repeatfile) {
    die "$repeatfile does not exist!\n"
}

if (!-e $resultfile) {
    die "$resultfile does not exist!\n"
}

if (!-e $sequencefile) {
    die "$sequencefile does not exist!\n"
}

if (-e $removed_repeats) {
    die "$removed_repeats already exists!\n"
}

if (!-e $step1) {
    die "$step1 does not exist!\n"
}


my (@elems, %datahashlLTRup, %datahashrLTRup, %datahashlLTRdown, %datahashrLTRdown);

my $repin = Bio::SeqIO->new ( -format => 'fasta', -file => $repeatfile); 
my $seqin = Bio::SeqIO->new ( -format => 'fasta', -file => $sequencefile);
my $seqout = Bio::SeqIO->new (-format => 'fasta', -file => ">$removed_repeats");

#---------------------------Gather IDs that passed LTR_Digest Filtering--------------------------------------------------#

open IN1, "$step1" or die "Couldn't open the Step1_Passed_Elements file\n"; 

my (@gff_array, $gff);

while ($gff = <IN1>) {
  chomp $gff; 
  push (@gff_array, $gff)
}

push (@gff_array, $gff);

#---------------------------Remove the LTRs with 50 or more N's-----------------------------------------------------------#

my $new_id;

while (my $seqobj1 = $repin->next_seq()) {
    my $repseq = $seqobj1->seq();
    my $n_count = 0;
    while ($repseq =~ /([Nn]+)/g) { 
	$n_count += length($1);
	if ($n_count >= 50) {
	    last;
	}
    }
    if ($n_count < 50) {
	my $repid = $seqobj1->display_id();
	my $repdesc = $seqobj1->desc();
	my $new_desc = $repdesc; 
	$new_desc =~ s/ /_/g;
	my $new_id = "$repid" . "_" . "$new_desc";
	my $blank_desc = " "; 
	$repdesc = $seqobj1->desc($blank_desc);
	$repid = $seqobj1->display_id($new_id);
	$seqout->write_seq($seqobj1);
    }
}

#--------------------------Find the sequence ID and LTR start position for each LTR-----------------------------------------#

my @art_array = ();
my @best_array = ();

my ($silly_index, $ltr_start, $seq_id, $remid, $artificial_key, %artificial_key_hash, $i2, $key3);

my $removed = Bio::SeqIO->new (-format => 'fasta', -file => $removed_repeats);

while (my $seqobj2 = $removed->next_seq()) {    
  $remid = $seqobj2->display_id();
  if ($remid =~ /^(.+)_\(/) {
    $seq_id = $1;
  }
  else {
  }
  if ($remid =~ /\(dbseq-nr_(\d+)\)_\[/) {
    $silly_index = $1;
  }
  if ($remid =~ /\[(\d+),/) {
    $ltr_start = $1;
  }

  $artificial_key = $silly_index . "_" . $ltr_start;
  push (@art_array, $artificial_key);
  $artificial_key_hash{$artificial_key} = $seq_id; 
}

push (@art_array, $artificial_key);


#--------------------Compare artificial_keys with those that passed LTR_Digest filtering-----------------------------------------#

foreach $key3 (keys %artificial_key_hash) { 
  for ($i2 = 0; $i2 < $#gff_array; $i2++) {
    if ($key3 eq $gff_array[$i2]) {
      push (@best_array, $key3);
      next;
    }
  }
}
push (@best_array, $key3);



#------------------Create new data structures with the sequence ids and start points for the lLTR and rLTR------------------------#

open IN2, "$resultfile" or die "Failed to load result file : $!\n";
while (my $line2 = <IN2>) {
    chomp $line2;
    if ($line2 =~ /^#/) {
	next;     
    }
    $line2 =~ s/\s+/ /g;
    my @elems = split " ", $line2;     
    my $key = $elems[10] . "_" . $elems[0];
    if (!exists($artificial_key_hash{$key})) {
	next;
    }
    $datahashlLTRup{$key} = $elems[3];
    $datahashrLTRup{$key} = $elems[6];
    $datahashlLTRdown{$key} = $elems[4];
    $datahashrLTRdown{$key} = $elems[7];
}

#--------------------Create new fasta files containing 50bp upstream or downstream of each lLTR and rLTR---------------------------#



my ($key2, $seqobj, $id, $i, %sequencehash, $length, $check, $lLTRup50, $rLTRup50, $lLTRdown50, $rLTRdown50, %lengthhash, %dump_hash);

$check = 1;

while ($seqobj = $seqin->next_seq()) {
  $id = $seqobj->display_id();
  $length = $seqobj->length();
  $sequencehash{$id} = $seqobj;
  $lengthhash{$id} = $length;
}

foreach $key2 (keys %sequencehash) { 
  for ($i= 0; $i < $#best_array; $i++) {
    if ($key2 eq $artificial_key_hash{$best_array[$i]}) {
      if (!exists $dump_hash{$artificial_key_hash{$best_array[$i]}}) {	
	if ($datahashlLTRup{$best_array[$i]}-50 < 0) {
	  $dump_hash{$check} = $artificial_key_hash{$best_array[$i]};
	  next;
	}
	if ($datahashrLTRdown{$best_array[$i]}+50-$lengthhash{$key2} > 0) {
	  $dump_hash{$check} = $artificial_key_hash{$best_array[$i]};
	  next;
	}
	open OUT2, ">Repeat_up$check.fasta";
	open OUT3, ">Repeat_down$check.fasta";
	$lLTRup50 = $sequencehash{$key2}->subseq($datahashlLTRup{$best_array[$i]}-49, $datahashlLTRup{$best_array[$i]}); 
	$rLTRup50 = $sequencehash{$key2}->subseq($datahashrLTRup{$best_array[$i]}-49, $datahashrLTRup{$best_array[$i]});
	$lLTRdown50 = $sequencehash{$key2}->subseq($datahashlLTRdown{$best_array[$i]}, $datahashlLTRdown{$best_array[$i]}+49); 
	$rLTRdown50 = $sequencehash{$key2}->subseq($datahashrLTRdown{$best_array[$i]}, $datahashrLTRdown{$best_array[$i]}+49);
	print OUT2 ">$artificial_key_hash{$best_array[$i]}" ."_". "$best_array[$i]" . "_" . "lLTR50bp" . "_" . "upstream" . "_" ."$check\n";
	print OUT2 "$lLTRup50\n";
	print OUT2 ">$artificial_key_hash{$best_array[$i]}" ."_". "$best_array[$i]" . "_" . "rLTR50bp" . "_" . "upstream" ."_" . "$check\n";
	print OUT2 "$rLTRup50\n";
	close OUT2;
	print OUT3 ">$artificial_key_hash{$best_array[$i]}" ."_". "$best_array[$i]" . "_" . "lLTR50bp" . "_" . "downstream" . "_" .  "$check\n";
	print OUT3 "$lLTRdown50\n";
	print OUT3 ">$artificial_key_hash{$best_array[$i]}" ."_". "$best_array[$i]" . "_" . "rLTR50bp" . "_" . "downstream" .  "_" . "$check\n";
	print OUT3 "$rLTRdown50\n";
	close OUT3;
	$dump_hash{$check} = $artificial_key_hash{$best_array[$i]};
	++$check;
      }
      else {
	next;
      }
    }
  }
}



