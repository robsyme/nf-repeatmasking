#!/usr/bin/perl -w
 
# CRL_Step4.pl
# This is the first step for creating exemplar elements. 

# Megan Bowman
# 31 January 2014

use strict;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;

my $usage = "\nUsage: $0 --step3 <path to CRL_Step3_Passed_Elements file> --resultfile <path to seqfile.result file> --innerfile <path to seqfile.outinner file> --sequencefile <path to whole genome sequence file>\n";

my ($resultfile, $sequencefile, $innerfile, $step3);

GetOptions('step3=s' => \$step3,
	   'resultfile=s' => \$resultfile,
	   'innerfile=s' => \$innerfile,
	   'sequencefile=s' => \$sequencefile);


if (!defined ($resultfile) || !defined ($innerfile) || !defined ($step3) || !defined ($sequencefile)) {
    die $usage;
}


if (!-e $resultfile) {
    die "$resultfile does not exist!\n"
}

if (!-e $innerfile) {
    die "$innerfile does not exist!\n"
}

if (!-e $step3) {
    die "$step3 does not exist!\n"
}

if (!-e $sequencefile) {
    die "$sequencefile does not exist!\n"
}


my $inner = Bio::SeqIO->new ( -format => 'fasta', -file => $innerfile);
my $step3in = Bio::SeqIO->new ( -format => 'fasta', -file => $step3);

my (%lLTRstart, %lLTRstop, %inner_id_hash, %inner_seq_hash, %inner_match_hash, %ltr_match_hash);



#-----------------------Create new sequence IDs based on the inner coordinates--------------------------------#

open OUT1, ">Inner_Seq_For_BLAST.fasta" or die "Failed to create the new OUT1 file\n";

open IN1, "$resultfile" or die "Failed to load result file : $!\n";


while (my $line2 = <IN1>) {
    chomp $line2;
    if ($line2 =~ /^#/) {
        next;     
    }
    $line2 =~ s/\s+/ /g;
    my @elems = split " ", $line2;  
    my $new_info = $elems[4] + 1;
    my $ID1 = $elems[10] . "_" . $elems[3];
    my $ID2 = $elems[10] . "_" ."$new_info";
    $inner_id_hash{$ID2} = $ID1;
    $lLTRstart{$ID1} = $elems[3];
    $lLTRstop{$ID1} = $elems[4];
}


my ($step3_seq, $step3_id, $step3_desc, $seq_id, $nr_index, $ltr_start, $inner_key, %result_seq_hash, $inner_start, $silly_key);

while (my $seqobj2 = $step3in->next_seq()) { 
    $step3_seq = $seqobj2->seq();
    $step3_id = $seqobj2->display_id();
    $step3_desc = $seqobj2->desc();
    if ($step3_id =~ /^(.+)_\(/) {
	$seq_id = $1;
    }
    if ($step3_id =~ /\(dbseq-nr_(\d+)\)_\[/) {
	$nr_index = $1;
    }
    if ($step3_id =~ /\[(\d+),/) {
	$ltr_start = $1;
    }
    $silly_key = $nr_index . "_" . $ltr_start;
    $result_seq_hash{$silly_key} = $seq_id;
}	

#-----------------------------Rename the sequence IDs for the inner sequences---------------------------------------#


while (my $seqobj1 = $inner->next_seq()) { 
    my $innerseq = $seqobj1->seq();
    my $innerid = $seqobj1->display_id();
    my $innerdesc = $seqobj1->desc();
    my $new_desc = $innerdesc; 
    $new_desc =~ s/ /_/g; 
    my $new_id = "$innerid" . "_" . "$new_desc";
    my $blank_desc = " "; 
    $innerdesc = $seqobj1->desc($blank_desc);
    $innerid = $seqobj1->display_id($new_id);
    if ($innerid =~ /^(.+)_\(/) {
	$seq_id = $1;
    }
    if ($innerid =~ /\(dbseq-nr_(\d+)\)_\[/) {
	$nr_index = $1;
    }
    if ($innerid =~ /\[(\d+),/) {
	$inner_start = $1;
    }
    $inner_key = $nr_index . "_" . $inner_start;
    $inner_seq_hash{$inner_key} = $innerseq;
}	

foreach my $key3 (keys %inner_id_hash) {
    foreach my $key4 (keys %inner_seq_hash) { 
	if ($key3 =~ /$key4/) { 
	    $inner_match_hash{$inner_id_hash{$key3}} = $inner_seq_hash{$key4};
	}
    }
}


#-------------------------Subseq the lLTR from the WGS and write to a new LTR file-----------------------------------#

my $seqin = Bio::SeqIO->new ( -format => 'fasta', -file => $sequencefile);

my (%sequencehash, $id, $lLTR);

while (my $seqobj3 = $seqin->next_seq()) {
    $id = $seqobj3->display_id();
    $sequencehash{$id} = $seqobj3;
}

open OUT2, ">lLTRs_Seq_For_BLAST.fasta";

  foreach my $key6 (keys %sequencehash) {
  foreach my $key7 (keys %result_seq_hash) {
      if ($key6 eq $result_seq_hash{$key7}) {
	  my $lLTR = $sequencehash{$key6}->subseq($lLTRstart{$key7}, $lLTRstop{$key7}); 
	  $ltr_match_hash{$key7} = $lLTR;
      }
  }
}


#-------------------------Create matching inner and LTR sequence files for BLAST ----------------------------------#

foreach my $key8 (keys %inner_match_hash) {
    foreach my $key9 (keys %ltr_match_hash) {
	if ($key8 eq $key9) {
	    print OUT1 ">$key8\n";
	    print OUT1 "$inner_match_hash{$key8}\n";
	    print OUT2 ">$key9\n";
	    print OUT2 "$ltr_match_hash{$key9}\n";
	}
    }
} 




