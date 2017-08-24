#!/usr/bin/perl -w
 
# ltr_library.pl
# This script is create a fasta file containing the lLTRs of candidate LTRs for the creation of a LTR + Mite.lib to mask the internal sequences found in the .outinner file from LTRHarvest.
# Megan Bowman
# 26 December 2013 

use strict;
use Getopt::Long();
use Bio::Seq;
use Bio::SeqIO;

my $usage = "\nUsage: $0 --resultfile <seq.result file name> --step3 <step 3 passed elements> --sequencefile <whole genome sequence fasta file> \n";

my ($resultfile, $sequencefile, $seqin, $step3); 

Getopt::Long::GetOptions('resultfile=s' => \$resultfile,
			 'sequencefile=s' => \$sequencefile,
			 'step3=s' => \$step3);


if (!defined ($resultfile) || !defined ($sequencefile) || !defined ($step3)) {
    die $usage;
}

if (!-e $resultfile) {
    die "$resultfile does not exist!\n";
}

if (!-e $sequencefile) {
    die "$sequencefile does not exist!\n";
}

if (!-e $step3) {
    die "$step3 does not exist!\n";
}


#---------------Get IDs of candidate LTRs without 50Ns------------------------#

my ($nr_index, $ltr_start, $seq_id, $remid, $artificial_key, %artificial_key_hash, @art_array);

my $removed = Bio::SeqIO->new (-format => 'fasta', -file => $step3);

while (my $seqobj2 = $removed->next_seq()) {    
    $remid = $seqobj2->display_id();
    if ($remid =~ /^(.+)_\(/) { 
        $seq_id = $1;
    }
    if ($remid =~ /\(dbseq-nr_(\d+)\)_\[/) {
        $nr_index = $1;
    }
    if ($remid =~ /\[(\d+),/) {
        $ltr_start = $1;
    }
    $artificial_key = $nr_index . "_" . $ltr_start;
    push (@art_array, $artificial_key); 
    $artificial_key_hash{$artificial_key} = $seq_id; 
}


#----------Get coordinates of the lLTR for each candidate LTR------------------#

my (%lLTRstart, %lLTRstop);

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
    $lLTRstart{$key} = $elems[3];
    $lLTRstop{$key} = $elems[4];
}


#----------Subseq the lLTR from the WGS and write to a new LTR file------------------#

my ($id, $seqobj, $i, %sequencehash, $lLTR);

$seqin = Bio::SeqIO->new ( -format => 'fasta', -file => $sequencefile);

while ($seqobj = $seqin->next_seq()) {
    $id = $seqobj->display_id();
    $sequencehash{$id} = $seqobj;
}

open OUT, ">lLTR_Only.lib";

foreach my $key (keys %sequencehash) {
    for ($i= 0; $i < $#art_array; $i++) {       
	if ($artificial_key_hash{$art_array[$i]} eq $key) {
	    $lLTR = $sequencehash{$key}->subseq($lLTRstart{$art_array[$i]}, $lLTRstop{$art_array[$i]}); 
	    print OUT ">$id" . "_" . "$lLTRstart{$art_array[$i]}\n";
	    print OUT "$lLTR\n";
	}
    }        
}    

