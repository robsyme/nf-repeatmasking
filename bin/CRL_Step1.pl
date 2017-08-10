#!/usr/bin/perl -w
 
# CRL_Step1.pl
# This script detects primer binding site (PBB) and RR_tracts from the outputs of LTR_Harvest and LTR_Digest. 

# Megan Bowman
# 13 January 2014

use strict;
use Getopt::Long();
use Bio::Seq;
use Bio::SeqIO;

my $usage = "\nUsage: $0 --gff <path to seqfile.gff.dgt file>\n";

my ($gff_file);

Getopt::Long::GetOptions('gff=s' => \$gff_file);

if (!defined ($gff_file)) { 
    die $usage;
}

if (!-e $gff_file) {
    die "$gff_file doesn't exist!\n"
}


#---------------------------Create files with or without specific protein domains-------------------------#

open IN, "$gff_file" or die "Can't open GFF3 file!\n";
open OUT1, ">NoPNoR.txt";
open OUT2, ">NoPYR.txt";
open OUT3, ">YPNoR.txt";
open OUT4, ">YPYR.txt";

my $primer_count = 0;
my $RR_count = 0;

my (@elems, $id, $info, $coord, $silly_index, $new, $data);

while (my $line = <IN>) {
    if ($line =~ /^#[^#]/) {
	next;
    }
    @elems = split "\t", $line;
    if ($elems[0] =~ /^seq/) {	
	$elems[0] =~ /^seq(\d+)$/;
	$info = $1; 
	if ($elems[2] =~ /LTR_retrotransposon/) {
	    $coord = $elems[3];
	    $silly_index = $info . "_" . $coord;
	}
	if ($elems[2] =~ /primer_binding_site/) {
	    $primer_count++;
	}
	if ($elems[2] =~ /RR_tract/) {
	    $RR_count++;
	}
	$data .= $line;
	next;
    }
    if ($elems[0] =~ /#{3}/) {
	$new = $primer_count ."_" . $RR_count;
	$primer_count = 0;
	$RR_count = 0;   
	if ($new =~ /0_0/) {
	    print OUT1 "$silly_index\n";
	    print OUT1"$data";
	    $data = "";
	    next;
	}
	if ($new =~ /0_1/) {
	    print OUT2 "$silly_index\n";
	    print OUT2 "$data";
	    $data = "";
	    next;
	}
	if ($new =~ /1_0/) {
	    print OUT3 "$silly_index\n";
	    print OUT3 "$data";
	    $data = "";
	    next;
	}
	if ($new =~ /1_1/) {
	    print OUT4 "$silly_index\n";
	    print OUT4 "$data";
	    $data = "";
	    next;
	}
    }        
}    
    
close OUT1;
close OUT2;
close OUT3;
close OUT4;

system ("rm NoPNoR.txt");

open IN1, "NoPYR.txt" or die "Can't open the No_Primer_Yes_RR file\n";

open OUT6, ">CRL_Step1_Passed_Elements.txt" or die "Didn't create the Final_Element File\n";

my $check1 = 0;

#-----------------------------Filter file without PPB but with RR_Tracts-----------------------------------#


my (%LTR1_start, %LTR1_end, %LTR1_S, %LTR1_E, %LTR2_start, %LTR2_S, %RR1_start, %RR1_end, %RR2_S, %RR2_E, %PBS1_start, %PBS1_end, %PBS2_S, %PBS2_E, $id1, $id2, $id3, @elems1, @elems2, @elems3);

while (my $line1 = <IN1>) {
    my @elems1 = split "\t", $line1;
    if ($elems1[0] =~ /^\d+_\d+$/) {
	$id1 = $elems1[0];
	next;
    }
    if ($elems1[2] =~ /long_terminal_repeat/) {
	if ($check1 == 0) {
	    $LTR1_start{$id1} = $elems1[3];
	    $LTR1_end{$id1} = $elems1[4];
	    $check1 = 1;
	    next;
	}
	if ($check1 == 1) {
	    $LTR2_start{$id1} = $elems1[3];
	    $check1 = 0;
	    next;
	}      
    }  
    if ($elems1[2] =~ /RR_tract/) {
	$RR1_start{$id1} = $elems1[3];
	$RR1_end{$id1} = $elems1[4];
	next;
    }
}    

foreach my $key1 (keys %RR1_start) {
    if ($RR1_start{$key1} - $LTR1_start{$key1} > $LTR2_start{$key1} - $RR1_start{$key1}) {
	if ($LTR2_start{$key1} - $RR1_end{$key1} <= 20) { 
	    if (($LTR2_start{$key1} -$RR1_start{$key1}) /($RR1_end{$key1} - $RR1_start{$key1}) * 100 >= 50) {
		print OUT6 "$key1\n";
		next;
	    }
	     if (($LTR2_start{$key1} -$RR1_start{$key1}) /($RR1_end{$key1} - $RR1_start{$key1}) * 100 < 50) { 
		 next;
	    }
	}
	if ($LTR2_start{$key1} - $RR1_end{$key1} > 20) {
	    next;
	}
    }
    if ($RR1_start{$key1} - $LTR1_start{$key1} < $LTR2_start{$key1} - $RR1_start{$key1}) {
	if ($RR1_start{$key1} - $LTR1_end{$key1} <= 20) { 
	    if (($RR1_end{$key1} -$LTR1_end{$key1}) /($RR1_end{$key1} - $RR1_start{$key1}) * 100 >= 50) { 
		print OUT6 "$key1\n";
		next;
	    }
	    if (($RR1_end{$key1} -$LTR1_end{$key1}) /($RR1_end{$key1} - $RR1_start{$key1}) * 100 < 50) { 
		next;
	    }
	}
	if ($RR1_start{$key1} - $LTR1_end{$key1} > 20) {
	    next;
	}
    }    
}
    
%LTR1_start = ();
%LTR1_end = ();
%LTR2_start = ();

#--------------------------------Filter file with PPB but without RR_Tracts---------------------------------#

open IN2, "YPNoR.txt" or die "Can't open the Yes_Primer_No_RR file\n";

my $check2 = 0;

while (my $line2 = <IN2>) {
    @elems2 = split "\t", $line2;
    if ($elems2[0] =~ /^\d+_\d+$/) {
	$id2 = $elems2[0];
	next;
    }
    if ($elems2[2] =~ /long_terminal_repeat/) {
      if ($check2 == 0) {
	  $LTR1_start{$id2} = $elems2[3];
	  $LTR1_end{$id2} = $elems2[4];
	  $check2 = 1;
	  next;
      }
      if ($check2 == 1) {
	  $LTR2_start{$id2} = $elems2[3];
	  $check2 = 0;
	  next;
      }      
  }  
  if ($elems2[2] =~ /primer_binding_site/) {
      $PBS1_start{$id2} = $elems2[3];
      $PBS1_end{$id2} = $elems2[4];
      next;
  }
}    
 
foreach my $key2 (keys %LTR1_start) {
    if ($PBS1_start{$key2} - $LTR1_start{$key2} > $LTR2_start{$key2} - $PBS1_start{$key2}) {
	if ($LTR2_start{$key2} - $PBS1_end{$key2} <= 20) { 
	    if (($LTR2_start{$key2} -$PBS1_start{$key2}) /($PBS1_end{$key2} - $PBS1_start{$key2}) * 100 >= 50) { 
		print OUT6 "$key2\n";
		next;
	    }
	     if (($LTR2_start{$key2} -$PBS1_start{$key2}) /($PBS1_end{$key2} - $PBS1_start{$key2}) * 100 < 50) { 
		 next;
	     }
	}
	if ($LTR2_start{$key2} - $PBS1_end{$key2} > 20) {
	    next;
	}
    }
    if ($PBS1_start{$key2} - $LTR1_start{$key2} < $LTR2_start{$key2} - $PBS1_start{$key2}) {
	if ($PBS1_start{$key2} - $LTR1_end{$key2} <= 20) { 
	    if (($PBS1_end{$key2} -$LTR1_end{$key2}) /($PBS1_end{$key2} - $PBS1_start{$key2}) * 100 >= 50) { 
		print OUT6 "$key2\n";
		next;
	    }
	    if (($PBS1_end{$key2} -$LTR1_end{$key2}) /($PBS1_end{$key2} - $PBS1_start{$key2}) * 100 < 50) {
		next;
	    }
	}
	if ($PBS1_start{$key2} - $LTR1_end{$key2} > 20) {
	    next;
	}
    }    
}

%LTR1_start = ();
%LTR1_end = ();
%LTR2_start = ();


#---------------------Filter file with PPB and with RR_Tracts----------------------------------#

open IN3, "YPYR.txt" or die "Can't open the Yes_Primer_Yes_RR file\n";

my $check3 = 0;

while (my $line3 = <IN3>) {
    @elems3 = split "\t", $line3;
    if ($elems3[0] =~ /^\d+_\d+$/) {
	$id3 = $elems3[0];
	next;
    }
    if ($elems3[2] =~ /long_terminal_repeat/) {
	if ($check3 == 0) {
	    $LTR1_S{$id3} = $elems3[3];
	    $LTR1_E{$id3} = $elems3[4];
	    $check3 = 1;
	    next;
	}
	if ($check3 == 1) {
	    $LTR2_S{$id3} = $elems3[3];
	    $check3 = 0;
	    next;
	}      
    }  
    if ($elems3[2] =~ /primer_binding_site/) {
	$PBS2_S{$id3} = $elems3[3];
	$PBS2_E{$id3} = $elems3[4];
	next;
    }
    if ($elems3[2] =~ /RR_tract/) {
	$RR2_S{$id3} = $elems3[3];
	$RR2_E{$id3} = $elems3[4];
	next;
    }
}  

my $key3;

foreach $key3 (keys %LTR1_S) {
    if (($PBS2_S{$key3} - $LTR1_S{$key3} < $LTR2_S{$key3} - $PBS2_S{$key3}) && ($RR2_S{$key3} - $LTR1_S{$key3} > $LTR2_S{$key3} - $RR2_S{$key3})) {
	if (($PBS2_S{$key3} - $LTR1_E{$key3}) > 20 || ($LTR2_S{$key3} - $RR2_E{$key3}) > 20) {
	    next;
	}
	if (($PBS2_S{$key3} - $LTR1_E{$key3} <= 20) && ($LTR2_S{$key3} - $RR2_E{$key3}) <= 20) {	    
	    if ((($LTR2_S{$key3} -$RR2_S{$key3}) /($RR2_E{$key3} - $RR2_S{$key3})) * 100 >= 50 && (($PBS2_E{$key3} -$LTR1_E{$key3}) /($PBS2_E{$key3} - $PBS2_S{$key3})) * 100 >= 50) { 
		print OUT6 "$key3\n";
		next;
	    }
	    if ((($LTR2_S{$key3} -$RR2_S{$key3}) /($RR2_E{$key3} - $RR2_S{$key3})) * 100 < 50 || (($PBS2_E{$key3} -$LTR1_E{$key3}) /($PBS2_E{$key3} - $PBS2_S{$key3})) * 100 < 50) {
		next;
	    }
	}
    }
    if (($PBS2_S{$key3} - $LTR1_S{$key3} > $LTR2_S{$key3} - $PBS2_S{$key3}) && ($RR2_S{$key3} - $LTR1_S{$key3} < $LTR2_S{$key3} - $RR2_S{$key3})) {
	if (($LTR2_S{$key3} - $PBS2_E{$key3}) > 20 || ($RR2_S{$key3} - $LTR1_E{$key3}) > 20) {
	    next;
	}
	if (($LTR2_S{$key3} - $PBS2_E{$key3}) <= 20 && ($RR2_S{$key3} - $LTR1_E{$key3}) <= 20) {
	    if ((($LTR2_S{$key3} -$PBS2_S{$key3}) /($PBS2_E{$key3} - $PBS2_S{$key3})) * 100 >= 50 && ($RR2_E{$key3} -$LTR1_E{$key3}) /($RR2_E{$key3} - $RR2_S{$key3}) * 100 >= 50) {
		print OUT6 "$key3\n";
		next;
	    }
	    if ((($LTR2_S{$key3} -$PBS2_E{$key3}) /($PBS2_E{$key3} - $PBS2_S{$key3})) * 100 < 50 || ($RR2_E{$key3} -$LTR1_E{$key3}) /($RR2_E{$key3} - $RR2_S{$key3}) * 100 < 50) { 
		next;
	    }
	}
    } 
}

close OUT6;	 

system ("rm NoPYR.txt YPNoR.txt YPYR.txt");

	    



