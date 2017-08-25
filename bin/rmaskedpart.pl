#! /usr/bin/perl

$usage = "rmaskedpart.pl maskedfile number-of-positions-per-line\n";

# this is for removing the masked portion or sequencing gaps 

if (@ARGV < 2) {die $usage;}
if ($ARGV[1] < 1) {die $usage;}

open(FA, "$ARGV[0]") || die $usage;

$seq = "";
while (<FA>) {
    if (/>\s*(.+)/) {
	if ($seq) {
	    @sym = split(//, $seq);
	    $ct = 0;
	    foreach $sym (@sym) {
		if ($sym ne "N") {
		print $sym;
		$ct ++;
		}
		if ( !($ct%$ARGV[1]) ) {print "\n";}
	    }
	    if ($ct%$ARGV[1]) {print "\n";}
	}
	printf ">%s\n", $1;
	$seq = "";
    } else {
	chomp;
	$seq .= $_;
    }
}
close FA;

@sym = split(//, $seq);
$ct = 0;
foreach $sym (@sym) {
    if ($sym ne "N") {
    print $sym;
    $ct ++;
    }
    if ( !($ct%$ARGV[1]) ) {print "\n";}
}
if ($ct%$ARGV[1]) {print "\n";}
