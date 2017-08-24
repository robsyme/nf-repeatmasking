#! /usr/bin/perl
 
$usage = "cleanRM.pl  seqfile.outinner99.out  seqfile.outinner99.masked\n";

if (@ARGV < 2) {die "$usage";}

open(OUT, "$ARGV[0]") || die "seqfile.outinner99.out does not exist!";

while (<OUT>){
if   (/^\s*(\d+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+/)  {
    if ($1 > 300 && $2 ne $lseq) {
	$seq{$ct} = $2;
	$ct ++;
	$lseq = $2;
    }
 }
}
close OUT;

open(FASTA, "$ARGV[1]") || die "seqfile.outinner99.masked does not exist!";
while (<FASTA>) {
  
if   (/^>(\S+)\s+/)  {
if (&comparison){$take = 0;}
    else {$take = 1;}
}
              if ($take){
		  print;
}
}
close FASTA;
                                                                                
sub comparison  {
    foreach $key (keys %seq){
  if ($seq{$key} eq $1)
  { return 1;}
}
}

