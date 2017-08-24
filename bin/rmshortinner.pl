#! /usr/bin/perl
 
$usage = "rmshortinner.pl  seqfile.outinner99.unmasked minimum_length_inner (suggested value 50 bp)\n";

if (@ARGV < 2) {die "$usage";}

if ($ARGV[1] <= 0) {die "minimum_length_inner must be more than zero!";}

open(FASTA, "$ARGV[0]") || die "seqfile.outinner99.unmasked does not exist!";

$seq = "";
while (<FASTA>){
  if (/>(\S+)\s*/) {
      if ($seq) {
	  @sym = split(//, $seq);
	  $ct = 0;
	  foreach $sym (@sym) {
	      $ct ++;
      }
      }
      $len{$name} = $ct;
      $name = $1;
      $seq = "";
  } else {
      chomp;
      $seq .= $_;
  }
}
close FASTA;

open(FASTA, "$ARGV[0]") || die "seqfile.outinner99.unmasked does not exist!";
while (<FASTA>) {
  
if   (/^>(\S+)\s+/)  {
if (&comparison){$take = 1;}
    else {$take = 0;}
}
              if ($take){
		  print;
}
}
close FASTA;
                                                                                
sub comparison  {
    foreach $key (keys %len){
  if ($key eq $1 && $len{$key} >= $ARGV[1])
  { return 1;}
}
}

