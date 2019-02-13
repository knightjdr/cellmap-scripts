package Interactions::Saint;

use strict;
use warnings;
 
use Exporter qw(import);
use List::Util qw(sum);
use Text::CSV_XS;
 
our @EXPORT_OK = qw(readSaintInteractions);

# Reads a SAINT file and parse interactions at a given FDR. For each bait-prey
# pair it will return the AvgSpec and PreySequenceLength. Use the third ($controlSubtract)
# argument to specify if control subtraction should be performed (0: no, 1: yes).
sub readSaintInteractions {
  my ($file, $cutoff, $controlSubtract) = @_;

  print STDERR "Reading SAINT file\n";
  my %interactions;
  my $tsv = Text::CSV_XS->new({
    sep_char => "\t"
  });
  open my $fh, '<', $file or die "Could not open $file: $!";
  $tsv->getline($fh); #discard header
  while(my $row = $tsv->getline($fh)) {
    my $fdr = @{$row}[15];
    my $avgSpec = @{$row}[5];
    if ($controlSubtract) {
      my @controlArray = split '\|', @{$row}[7];
      my $controlAvg = sum(@controlArray) / @controlArray;
      $avgSpec -= $controlAvg;
    }
	  my $bait = @{$row}[0];
    my $length = @{$row}[20];
    my $prey = @{$row}[2];
    if ($fdr <= $cutoff) {
      $interactions{$bait}{$prey}{'avgspec'} = $avgSpec;
      $interactions{$bait}{$prey}{'length'} = $length;
    }
  }
  close $fh;

  return \%interactions;
}
 