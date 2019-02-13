package Bait::Gene;

use strict;
use warnings;
 
use Exporter qw(import);
use Text::CSV_XS;
 
our @EXPORT_OK = qw(parseMap);

# Reads a text file with a header and parses the gene names for each bait.
# The bait names should be in column one and the gene names in column two.
# It will return a map (baitMap) for bait to gene names, and a map (geneMap)
# of the gene name in lowercase to the actual Bait name.
sub parseMap {
  my ($file) = @_;

  print STDERR "Reading bait-gene map\n";

  my %baitMap;
  my %geneMap;

  my $tsv = Text::CSV_XS->new({
    sep_char => "\t"
  });
  open my $fh, '<', $file or die "Could not open $file: $!";
  $tsv->getline($fh); #discard header
  while(my $row = $tsv->getline($fh)) {
	  my $bait = @{$row}[0];
	  my $gene = @{$row}[1];
	  $baitMap{$bait} = $gene;
    $geneMap{lc $gene} = $bait;
  }
  close($fh);

  return (\%baitMap, \%geneMap);
}
 