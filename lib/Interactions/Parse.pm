package Interactions::Parse;

use strict;
use warnings;
 
use Exporter qw(import);
use Text::CSV_XS;
 
our @EXPORT_OK = qw(readInteractions);

# Reads a text file with a header and parses known interactions for requested genes.
# The input should be of the format "source\ttarget\tmethod". Gene names will be converted
# to lowercase.
sub readInteractions {
  my $file = $_[0];
  my %genes = %{$_[1]};

  print STDERR "Getting list of interactors\n";
  my %interactions;

  my $tsv = Text::CSV_XS->new({
    sep_char => "\t"
  });
  open my $fh, '<', $file or die "Could not open $file: $!";
  $tsv->getline($fh); #discard header
  while(my $row = $tsv->getline($fh)) {
	  my $source = lc @{$row}[0];
	  my $target = lc @{$row}[1];
    if (exists $genes{$source}) {
      my %currHash = map { $_ => 1 } @{$interactions{$source}};
      if (!(exists $currHash{$target})) {
        push @{$interactions{$source}}, $target;
      }
    }
    if (exists $genes{$target}) {
      my %currHash = map { $_ => 1 } @{$interactions{$target}};
      if (!(exists $currHash{$source})) {
        push @{$interactions{$target}}, $source;
      }
    }
  }
  close($fh);

  return \%interactions;
}
 