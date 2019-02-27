package Go::Gaf;

use strict;
use warnings;

use Exporter qw(import);
use Text::CSV_XS;
 
our @EXPORT_OK = qw(readGaf);

# Reads GO annotation file.
sub readGaf {
  my $file = $_[0];
  my $namespace = $_[1];

  print STDERR "Reading known localizations\n";
  my $tsv = Text::CSV_XS->new({
    binary => 1,
    sep_char => "\t",
    quote_char => undef,
    escape_char => undef,
    allow_loose_quotes => 1,
    allow_loose_escapes => 1,
  });
  open my $fh, '<', $file or die "Could not open $file: $!";
  my %geneAnnotations;
  while(my $row = $tsv->getline($fh)) {
    my $ns = @{$row}[8];
    if ($ns eq $namespace) {
      my $gene = @{$row}[2];
      my $term = @{$row}[4];
      push @{$geneAnnotations{$gene}}, $term;
    }
  }
  close $fh;

  return \%geneAnnotations;
}
 