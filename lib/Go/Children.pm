package Go::Children;

use strict;
use warnings;

use Exporter qw(import);
use Text::CSV_XS;
 
our @EXPORT_OK = qw(readChildren);

# Reads GO children from a txt file.
sub readChildren {
  my $file = $_[0];

  print STDERR "Reading hierarchy file\n";
  open my $fh, '<', $file or die "Could not open $file: $!";
  my $tsv = Text::CSV_XS->new({
    binary => 1,
    sep_char => "\t",
    quote_char => undef,
    escape_char => undef,
    allow_loose_quotes => 1,
    allow_loose_escapes => 1,
  });
  my %children;
  while(my $row = $tsv->getline($fh)) {
    @{$children{@{$row}[0]}} = split ',', @{$row}[1];
  }
  close $fh;

  return \%children;
}
 