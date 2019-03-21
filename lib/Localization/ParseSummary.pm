package Localization::ParseSummary;

use strict;
use warnings;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(parseSummary);

my %tsvParams = (
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
);

# Parse localizations from NMF or SAFE summary file
sub parseSummary {
  my $file = $_[0];

  my $noComparments = 0;
  my %compartmentToGo;

  open my $fh, '<', $file or die "Could not open $file: $!";
  my $tsv = Text::CSV_XS->new(\%tsvParams);
  $tsv->getline($fh); # discard header
  while(my $row = $tsv->getline($fh)) {
    $noComparments++;
    my $idString = @{$row}[3];
    $idString =~ s/[\[\]"]//g;
    my $termString = @{$row}[1];
    $termString =~ s/[\[\]"]//g;
    @{$compartmentToGo{'term'}[@{$row}[0]]} = split /, /, $termString;
    @{$compartmentToGo{'go'}[@{$row}[0]]} = split /, /, $idString;
  }
  close $fh;

  return($noComparments, \%compartmentToGo);
}
