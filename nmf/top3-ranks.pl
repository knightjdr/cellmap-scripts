#!/usr/bin/perl

# 27/2/2019

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # NMF basis file
my $rfile = ''; # NMF summary file

if ($#ARGV==0){
	print "\nTakes an NMF basis.csv file and returns the top three ranks for each prey.\n";
	print "\nusage:\n $0\n";
	print "-b [NMF basis file]\n";
  print "-r [rank summary file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-r'){
      $i++;
      $rfile = $ARGV[$i];
    }
		$i++;
	}
}

# Read basis matrix
print STDERR "Reading basis matrix\n";

## Open file and discard header
open my $basisFH, '<', $bfile or die "Could not open $bfile: $!";
my $basisTsv = Text::CSV_XS->new({ sep_char => "," });
$basisTsv -> getline($basisFH);

## Read NMF profile for each prey
my %preys;
while(my $row = $basisTsv->getline($basisFH)) {
  my $gene = @{$row}[0];
  shift @{$row};
  @{$preys{$gene}} = @{$row};
}
close $basisFH;

# Read rank details
print STDERR "Reading rank details\n";

## open file and discard header
my $rankTSV = Text::CSV_XS->new({
  binary => 1,
  sep_char => "\t",
  escape_char => undef,
  allow_loose_escapes => 1,
});
open my $rankfh, '<', $rfile or die "Could not open $rfile: $!";
$rankTSV->getline($rankfh); # discard header

## Read details for each rank
my %rankDetails;
while(my $row = $rankTSV->getline($rankfh)) {
  my $rank = @{$row}[0];
  my $termString = @{$row}[1];
  $termString =~ s/[\[\]"]//g;
  $rankDetails{$rank} = $termString;
}
close $rankfh;

# Open file for output
open my $outfile, '>', 'top3-ranks.txt';
print $outfile "prey\tprimary\tsecondary\ttertiary\tprimary score\tsecondary score\ttertiary score\tsecondary / primary\ttertiary / primary\n";

# Iterate over preys and return top ranks
foreach my $gene (sort {lc $a cmp lc $b} keys %preys) {
  my %primary;
  $primary{'rank'} = 0;
  $primary{'score'} = 0;
  my %secondary;
  $secondary{'rank'} = 0;
  $secondary{'score'} = 0;
  my %tertiary;
  $tertiary{'rank'} = 0;
  $tertiary{'score'} = 0;

  while (my ($index, $score) = each @{$preys{$gene}}) {
    if ($score > $primary{'score'}) {
      $tertiary{'rank'} = $secondary{'rank'};
      $tertiary{'score'} = $secondary{'score'};
      $secondary{'rank'} = $primary{'rank'};
      $secondary{'score'} = $primary{'score'};
      $primary{'rank'} = $index + 1;
      $primary{'score'} = $score;
    } elsif ($score > $secondary{'score'}) {
      $tertiary{'rank'} = $secondary{'rank'};
      $tertiary{'score'} = $secondary{'score'};
      $secondary{'rank'} = $index + 1;
      $secondary{'score'} = $score;
    } elsif ($score > $tertiary{'score'}) {
      $tertiary{'rank'} = $index + 1;
      $tertiary{'score'} = $score;
    }
  }
  
  print $outfile "$gene\t$primary{'rank'} - $rankDetails{$primary{'rank'}}\t";
  my $secondaryName = '';
  if ($secondary{'rank'} > 0) {
    $secondaryName = " - $rankDetails{$secondary{'rank'}}";
  }
  my $tertiaryName = '';
  if ($tertiary{'rank'} > 0) {
    $tertiaryName = " - $rankDetails{$tertiary{'rank'}}";
  }
  print $outfile "$secondary{'rank'}$secondaryName\t";
  print $outfile "$tertiary{'rank'}$tertiaryName\t";

  my $primaryScore = sprintf("%.3f", $primary{'score'});
  my $secondaryScore = sprintf("%.3f", $secondary{'score'});
  my $tertiaryScore = sprintf("%.3f", $tertiary{'score'});
  print $outfile "$primaryScore\t$secondaryScore\t$tertiaryScore\t";

  my $secondaryRatio = sprintf("%.3f", $secondary{'score'} / $primary{'score'});
  my $tertiaryRatio = sprintf("%.3f", $tertiary{'score'} / $primary{'score'});
  print $outfile "$secondaryRatio\t$tertiaryRatio\n"
}

close $outfile;
