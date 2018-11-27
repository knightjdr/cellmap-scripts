#!/usr/bin/perl

# 22/11/2018

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # NMF basis file
my @poolA;
my @poolB;
my $threshold = 0.3;

if ($#ARGV==0){
	print "\nTakes an NMF basis.csv file, a moonlighting threshold and two pools of NMF ranks\n";
  print "and subsets the basis matrix to only include preys that are moonlighting in the supplied\n";
  print "pools.\n";
	print "\nusage:\n $0\n";
	print "-b [NMF basis file]\n";
  print "-poola [list of ranks to include in pool A; comma separated list]\n";
  print "-poolb [list of ranks to include in pool A; comma separated list]\n";
  print "-t [moonlighting threshold]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-poola'){
      $i++;
      @poolA = split(/,/, $ARGV[$i]);
    } elsif ($ARGV[$i] eq '-poolb'){
      $i++;
      @poolB = split(/,/, $ARGV[$i]);
    } elsif ($ARGV[$i] eq '-t'){
      $i++;
      $threshold = $ARGV[$i];
    } 
		$i++;
	}
}

if (scalar @poolA == 0) {
  die "Supply a list of ranks to use for pool A: $!";
} elsif (scalar @poolB == 0) {
  die "Supply a list of ranks to use for pool B: $!";
}

# Read basis matrix
print STDERR "Reading basis matrix\n";

open my $basisFH, '<', $bfile or die "Could not open $bfile: $!";
my $basisTsv = Text::CSV_XS->new({ sep_char => "," });

## Discard header
my $header = $basisTsv -> getline($basisFH);

## Read NMF profile for each prey
my %preys;
while(my $row = $basisTsv->getline($basisFH)) {
  my $gene = @{$row}[0];
  shift @{$row};
  @{$preys{$gene}} = @{$row};
}
close $basisFH;

# Open file for output
open my $outfile, '>', 'basis-moonlighting.csv';
my $headerString = join(',', @{$header});
print $outfile "$headerString\n";

# Calculate preys moonlighting
my $totalPreys = scalar(keys %preys);
$threshold = $threshold / 100;
foreach my $prey (keys %preys) {
  my $index = 0;
  my $primary = 0;
  my $primaryRank = 0;
  my $secondary = 0;
  my $secondaryRank = 0;
  foreach my $value (@{$preys{$prey}}) {
    if ($value > $primary) {
      $secondary = $primary;
      $secondaryRank = $primaryRank;
      $primary = $value;
      $primaryRank = $index;
    }
    $index++;
  }
  $primaryRank += 1;
  $secondaryRank += 1;
  my $diff = ($primary - $secondary) / $primary;
  if (
    $diff <= $threshold
    && (
      (grep( /^$primaryRank$/, @poolA) && grep( /^$secondaryRank$/, @poolB))
      || (grep( /^$secondaryRank$/, @poolA) && grep( /^$primaryRank$/, @poolB))
    )
  ) {
    my $rowString = join(',', @{$preys{$prey}});
    print $outfile "$prey,$rowString\n";
  }
}

close $outfile;
