#!/usr/bin/perl

# 22/11/2018

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # NMF basis file

if ($#ARGV==0){
	print "\nTakes an NMF basis.csv file and determines prey moonlighting. A moonlighting\n";
  print "prey is one with a secondary localization within X% of its primary. X is set\n";
  print "beteen 0-1 and incremented by 0.01.\n";
	print "\nusage:\n $0\n";
	print "-b [NMF basis file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		}
		$i++;
	}
}

# Read basis matrix
print STDERR "Reading basis matrix\n";

open my $basisFH, '<', $bfile or die "Could not open $bfile: $!";
my $basisTsv = Text::CSV_XS->new({ sep_char => "," });

## Discard header
$basisTsv -> getline($basisFH);

## Read NMF profile for each prey
my %preys;
while(my $row = $basisTsv->getline($basisFH)) {
  my $gene = @{$row}[0];
  shift @{$row};
  @{$preys{$gene}} = @{$row};
}
close $basisFH;

# Open file for output
open my $outfile, '>', 'moonlighting.tsv';
print $outfile "percent threshold\tfraction\n";

# Calculate %preys moonlighting.
my $totalPreys = scalar(keys %preys);
for (my $i = 100; $i >= 0; $i--) {
  my $moonlighting = 0;
  my $threshold = $i / 100;
  foreach my $prey (keys %preys) {
    my $primary = 0;
    my $secondary = 0;
    foreach my $value (@{$preys{$prey}}) {
      if ($value > $primary) {
        $secondary = $primary;
        $primary = $value;
      }
    }
    my $diff = ($primary - $secondary) / $primary;
    if ($diff <= $threshold) {
      $moonlighting++;
    }
  }
  my $frac = sprintf("%.2f", $moonlighting / $totalPreys);
  print $outfile "$i\t$frac\n";
}

close $outfile;
