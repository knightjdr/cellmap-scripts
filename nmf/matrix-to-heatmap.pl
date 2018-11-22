#!/usr/bin/perl

# 22/11/2018

use strict;
use warnings;

# modules
use Math::Round;
use Text::CSV_XS;

# paremeters

# command line parameters
my $mfile = '';	# matrix file in csv format

if ($#ARGV == 0) {
	print "\nTakes a matrix and converts it to ProHits-viz compatible format.\n";
	print "\nusage:\n $0\n";
	print "-m [matrix in csv format]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-m'){
			$i++;
			$mfile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

my $formatDetails = '{"type": "heatmap", "kind": "rank vs prey", "xAxis": "Rank", "yAxis": "Prey", "filterType": 0, "primary": 0.01, "secondary": 0.05, "score": "N/A", "abundance": "N/A"}';

# basis matrix
open my $matrixFH, '>', 'matrix.tsv';
print $matrixFH "row\tcolumn\tvalue\tparams\n";
open my $matrixInputFH, '<', $mfile or die "Could not open $mfile: $!";
my $matrixTSV= Text::CSV_XS->new({});
my $header = $matrixTSV->getline($matrixInputFH); # discard header
my $noRanks = scalar @{$header} - 1;
my $formatPrinted = 0;
while(my $row = $matrixTSV->getline($matrixInputFH)) {
  for(my $i = 1; $i <= $noRanks; $i++) {
    my $value = sprintf "%.4f", @{$row}[$i];
    print $matrixFH "@{$row}[0]\t@{$header}[$i]\t$value";
    if (!$formatPrinted) {
      $formatPrinted = 1;
      print $matrixFH "\t$formatDetails\n"
    } else {
      print $matrixFH "\t \n";
    }
  }
}
close $matrixFH;
close $matrixInputFH;
