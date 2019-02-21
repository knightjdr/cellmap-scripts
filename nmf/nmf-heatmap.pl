#!/usr/bin/perl

# 7/7/2017

use strict;
use warnings;

# modules
use Math::Round;
use Text::CSV_XS;

# paremeters
my $similarity = 0.05; # how close do two values need to be to be an edge
my @colors = (
	"#FFFFFF", "#1CE6FF", "#FF4A46", "#008941", "#006FA6", "#8FB0FF", "#FFB500",
	"#809693", "#4FC601", "#3B5DFF", "#00C2A0", "#FFAA92", "#D16100", "#7B4F4B",
	"#A1C299", "#0AA6D8", "#00846F", "#997D87", "#A079BF", "#C0B9B2", "#C2FF99",
	"#00489C", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
	"#885578", "#FAD09F", "#FF8A9A", "#BEC459", "#456648", "#0086ED", "#886F4C",
	"#B4A8BD", "#00A6AA", "#636375", "#A3C8C9", "#FF913F", "#938A81", "#00FECF",
	"#B05B6F", "#8CD0FF", "#C8A1A1", "#A77500", "#6367A9", "#A05837", "#D790FF",
	"#9B9700", "#549E79", "#72418F", "#99ADC0", "#0089A3", "#CB7E98", "#324E72",
	"#83AB58", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#BF5650", "#66796D",
	"#8ADBB4", "#C895C5", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379"
);
my %transparency = (
	"100" => "FF",
	"95" => "F2",
	"90" => "E6",
	"85" => "D9",
	"80" => "CC",
	"75" => "BF",
	"70" => "B3",
	"65" => "A6",
	"60" => "99",
	"55" => "8C",
	"50" => "80",
	"45" => "73",
	"40" => "66",
	"35" => "59",
	"30" => "4D",
	"25" => "40",
	"20" => "33",
	"15" => "26",
	"10" => "1A",
	"5" => "0D",
	"0" => 00,
);

# command line parameters
my $bfile = '';	# basis.csv
my $rfile = '';	# rank summary file
my $sfile = ''; # scores.csv

if ($#ARGV == 0) {
	print "\nTakes NMF output files and converts them to ProHits-viz compatible format.\n";
	print "Also produces a cytoscape compatible file.\n";
	print "\nusage:\n $0\n";
	print "-b [basis.csv]\n";
  print "-r [rank summary.txt]\n";
	print "-s [scores.csv]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-r'){
      $i++;
      $rfile = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Read rank names.
my @rankNames;
push @rankNames, '';
open my $rankFH, '<', $rfile or die "Could not open $rfile: $!";
my $rankTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$rankTSV->getline($rankFH); # discard header
while(my $row = $rankTSV->getline($rankFH)) {
	my ($rank) = @{$row}[1] =~ /\[(.+)\]/;
	push  @rankNames, $rank;
}
close $rankFH;

my $formatDetails = '{"type": "heatmap", "kind": "rank vs prey", "xAxis": "Rank", "yAxis": "Prey", "filterType": 0, "primary": 0.01, "secondary": 0.05, "score": "N/A", "abundance": "N/A"}';

# basis matrix
open my $basisFH, '>', 'pv_basis.tsv';
print $basisFH "row\tcolumn\tvalue\tparams\n";
open my $rankOutputFH, '>', 'rank-color.tsv';
print $rankOutputFH "name\trank\trank name\tcolor\ttransparency\n";
open my $basisInputFH, '<', $bfile or die "Could not open $bfile: $!";
my $basisTSV = Text::CSV_XS->new({});
my $header = $basisTSV->getline($basisInputFH); # discard header
my $noRanks = scalar @{$header} - 1;
my $formatPrinted = 0;
my %genes;
while(my $row = $basisTSV->getline($basisInputFH)) {
	my $bestRank = 0;
	my $bestValue = 0;
  for(my $i = 1; $i <= $noRanks; $i++) {
    my $value = sprintf "%.4f", @{$row}[$i];
    print $basisFH "@{$row}[0]\t@{$header}[$i]\t$value";
    if (!$formatPrinted) {
      $formatPrinted = 1;
      print $basisFH "\t$formatDetails\n"
    } else {
      print $basisFH "\t \n";
    }
		$genes{@{$row}[0]}[$i] = $value;
		if ($value > $bestValue) {
			$bestValue = $value;
			$bestRank = @{$header}[$i];
		}
  }
	# calculate transparency based on NMF value of best rank
	# my $transparencyValue = nearest(5, $bestValue * 100);
	# create color based on transparency and color for best rank
	#my $color = '#' . $transparency{$transparencyValue} . substr $colors[$bestRank], 1;
	my $currTransparency = round $bestValue * 255;
	print $rankOutputFH "@{$row}[0]\t$bestRank\t$rankNames[$bestRank]\t$colors[$bestRank]\t$currTransparency\n";
}
close $rankOutputFH;
close $basisFH;
close $basisInputFH;

# scores matrix
open my $scoresFH, '>', 'pv_scores.tsv';
print $scoresFH "row\tcolumn\tvalue\tparams\n";
open my $scoresInputFH, '<', $sfile or die "Could not open $sfile: $!";
my $scoresTSV = Text::CSV_XS->new({});
$header = $scoresTSV->getline($scoresInputFH); # discard header
$noRanks = scalar @{$header} - 1;
$formatPrinted = 0;
while(my $row = $scoresTSV->getline($scoresInputFH)) {
  for(my $i = 1; $i <= $noRanks; $i++) {
    my $value = sprintf "%.4f", @{$row}[$i];
    print $scoresFH "@{$row}[0]\t@{$header}[$i]\t$value";
    if (!$formatPrinted) {
      $formatPrinted = 1;
      print $scoresFH "\t$formatDetails\n"
    } else {
      print $scoresFH "\t \n";
    }
  }
}
close $scoresFH;
close $scoresInputFH;

# print cytoscape files
my %nodes;
print STDERR "Generating cytoscape file\n";
foreach my $source (keys %genes) {
	foreach my $target (keys %genes) {
		if ($source ne $target) {
			my $bestRank = 0;
			my $bestValue = 0;
			for(my $i = 1; $i <= $noRanks; $i++) {
				if (
					$genes{$source}[$i] > 0 &&
					$genes{$target}[$i] > 0
				) {
					my $diff = abs($genes{$target}[$i] - $genes{$source}[$i]) / $genes{$source}[$i];
					if (
						$diff <= $similarity &&
						$genes{$source}[$i] >= $bestValue
					) {
						$bestValue = $genes{$source}[$i];
						$bestRank = $i;
					}
				}
			}
			if ($bestValue > 0) {
				my @sorted = sort ($source, $target);
				if (
					exists $nodes{$sorted[0]} &&
					exists $nodes{$sorted[0]}{$sorted[1]} &&
					$bestValue > $nodes{$sorted[0]}{$sorted[1]}{'value'}
				) {
					$nodes{$sorted[0]}{$sorted[1]}{'color'} = $colors[$bestRank];
					$nodes{$sorted[0]}{$sorted[1]}{'value'} = $bestValue;
				} else {
					$nodes{$sorted[0]}{$sorted[1]}{'color'} = $colors[$bestRank];
					$nodes{$sorted[0]}{$sorted[1]}{'value'} = $bestValue;
				}
			}
		}
	}
}
open my $cytoscapeFH, '>', 'cytoscape_nmf.tsv';
print $cytoscapeFH "source\ttarget\tedge\tedge color\n";
foreach my $source (keys %nodes) {
	foreach my $target (keys %{$nodes{$source}}) {
		print $cytoscapeFH "$source\t$target\t$nodes{$source}{$target}{'value'}\t$nodes{$source}{$target}{'color'}\n"
	}
}
close $cytoscapeFH;
