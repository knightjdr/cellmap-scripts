#!/usr/bin/perl

# 10/7/2017

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# parameters
my $fdr = 0.01;

# command line parameters
my $mfile = ''; # bait gene map file
my $sfile = ''; # SAINT file

if ($#ARGV==0) {
	print "\nTakes in a SAINT file and calculates summary statistics at a 0.01 FDR cutoff.\n";
	print "\nusage:\n $0\n";
	print "-f [FDR cutoff]\n";
	print "-m [bait gene map]\n";
  print "-s [SAINT file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-f'){
			$i++;
			$fdr = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-m'){
			$i++;
			$mfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read bait map
print STDERR "Creating map\n";
my %baitMap;
my $mapTSV = Text::CSV_XS->new({ sep_char => "\t" });
open my $mapFH, '<', $mfile or die "Could not open $mfile: $!";
$mapTSV->getline($mapFH); #discard header
while(my $row = $mapTSV->getline($mapFH)) {
	my $bait = lc @{$row}[0];
	my $gene = lc @{$row}[1];
	$baitMap{$bait} = $gene;
}
close($mapFH);

# read SAINT file
print STDERR "Reading saint file\n";
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
my $saintTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
my @baits;
my @interactionPairs;
my @preys;
my $totalInteractions = 0;
$saintTSV->getline($saintFH); # discard header
while(my $row = $saintTSV->getline($saintFH)) {
	if (@{$row}[15] <= $fdr) {
		push @baits, @{$row}[0];
		push @preys, @{$row}[2];
		$totalInteractions++;
		my @pair = (lc @{$row}[0], lc @{$row}[2]);
		my $sortedPair = join '_', sort @pair;
		push @interactionPairs, $sortedPair;
	}
}
close $saintFH;

# remove duplicates
@baits = uniq @baits;
@interactionPairs = uniq @interactionPairs;
@preys = uniq @preys;

# For each bait, get it's gene name so that I can count unique genes.
my @baitgenes;
foreach my $bait (@baits) {
	push @baitgenes, $baitMap{lc $bait};
}
@baitgenes = uniq @baitgenes;

# count
my $noBaits = scalar @baits;
my $noGenes = scalar @baitgenes;
my $uniqueInteractions = scalar @interactionPairs;
my $noPreys = scalar @preys;

# output
open my $summaryFH, '>', 'summary.txt';
print $summaryFH "baits\t$noBaits\n";
print $summaryFH "unique bait genes\t$noGenes\n";
print $summaryFH "total interactions\t$totalInteractions\n";
print $summaryFH "unique interactions\t$uniqueInteractions\n";
print $summaryFH "unique preys\t$noPreys\n";
close $summaryFH;
