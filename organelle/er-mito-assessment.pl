#!/usr/bin/perl

# 8/6/2017

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV;

# parameters
my $fdr = 0.01;
my @erCytosol = ('ATP2A1', 'DERL1', 'STIM1', 'ELOVL5'); # omit CYP2C1_sigseq, future add STIM1 and ELOVL5
my @mitoCytosol = ('AKAP1', 'PEX3');

# command line parameters
my $efile = ''; # list of genes to highlight
my $sfile = ''; # SAINT file

if($#ARGV==0){
	print "\nTakes a SAINT file, and a list of ER or mito proteins, and calculates the overlap\n";
	print "between the cytoplasmic faces of each.\n\n";
	print "\nusage:\n $0\n";
	print "-e [genes to highlight]\n";
  print "-f [FDR cutoff]\n";
	print "-s [SAINT file]\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if($ARGV[$i] eq '-e'){
			$i++;
			$efile = $ARGV[$i];
		} if($ARGV[$i] eq '-f'){
      $i++;
      $fdr = $ARGV[$i];
    } elsif($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read list of genes to highlight
my @highlightGenes;
if ($efile) {
	print STDERR "Reading list of genes to highlight\n";
	open my $geneFH, '<', $efile or die "Could not open $efile: $!";
	my $tsv = Text::CSV->new({
		sep_char => "\t",
	});
	while(my $row = $tsv->getline($geneFH)) {
		push @highlightGenes, @{$row}[0];
	}
	close $geneFH;
}

# read SAINT file
print STDERR "Reading SAINT file\n";
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
my $tsv = Text::CSV->new({
	sep_char => "\t",
});
$tsv->getline($saintFH); # discard header
my %baitPreyList;
my @preyList;
while(my $row = $tsv->getline($saintFH)) {
	my $currBait = @{$row}[0];
	my $currPrey = @{$row}[2];
	my $currFDR = @{$row}[15];
	if ($currFDR <= $fdr &&
		(grep(/^$currBait$/, @erCytosol) || grep (/^$currBait$/, @mitoCytosol))
	) {
		push @{$baitPreyList{$currBait}}, $currPrey;
		push @preyList, $currPrey;
	}
}
close $saintFH;
@preyList = uniq @preyList;

# assigning preys to localization
print STDERR "Assinging localization\n";
my %localized;
for(my $i = 0; $i < scalar @preyList; $i++) {
	my @matchedBaits;
	for my $bait (keys %baitPreyList) {
		if (grep(/^$preyList[$i]/, @{$baitPreyList{$bait}})) {
			push @matchedBaits, $bait;
		}
	}
	my $erBoolean = 0;
	my $erCount = 0;
	my $mitoBoolean = 0;
	my $mitoCount = 0;
	for(my $j = 0; $j < scalar @matchedBaits; $j++) {
		if (grep(/^$matchedBaits[$j]$/, @erCytosol)) {
			$erBoolean = 1;
			$erCount++;
		} elsif (grep(/^$matchedBaits[$j]$/, @mitoCytosol)) {
			$mitoBoolean = 1;
			$mitoCount++;
		}
	}
	if ($erBoolean && $mitoBoolean) {
		my $matched = $erCount + $mitoCount;
		push @{$localized{'shared'}{$matched}}, $preyList[$i];
	} elsif($erBoolean) {
		push @{$localized{'er'}{$erCount}}, $preyList[$i];
	} elsif($mitoBoolean) {
		push @{$localized{'mito'}{$mitoCount}}, $preyList[$i];
	}
}

# Output
open my $detailsFH, '>', 'er_mito_details.txt';
my $totalHighlight = scalar @highlightGenes;
print $detailsFH "Total query\t$totalHighlight\n\n";
print $detailsFH "ER\n";
for my $matches (keys %{$localized{'er'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'er'}{$matches}}; $i++) {
			my $currGene = $localized{'er'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'er'}{$matches}[$i] = $localized{'er'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'er'}{$matches}} = sort @{$localized{'er'}{$matches}};
	my $joinedList = join ';', @{$localized{'er'}{$matches}};
	my $length = scalar @{$localized{'er'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nShared\n";
for my $matches (keys %{$localized{'shared'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'shared'}{$matches}}; $i++) {
			my $currGene = $localized{'shared'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'shared'}{$matches}[$i] = $localized{'shared'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'shared'}{$matches}} = sort @{$localized{'shared'}{$matches}};
	my $joinedList = join ';', @{$localized{'shared'}{$matches}};
	my $length = scalar @{$localized{'shared'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nMitochondria\n";
for my $matches (keys %{$localized{'mito'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'mito'}{$matches}}; $i++) {
			my $currGene = $localized{'mito'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'mito'}{$matches}[$i] = $localized{'mito'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'mito'}{$matches}} = sort @{$localized{'mito'}{$matches}};
	my $joinedList = join ';', @{$localized{'mito'}{$matches}};
	my $length = scalar @{$localized{'mito'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
close $detailsFH;
