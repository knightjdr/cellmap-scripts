#!/usr/bin/perl

# 8/6/2017

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV;

# parameters
my $fdr = 0.01;
my @erCytosol = (
  'ATP2A1',
  'CKAP4',
  'CYP2C1_sigseq',
  'CYP2C9',
  'DERL1',
  'ELOVL5',
  'LNPK_Cterm',
  'LNPK_Nterm',
  'LRRC59_Nterm',
  'REEP5',
  'RPN1',
  'RPN2',
  'SEC61B_Nterm',
  'SEC62',
  'SSR1',
  'STIM1'
);
my %erCytosolHash = map { $_ => 1 } @erCytosol;
my @erLumen = (
  'BCAP31',
  'CALR3',
  'CALU',
  'LRRC59_Cterm',
  'PDIA4',
  'SEC61B_Cterm',
  'VAPA'
);
my %erLumenHash = map { $_ => 1 } @erLumen;

# command line parameters
my $efile = ''; # list of genes to highlight
my $sfile = ''; # SAINT file

if($#ARGV==0){
	print "Takes a SAINT file, and a list of ER proteins, and calculates the overlap\n";
	print "between the cytosol and lumen.\n\n";
	print "\nusage:\n $0\n";
	print "-e [genes to highlight]\n";
  print "-f [FDR cutoff]\n";
	print "-s [SAINT file]\n\n";
	die;
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-e'){
			$i++;
			$efile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-f'){
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
my %highlightGenes;
if ($efile) {
	print STDERR "Reading list of genes to highlight\n";
	open my $geneFH, '<', $efile or die "Could not open $efile: $!";
	my $tsv = Text::CSV->new({
		sep_char => "\t",
	});
	while(my $row = $tsv->getline($geneFH)) {
    $highlightGenes{@{$row}[0]} = 1;
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
	if (
    $currFDR <= $fdr
    && (
      grep(/^$currBait$/, @erCytosol)
      || grep (/^$currBait$/, @erLumen)
    )
	) {
		push @{$baitPreyList{$currBait}}, $currPrey;
		push @preyList, $currPrey;
	}
}
close $saintFH;
@preyList = uniq @preyList;

# assigning preys to localization
print STDERR "Assigning localization\n";
my %localized;
for(my $i = 0; $i < scalar @preyList; $i++) {
	my @matchedBaits;
	for my $bait (keys %baitPreyList) {
		if (grep(/^$preyList[$i]/, @{$baitPreyList{$bait}})) {
			push @matchedBaits, $bait;
		}
	}
	my $cytosolBoolean = 0;
	my $cytosolCount = 0;
	my $lumenBoolean = 0;
	my $lumenCount = 0;
	for(my $j = 0; $j < scalar @matchedBaits; $j++) {
		if (grep(/^$matchedBaits[$j]$/, @erCytosol)) {
			$cytosolBoolean = 1;
			$cytosolCount++;
		} elsif (grep(/^$matchedBaits[$j]$/, @erLumen)) {
			$lumenBoolean = 1;
			$lumenCount++;
		}
	}
	if ($cytosolBoolean && $lumenBoolean) {
		my $matched = $cytosolCount + $lumenCount;
		push @{$localized{'transmembrane'}{$matched}}, $preyList[$i];
	} elsif($cytosolBoolean) {
		push @{$localized{'cytosol'}{$cytosolCount}}, $preyList[$i];
	} elsif($lumenBoolean) {
		push @{$localized{'lumen'}{$lumenCount}}, $preyList[$i];
	}
}

# Output
open my $detailsFH, '>', 'er-details.txt';
my $totalHighlight = keys %highlightGenes;
my $totalMatched = 0;
print $detailsFH "Cytosol\n";
for my $matches (keys %{$localized{'cytosol'}}) {
	my $highlighted = 0;
	if ($totalHighlight > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'cytosol'}{$matches}}; $i++) {
			my $currGene = $localized{'cytosol'}{$matches}[$i];
			if (exists $highlightGenes{$currGene}) {
        delete $highlightGenes{$currGene}; 
				$highlighted++;
        $totalMatched++;
				$localized{'cytosol'}{$matches}[$i] = $localized{'cytosol'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'cytosol'}{$matches}} = sort @{$localized{'cytosol'}{$matches}};
	my $joinedList = join ';', @{$localized{'cytosol'}{$matches}};
	my $length = scalar @{$localized{'cytosol'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nTransmembrane\n";
for my $matches (keys %{$localized{'transmembrane'}}) {
	my $highlighted = 0;
	if ($totalHighlight > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'transmembrane'}{$matches}}; $i++) {
			my $currGene = $localized{'transmembrane'}{$matches}[$i];
			if (exists $highlightGenes{$currGene}) {
        delete $highlightGenes{$currGene}; 
				$highlighted++;
        $totalMatched++;
				$localized{'transmembrane'}{$matches}[$i] = $localized{'transmembrane'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'transmembrane'}{$matches}} = sort @{$localized{'transmembrane'}{$matches}};
	my $joinedList = join ';', @{$localized{'transmembrane'}{$matches}};
	my $length = scalar @{$localized{'transmembrane'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nLumen\n";
for my $matches (keys %{$localized{'lumen'}}) {
	my $highlighted = 0;
	if ($totalHighlight > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'lumen'}{$matches}}; $i++) {
			my $currGene = $localized{'lumen'}{$matches}[$i];
			if (exists $highlightGenes{$currGene}) {
        delete $highlightGenes{$currGene}; 
				$highlighted++;
        $totalMatched++;
				$localized{'lumen'}{$matches}[$i] = $localized{'lumen'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'lumen'}{$matches}} = sort @{$localized{'lumen'}{$matches}};
	my $joinedList = join ';', @{$localized{'lumen'}{$matches}};
	my $length = scalar @{$localized{'lumen'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\nTotal query\t$totalHighlight\n";
print $detailsFH "total matches\t$totalMatched\n";
my $missed = join ';', keys %highlightGenes;
print $detailsFH "missed\t$missed\n";
close $detailsFH;
