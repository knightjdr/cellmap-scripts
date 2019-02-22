#!/usr/bin/perl

# 8/6/2017

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV;

# parameters
my $fdr = 0.01;
my @mitoCytosol = ('AKAP1', 'AKAP1_sigseq');
my @mitoInter = ('AIFM1', 'COX4I1', 'IMMT');
my @mitoMatrix = ('COX8A_target', 'CS', 'PDHA1');

# command line parameters
my $efile = ''; # list of genes to highlight
my $sfile = ''; # SAINT file

if($#ARGV==0){
	print "Takes a SAINT file, and a list of mito proteins, and calculates the overlap\n";
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
		(	grep(/^$currBait$/, @mitoCytosol) ||
			grep (/^$currBait$/, @mitoInter) ||
			grep (/^$currBait$/, @mitoMatrix)
	)) {
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
	my $cytosolBoolean = 0;
	my $cytosolCount = 0;
	my $interBoolean = 0;
	my $interCount = 0;
	my $matrixBoolean = 0;
	my $matrixCount = 0;
	for(my $j = 0; $j < scalar @matchedBaits; $j++) {
		if (grep(/^$matchedBaits[$j]$/, @mitoCytosol)) {
			$cytosolBoolean = 1;
			$cytosolCount++;
		} elsif (grep(/^$matchedBaits[$j]$/, @mitoInter)) {
			$interBoolean = 1;
			$interCount++;
		} elsif (grep(/^$matchedBaits[$j]$/, @mitoMatrix)) {
			$matrixBoolean = 1;
			$matrixCount++;
		}
	}
	if (($cytosolBoolean && $interBoolean && $matrixBoolean) ||
		($cytosolBoolean && $matrixBoolean)
	) {
		my $matched = $cytosolCount + $interCount + $matrixCount;
		push @{$localized{'inter'}{$matched}}, $preyList[$i];
	} elsif($cytosolBoolean && $interBoolean) {
		my $matched = $cytosolCount + $interCount;
		push @{$localized{'outermembrane'}{$matched}}, $preyList[$i];
	} elsif($interBoolean && $matrixBoolean) {
		my $matched = $interCount + $matrixCount;
		push @{$localized{'innermembrane'}{$matched}}, $preyList[$i];
	} elsif($cytosolBoolean) {
		push @{$localized{'cytosol'}{$cytosolCount}}, $preyList[$i];
	} elsif($interBoolean) {
		push @{$localized{'inter'}{$interCount}}, $preyList[$i];
	} elsif($matrixBoolean) {
		push @{$localized{'matrix'}{$matrixCount}}, $preyList[$i];
	}
}

# Output
open my $detailsFH, '>', 'mito_details.txt';
my $totalHighlight = scalar @highlightGenes;
print $detailsFH "Total query\t$totalHighlight\n\n";
print $detailsFH "Cytosol\n";
for my $matches (keys %{$localized{'cytosol'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'cytosol'}{$matches}}; $i++) {
			my $currGene = $localized{'cytosol'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'cytosol'}{$matches}[$i] = $localized{'cytosol'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'cytosol'}{$matches}} = sort @{$localized{'cytosol'}{$matches}};
	my $joinedList = join ';', @{$localized{'cytosol'}{$matches}};
	my $length = scalar @{$localized{'cytosol'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nOutermembrane\n";
for my $matches (keys %{$localized{'outermembrane'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'outermembrane'}{$matches}}; $i++) {
			my $currGene = $localized{'outermembrane'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'outermembrane'}{$matches}[$i] = $localized{'outermembrane'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'outermembrane'}{$matches}} = sort @{$localized{'outermembrane'}{$matches}};
	my $joinedList = join ';', @{$localized{'outermembrane'}{$matches}};
	my $length = scalar @{$localized{'outermembrane'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nIntermembrane space\n";
for my $matches (keys %{$localized{'inter'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'inter'}{$matches}}; $i++) {
			my $currGene = $localized{'inter'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'inter'}{$matches}[$i] = $localized{'inter'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'inter'}{$matches}} = sort @{$localized{'inter'}{$matches}};
	my $joinedList = join ';', @{$localized{'inter'}{$matches}};
	my $length = scalar @{$localized{'inter'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nInnermembrane\n";
for my $matches (keys %{$localized{'innermembrane'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'innermembrane'}{$matches}}; $i++) {
			my $currGene = $localized{'innermembrane'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'innermembrane'}{$matches}[$i] = $localized{'innermembrane'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'innermembrane'}{$matches}} = sort @{$localized{'innermembrane'}{$matches}};
	my $joinedList = join ';', @{$localized{'innermembrane'}{$matches}};
	my $length = scalar @{$localized{'innermembrane'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
print $detailsFH "\n\nMatrix\n";
for my $matches (keys %{$localized{'matrix'}}) {
	my $highlighted = 0;
	if (scalar @highlightGenes > 0) {
		for(my $i = 0 ; $i < scalar @{$localized{'matrix'}{$matches}}; $i++) {
			my $currGene = $localized{'matrix'}{$matches}[$i];
			if (grep(/^$currGene$/, @highlightGenes)) {
				$highlighted++;
				$localized{'matrix'}{$matches}[$i] = $localized{'matrix'}{$matches}[$i] . '*';
			}
		}
	}
	@{$localized{'matrix'}{$matches}} = sort @{$localized{'matrix'}{$matches}};
	my $joinedList = join ';', @{$localized{'matrix'}{$matches}};
	my $length = scalar @{$localized{'matrix'}{$matches}};
	print $detailsFH "$matches\t$length\t$highlighted\t$joinedList\n";
}
close $detailsFH;
