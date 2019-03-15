#!/usr/bin/perl

# 7/6/2018

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use List::Util qw(min);
use POSIX qw(floor);
use Statistics::Basic qw(mean);
use String::Util qw(trim);
use Text::CSV_XS;

# parameters
my $baitBinSize = 5;
my $limitBaits = 15;
my $baitBins = 4; # 1-5, 6-10, 11-15, 15+

my $specBinSize = 5;
my $specBins = 9; # 0 - 9, 10 - 19, 20 - 29, 30 - 39, 40+

my $fdr = 0.01;
my $fileType = 'n';
my $namespace = 'C';

my %tsvParams = (
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
);

# command line parameters
my $hfile = ''; # hierarchy file
my $gfile = ''; # tsv file (with header), two columns with gene and it's rank/domain/GO term
my $lfile = ''; # GO localization file (goa_human.gaf)
my $rfile = ''; # for NMF and SAFE, file with rank/domain details
my $sfile = ''; # SAINT file.

if ($#ARGV==0){
	print "\nTakes GO hierarchy from get_children.pl, a list of official localizations for genes,\n";
	print "assigned localization from NMF or SAFE, a SAINT file, assesses the assigned localizations, and\n";
	print "also returns assessments on prediction quality relative to bait number and spectral count.\n";
	print "\nusage:\n $0\n";
	print "-g [gene localization file]\n";
  print "-h [hierarcy file procduced by get_children.pl]\n";
	print "-l [GO localization file (goa_human.gaf)]\n";
	print "-r [file with rank/domain details]\n";
	print "-s [SAINT file]\n";
	print "-t [type of assessment analyzed. n = NMF (default), h = HPA, o = other, s = SAFE]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-h'){
			$i++;
			$hfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-r'){
			$i++;
			$rfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
			$i++;
			$fileType = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# get child terms for each GO term
print STDERR "Reading hierarchy file\n";
open my $hierarchyFH, '<', $hfile or die "Could not open $hfile: $!";
my $tsv = Text::CSV_XS->new(\%tsvParams);
my %children;
while(my $row = $tsv->getline($hierarchyFH)) {
	@{$children{@{$row}[0]}} = split ',', @{$row}[1];
}
close $hierarchyFH;

# reading known localizations
print STDERR "Reading known localizations\n";
$tsv = Text::CSV_XS->new(\%tsvParams);
open my $localizationfh, '<', $lfile or die "Could not open $lfile: $!";
my %geneLocalizations;
while(my $row = $tsv->getline($localizationfh)) {
	if (@{$row}[8] eq $namespace) {
		push @{$geneLocalizations{@{$row}[2]}}, @{$row}[4];
	}
}
close $localizationfh;

# get rank details
my $noRanks = 0;
my %rankToGo;
if ($fileType eq 'n' || $fileType eq 's') {
	print STDERR "Reading rank details\n";
	my $tsv = Text::CSV_XS->new({
		binary => 1,
		sep_char => "\t",
		escape_char => undef,
		allow_loose_escapes => 1,
	});
	open my $rankfh, '<', $rfile or die "Could not open $rfile: $!";
	$tsv->getline($rankfh); # discard header
	while(my $row = $tsv->getline($rankfh)) {
		$noRanks++;
		my $idString = @{$row}[3];
		$idString =~ s/[\[\]"]//g;
		my $listedString = @{$row}[2];
		$listedString =~ s/[\[\]"]//g;
		my $termString = @{$row}[1];
		$termString =~ s/[\[\]"]//g;
		@{$rankToGo{'term'}[@{$row}[0]]} = split /, /, $termString;
		@{$rankToGo{'displayname'}[@{$row}[0]]} = split /, /, $listedString;
		@{$rankToGo{'go'}[@{$row}[0]]} = split /, /, $idString;
	}
	close $rankfh;
}

# get gene localizations
print STDERR "Reading assigned gene localizations\n";
my %actualLocalizationName;
my %assignedLocalization;
my %indexMap;
my @localizationsPerGene;
open my $assignedlocalizationfh, '<', $gfile or die "Could not open $gfile: $!";
if ($fileType eq 'n') {
	$tsv = Text::CSV_XS->new(\%tsvParams);
	$tsv->getline($assignedlocalizationfh); # discard header
	while(my $row = $tsv->getline($assignedlocalizationfh)) {
		$assignedLocalization{@{$row}[0]} = @{$row}[1];
	}
} elsif ($fileType eq 's') {
	$tsv = Text::CSV_XS->new(\%tsvParams);
	$tsv->getline($assignedlocalizationfh); # discard header
	while(my $row = $tsv->getline($assignedlocalizationfh)) {
		$assignedLocalization{@{$row}[1]} = @{$row}[2];
	}
}
close $assignedlocalizationfh;

# Read SAINT file and for each significant prey, count number of baits store avgspec.
my %preyInfo;
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
$tsv = Text::CSV_XS->new(\%tsvParams);
$tsv->getline($saintFH); # discard header
while(my $row = $tsv->getline($saintFH)) {
	my $avgspec = @{$row}[5];
	my $prey = @{$row}[2];
	if (@{$row}[15] <= $fdr) {
		if (exists $preyInfo{'baits'}{$prey}) {
			$preyInfo{'baits'}{$prey}++;
			push @{$preyInfo{'avgspec'}{$prey}}, $avgspec;
		} else {
			$preyInfo{'baits'}{$prey} = 1;
			@{$preyInfo{'avgspec'}{$prey}} = [$avgspec];
		}
	}
}
close $saintFH;

# Calculate average spectral counts for each prey across the baits it was found with.
# Average of the AvgSpec. Also find the maximum bait number.
my $maxBait = 0;
foreach my $prey (keys %{$preyInfo{'baits'}}) {
	$preyInfo{'avg'}{$prey} = mean @{$preyInfo{'avgspec'}{$prey}};
	if ($preyInfo{'baits'}{$prey} > $maxBait) {
		$maxBait = $preyInfo{'baits'}{$prey};
	}
}

# Init arrays for storing the number of preys per bait number and the
# number of correct prey localization per bait number.
my %preyBaitTotals;
for (my $i = 0; $i <= $limitBaits; $i++) {
	$preyBaitTotals{$i}{'total'} = 0;
	$preyBaitTotals{$i}{'known'} = 0;
}

# Init arrays for storing the number of preys per average spec number and the
# number of correct prey localization per avgerage spec number.
my %preyAvgTotals;
for (my $i = 0; $i < $specBins; $i++) {
	$preyAvgTotals{$i}{'total'} = 0;
	$preyAvgTotals{$i}{'known'} = 0;
}

# Init arrays for storing the number of preys per avgerage spec number and bait number.
my %preyCombinedTotals;
for (my $i = 0; $i < $baitBins * $specBins; $i++) {
	$preyCombinedTotals{$i}{'total'} = 0;
	$preyCombinedTotals{$i}{'known'} = 0;
}

# assess localization
print STDERR "Assessing assigned localizations\n";
open my $outputFH, '>', 'metrics.txt';
my $correctAssignments = 0;
my $totalGenes = 0;
if ($fileType eq 'n' || $fileType eq 's') {
	while( my( $gene, $rank ) = each %assignedLocalization ) {
		my $noBaits = min($preyInfo{'baits'}{$gene}, $limitBaits);
		my $avgBin = min(floor($preyInfo{'avg'}{$gene} / $specBinSize), $specBins - 1);
		my $combinedBin = ($avgBin * $baitBins) + min(floor(($preyInfo{'baits'}{$gene} - 1)/ 5), 3);
		$totalGenes++;
		my @possibleTerms = @{$rankToGo{'go'}[$rank]};
		my $noInitialTerms = scalar @possibleTerms;
		push @localizationsPerGene, $noInitialTerms;
		foreach my $term (@possibleTerms) {
			if (exists $children{$term}) {
				push @possibleTerms, @{$children{$term}};
			}
		}
		@possibleTerms = uniq @possibleTerms;
		if (scalar @possibleTerms > 0 && exists $geneLocalizations{$gene}) {
			foreach my $term (@{$geneLocalizations{$gene}}) {
				if (grep /^$term$/, @possibleTerms) {
					$correctAssignments++;
					$preyAvgTotals{$avgBin}{'known'}++;
					$preyCombinedTotals{$combinedBin}{'known'}++;
					$preyBaitTotals{$noBaits}{'known'}++;
					last;
				}
			}
		}
		$preyAvgTotals{$avgBin}{'total'}++;
		$preyCombinedTotals{$combinedBin}{'total'}++;
		$preyBaitTotals{$noBaits}{'total'}++;
	}
	my $totalFrac = sprintf "%.3f", $correctAssignments / $totalGenes;
	print $outputFH "total\t$totalFrac\t$totalGenes\t$correctAssignments\n";
	# Fraction of known localizations per number of baits a prey was seen with.
	print $outputFH "\nPreys per bait summary\n";
	print $outputFH "baits\tfraction\ttotal localizations\tknown localizations\n";
	for (my $i = 1; $i <= $limitBaits; $i++) {
		my $baitFrac = sprintf "%.3f", $preyBaitTotals{$i}{'known'} / $preyBaitTotals{$i}{'total'};
		print $outputFH "$i\t$baitFrac\t$preyBaitTotals{$i}{'total'}\t$preyBaitTotals{$i}{'known'}\n";
	}
	# Fraction of known localizations per avg spec a prey was seen with.
	print $outputFH "\nPreys per avg spec summary\n";
	print $outputFH "spec\tfraction\ttotal localizations\tknown localizations\n";
	for (my $i = 0; $i < $specBins; $i++) {
		my $avgFrac = sprintf "%.3f", $preyAvgTotals{$i}{'known'} / $preyAvgTotals{$i}{'total'};
		my $bottom = ($i * $specBinSize);
		my $top = $bottom + $specBinSize - 1;
		print $outputFH "$bottom - $top\t$avgFrac\t$preyAvgTotals{$i}{'total'}\t$preyAvgTotals{$i}{'known'}\n";
	}
	# Fraction of known localizations per baits and avg spec.
	print $outputFH "\nCombined spec summary\n";
	print $outputFH "spec\tbait\tfraction\ttotal localizations\tknown localizations\n";
	for (my $i = 0; $i < $baitBins * $specBins; $i++) {
		my $bottom = floor($i / 4) * $specBinSize;
		my $top = $bottom + $specBinSize - 1;
		my $avgFrac = sprintf "%.3f", $preyCombinedTotals{$i}{'known'} / $preyCombinedTotals{$i}{'total'};
		print $outputFH "$bottom - $top\t$i\t$avgFrac\t$preyCombinedTotals{$i}{'total'}\t$preyCombinedTotals{$i}{'known'}\n";
	}
}

close $outputFH;
