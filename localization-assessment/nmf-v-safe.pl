#!/usr/bin/perl

# 22/8/2017

use strict;
use warnings;

# libraries
use Array::Utils qw(:all);
use Text::CSV_XS;

# command line parameters
my $gfile = ''; # file with lost of children per GO term
my $nfile = ''; # nmf summary file
my $nmfAssignedFile = ''; # file with nmf rank assigned to each gene
my $sfile = ''; # safe summary file
my $safeAssignedFile = ''; # file with safe assigned to each gene

if ($#ARGV==0) {
	print "\nTakes NMF and SAFE summaries and summarizes their overlap, etc.\n";
	print "\nusage:\n $0\n";
	print "-g [File with list of children for every GO term (from get_children.pl)]\n";
  print "-n [nmf_summary.txt]\n";
  print "-na [file with nmf rank assigned to each gene (gene-localizations.txt)]\n";
	print "-s [safe_summary.txt]\n";
  print "-sa [file with safe domain assigned to each gene (node_properties_annotation-highest.txt )]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV) {
		if ($ARGV[$i] eq '-g') {
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-n') {
			$i++;
			$nfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-na') {
			$i++;
			$nmfAssignedFile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s') {
			$i++;
			$sfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-sa') {
			$i++;
			$safeAssignedFile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# get child terms for each GO term
print STDOUT "Reading hierarchy file\n";
open my $gochildrenFH, '<', $gfile or die "Could not open $gfile: $!";
my $gochildrenTSV = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
my %children;
while(my $row = $gochildrenTSV->getline($gochildrenFH)) {
	@{$children{@{$row}[0]}} = split ',', @{$row}[1];
}
close $gochildrenFH;

# read NMF summary terms
my $noRanks = 0;
my %rankToGo;
print STDOUT "Reading NMF details\n";
open my $nmfFH, '<', $nfile or die "Could not open $nfile: $!";
my $nmfTSV = Text::CSV_XS->new({
  binary => 1,
  sep_char => "\t",
  escape_char => undef,
  allow_loose_escapes => 1,
});
$nmfTSV->getline($nmfFH); # discard header
while(my $row = $nmfTSV->getline($nmfFH)) {
  $noRanks++;
  my $idString = @{$row}[3];
  $idString =~ s/[\[\]"]//g;
  my $termString = @{$row}[1];
  $termString =~ s/[\[\]"]//g;
  @{$rankToGo{'term'}[@{$row}[0]]} = split /, /, $termString;
	@{$rankToGo{'go'}[@{$row}[0]]} = split /, /, $idString;
}
close $nmfFH;

# read SAFE summary terms
my $noDomains = 0;
my %domainToGo;
print STDOUT "Reading SAFE details\n";
open my $safeFH, '<', $sfile or die "Could not open $sfile: $!";
my $safeTSV = Text::CSV_XS->new({
  binary => 1,
  sep_char => "\t",
  escape_char => undef,
  allow_loose_escapes => 1,
});
$safeTSV->getline($safeFH); # discard header
while(my $row = $safeTSV->getline($safeFH)) {
  $noDomains++;
  my $idString = @{$row}[3];
  $idString =~ s/[\[\]"]//g;
  my $termString = @{$row}[1];
  $termString =~ s/[\[\]"]//g;
  @{$domainToGo{'term'}[@{$row}[0]]} = split /, /, $termString;
	@{$domainToGo{'go'}[@{$row}[0]]} = split /, /, $idString;
}
close $safeFH;

# read assgined NMF ranks
print STDOUT "Reading assigned NMF gene localizations\n";
my %assignedNMFLocalization;
open my $assignedNmfFH, '<', $nmfAssignedFile or die "Could not open $nmfAssignedFile: $!";
my $nmfAssignedTSV = Text::CSV_XS->new({
  sep_char => "\t",
});
$nmfAssignedTSV->getline($assignedNmfFH); # discard header
while(my $row = $nmfAssignedTSV->getline($assignedNmfFH)) {
  $assignedNMFLocalization{@{$row}[0]} = @{$row}[1];
}
close $assignedNmfFH;

# read assgined SAFE ranks
print STDOUT "Reading assigned SAFE gene localizations\n";
my %assignedSAFELocalization;
open my $assignedSafeFH, '<', $safeAssignedFile or die "Could not open $safeAssignedFile: $!";
my $safeAssignedTSV = Text::CSV_XS->new({
  binary => 1,
  sep_char => "\t",
  quote_char => undef,
  escape_char => undef,
  allow_loose_quotes => 1,
  allow_loose_escapes => 1,
});
$safeAssignedTSV->getline($assignedSafeFH); # discard header
while(my $row = $safeAssignedTSV->getline($assignedSafeFH)) {
  $assignedSAFELocalization{@{$row}[0]} = @{$row}[2];
}
close $assignedSafeFH;

# compare NMF and SAFE
print STDOUT "Comparing NMF and SAFE\n";
my %details;
my $matchedAlready = 0;
my $nmfMatches = 0;
my $safeMatches = 0;
my $totalMatches = 0;
foreach my $gene (keys %assignedNMFLocalization) {
  # get NMF  info
  my $rank = $assignedNMFLocalization{$gene};
  my @nmfGo = @{$rankToGo{'go'}[$rank]};
  my @nmfTerm = @{$rankToGo{'term'}[$rank]};
  my @possibleNmfGo = @{$rankToGo{'go'}[$rank]};
  # get all NMF possibilities
  foreach my $term (@nmfGo) {
    if (exists $children{$term}) {
      push @possibleNmfGo, @{$children{$term}};
    }
  }
  # get SAFE info
  my $domain;
  my @safeGo;
  my @safeTerm;
  my @possibleSafeGo;
  if (exists $assignedSAFELocalization{$gene}) {
    $domain = $assignedSAFELocalization{$gene};
    @safeGo = @{$domainToGo{'go'}[$domain]};
    @safeTerm = @{$domainToGo{'term'}[$domain]};
    @possibleSafeGo = @{$domainToGo{'go'}[$domain]};
    # get all SAFE possibilities
    foreach my $term (@safeGo) {
      if (exists $children{$term}) {
        push @possibleSafeGo, @{$children{$term}};
      }
    }
  }
  # is NMF term in SAFE?
  my @nmfIsect = intersect(@nmfGo, @possibleSafeGo);
  # is SAFE in NMF?
  my @safeIsect = intersect(@safeGo, @possibleNmfGo);
  # details
  $details{$gene}{'rank'} = $rank;
  $details{$gene}{'nmfTerm'} = join ', ', @nmfTerm;
  if ($domain) {
    $details{$gene}{'domain'} = $domain;
  } else {
    $details{$gene}{'domain'} = '';
  }
  if (scalar @safeTerm) {
    $details{$gene}{'safeTerm'} = join ', ', @safeTerm;
  } else {
    $details{$gene}{'safeTerm'} = '';
  }
  if (scalar @nmfIsect > 0) {
    $nmfMatches++;
    $totalMatches++;
    $matchedAlready = 1;
    $details{$gene}{'nmfAgree'} = 'true';
  } else {
    $details{$gene}{'nmfAgree'} = 'false';
  }
  if (scalar @safeIsect > 0) {
    $safeMatches++;
    if (!$matchedAlready) {
      $totalMatches++;
    }
    $details{$gene}{'safeAgree'} = 'true';
  } else {
    $details{$gene}{'safeAgree'} = 'false';
  }
  $matchedAlready = 0;
}

# print output
open my $outputFH, '>', 'nmf-v-safe.txt';
print $outputFH "gene\trank\tNMF term(s)\tdomain\tSAFE term(s)\tNMF in SAFE?\tSAFE in NMF?\n";
foreach my $gene (sort keys %details) {
  print $outputFH "$gene\t$details{$gene}{'rank'}\t$details{$gene}{'nmfTerm'}\t$details{$gene}{'domain'}\t$details{$gene}{'safeTerm'}\t$details{$gene}{'nmfAgree'}\t$details{$gene}{'safeAgree'}\n";
}
close $outputFH;

# print summary
my $nmfGenes = scalar (keys %assignedNMFLocalization);
my $safeGenes = scalar (keys %assignedSAFELocalization);
my $totalGenes;
if ($nmfGenes > $safeGenes) {
  $totalGenes = $nmfGenes;
} else {
  $totalGenes = $safeGenes;
}
my $nmfFraction = sprintf "%.2f", 100 * $nmfMatches / $nmfGenes;
my $safeFraction = sprintf "%.2f", 100 * $safeMatches / $safeGenes;
my $totalFraction = sprintf "%.2f", 100 * $totalMatches / $totalGenes;
print STDOUT "NMF: $nmfMatches of $nmfGenes, $nmfFraction\n";
print STDOUT "SAFE: $safeMatches of $safeGenes, $safeFraction\n";
print STDOUT "Total: $totalMatches of $totalGenes, $totalFraction\n";
