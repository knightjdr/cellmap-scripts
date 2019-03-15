#!/usr/bin/perl

# 13/6/2018

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# Variables
my $secondaryThreshold = 0.75;
my $shortNameSpace = "C";
my $longNameSpace = "cellular_component";
my @excludeLocalizations = [
  'GO:0005575', # 'cellular_component',
  'GO:0005623', # 'cell',
  'GO:0044464', # 'cell part',
  'GO:0005622', # 'intracellular',
  'GO:0044424', # 'intracellular part',
  'GO:0043226', # 'organelle',
  'GO:0043229', # 'intracellular organelle',
  'GO:0044422', # 'organelle part',
];

# command line parameters
my $bfile = ''; # basis.csv file
my $hfile = ''; # hierarchy file
my $lfile = ''; # GO localization file (goa_human.gaf)
my $rfile = ''; # for NMF the file with rank details

if ($#ARGV==0){
	print "\nTakes GO hierarchy from get_children.pl, a list of official localizations for genes,\n";
	print "an NMF basis matrix and reports on the number of genes whose primary and secondary\n";
  print "assignments are previously known.\n";
	print "\nusage:\n $0\n";
	print "-b [basis.csv]\n";
  print "-h [hierarcy file procduced by get_children.pl]\n";
	print "-l [GO localization file (goa_human.gaf)]\n";
	print "-r [file with rank details]\n";
  print "-s [secondary threshold]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
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
      $secondaryThreshold = $ARGV[$i];
    } else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# get child terms for each GO term
print STDERR "Reading hierarchy file\n";
open my $hierarchyFH, '<', $hfile or die "Could not open $hfile: $!";
my $tsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
my %children;
while(my $row = $tsv->getline($hierarchyFH)) {
	@{$children{@{$row}[0]}} = split ',', @{$row}[1];
}
close $hierarchyFH;

# reading known localizations
print STDERR "Reading known localizations\n";
$tsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
open my $localizationfh, '<', $lfile or die "Could not open $lfile: $!";
my %geneLocalizations;
while(my $row = $tsv->getline($localizationfh)) {
	if (@{$row}[8] eq $shortNameSpace) {
		push @{$geneLocalizations{@{$row}[2]}}, @{$row}[4];
	}
}
close $localizationfh;

# get rank details
my $noRanks = 0;
my %rankToGo;
print STDERR "Reading rank details\n";
$tsv = Text::CSV_XS->new({
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

# Get primary and secondary gene localizations.
print STDERR "Reading assigned gene localizations\n";
my %primaryLocalization;
my %secondaryLocalization;
open my $assignedlocalizationfh, '<', $bfile or die "Could not open $bfile: $!";
$tsv = Text::CSV_XS->new({
  sep_char => ",",
});
$tsv->getline($assignedlocalizationfh); # discard header
while(my $row = $tsv->getline($assignedlocalizationfh)) {
  my $gene = @{$row}[0];
  my @values = splice @{$row}, 1;
  my $max = 0;
  my $primary = -1;
  for (my $i = 0, my $iLen = scalar(@values); $i < $iLen; $i++) {
    if ($values[$i] > $max) {
      $max = $values[$i];
      $primary = $i;
    }
  }
  my $secondary = -1;
  my $secondaryMax = 0;
  my $secondaryMin = $max * $secondaryThreshold;
  for (my $i = 0, my $iLen = scalar(@values); $i < $iLen; $i++) {
    if (
      $i != $primary &&
      $values[$i] >= $secondaryMin &&
      $values[$i] > $secondaryMax
    ) {
      $secondaryMax = $values[$i];
      $secondary = $i;
    }
  }
  $primaryLocalization{$gene} = $primary + 1;
  $secondaryLocalization{$gene} = $secondary + 1;
}
close $assignedlocalizationfh;

# check for genes with unknown or excluded localizations.
print STDERR "Checking for unknown localizations\n";
my $unknown = 0; # Count genes with no known localizations
while ( my( $gene ) = each %primaryLocalization ) {
  if (!(exists $geneLocalizations{$gene})) {
    $unknown++;
  } else {
    my $notExcluded = 0;
    foreach my $term (@{$geneLocalizations{$gene}}) {
      if (!(grep /^$term$/, @excludeLocalizations)) {
        $notExcluded = 1;
        last;
      }
    }
    if ($notExcluded == 0) {
      $unknown++;
    }
  }
}
print STDOUT "genes with no known localization: $unknown\n";

# Assess primary localization
print STDERR "Assessing localizations\n";
my %correctAssignments = (
	'primary' => 0,
  'secondary' => 0,
  'total' => 0,
);
my %totalGenes = (
  'primary' => 0,
  'secondary' => 0,
  'total' => 0,
);
my %omit; # genes with a correct primary localization. Want to omit from secondary search.
while ( my( $gene, $rank ) = each %primaryLocalization ) {
  $totalGenes{'primary'}++;
  $totalGenes{'total'}++;
  my @possibleTerms = @{$rankToGo{'go'}[$rank]};
  my $noInitialTerms = scalar @possibleTerms;
  foreach my $term (@possibleTerms) {
    if (exists $children{$term}) {
      push @possibleTerms, @{$children{$term}};
    }
  }
  @possibleTerms = uniq @possibleTerms;
  if (scalar @possibleTerms > 0 && exists $geneLocalizations{$gene}) {
    foreach my $term (@{$geneLocalizations{$gene}}) {
      if (grep /^$term$/, @possibleTerms) {
        $correctAssignments{'primary'}++;
        $correctAssignments{'total'}++;
        $omit{$gene} = 1;
        last;
      }
    }
  }
}
my $primaryFrac = sprintf "%.3f", $correctAssignments{'primary'} / $totalGenes{'primary'};
my $primaryFracWoUnknown = sprintf "%.3f", $correctAssignments{'primary'} / ($totalGenes{'primary'} - $unknown);
print STDOUT "primary\t$primaryFrac\t$totalGenes{'primary'}\t$correctAssignments{'primary'}\two unknown: $primaryFracWoUnknown\n";

# Assess secondary localization
while ( my( $gene, $rank ) = each %secondaryLocalization ) {
  if ($rank > 0) {
    $totalGenes{'secondary'}++;
    my @possibleTerms = @{$rankToGo{'go'}[$rank]};
    my $noInitialTerms = scalar @possibleTerms;
    foreach my $term (@possibleTerms) {
      if (exists $children{$term}) {
        push @possibleTerms, @{$children{$term}};
      }
    }
    @possibleTerms = uniq @possibleTerms;
    if (scalar @possibleTerms > 0 && exists $geneLocalizations{$gene}) {
      foreach my $term (@{$geneLocalizations{$gene}}) {
        if (grep /^$term$/, @possibleTerms) {
          $correctAssignments{'secondary'}++;
          if (!(exists $omit{$gene})) {
            $correctAssignments{'total'}++;
          }
          last;
        }
      }
    }
  }
}
my $secondaryFrac = sprintf "%.3f", $correctAssignments{'secondary'} / $totalGenes{'secondary'};
my $totalFrac = sprintf "%.3f", $correctAssignments{'total'} / $totalGenes{'total'};
my $totalFracWoUnknown = sprintf "%.3f", $correctAssignments{'total'} / ($totalGenes{'total'} - $unknown);
print STDOUT "secondary\t$secondaryFrac\t$totalGenes{'secondary'}\t$correctAssignments{'secondary'}\n";
print STDOUT "total\t$totalFrac\t$totalGenes{'total'}\t$correctAssignments{'total'}\two unknown: $totalFracWoUnknown\n";
