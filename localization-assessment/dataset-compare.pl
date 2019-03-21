#!/usr/bin/perl

# 21/3/2019

use strict;
use warnings;

# libraries
use FindBin;
use lib "$FindBin::RealBin/../lib"; 

use Array::Utils qw(:all);
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(uniq);
use Localization::ParseSummary qw(parseSummary);
use Text::CSV_XS;

# command line parameters
my $gfile = ''; # file with list of children per GO term
my $lfile = ''; # goa_human annotation file
my $mfile = ''; # map of GO ID to term
my $nfile = ''; # nmf summary file
my $sfile = ''; # safe summary file
my $file1 = ''; # file with localization for each gene in dataset 1
my $file2 = ''; # file with localization for each gene in dataset 2
my $file1Type = '';
my $file2Type = '';
my $outPrefix = 'comparison'; # outfile prefix

if ($#ARGV==0) {
	print "\nTakes two datasets and summarizes their localization overlap.\n";
	print "\nusage:\n $0\n";
	print "-g [File with list of children for every GO term (from get_children.pl)]\n";
  print "-l [GO localization file (goa_human.gaf)]\n";
  print "-m [map of GO ID to term]\n";
  print "-n [nmf_summary.txt]\n";
	print "-s [safe_summary.txt]\n";
  print "-f1 [file with first dataset]\n";
  print "-f2 [file with second dataset]\n";
  print "-f1t [first file type: n = NMF, s = SAFE, h = HPA, o = Other]\n";
  print "-f2t [second file type]\n";
  print "-o [output file prefix]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV) {
		if ($ARGV[$i] eq '-g') {
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-l') {
      $i++;
      $lfile = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-m') {
      $i++;
      $mfile = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-n') {
			$i++;
			$nfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s') {
			$i++;
			$sfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-f1') {
      $i++;
      $file1 = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-f2') {
			$i++;
			$file2 = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-f1t') {
      $i++;
      $file1Type = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-f2t') {
      $i++;
      $file2Type = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-o') {
      $i++;
      $outPrefix = $ARGV[$i];
    } else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

my %tsvParams = (
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
);

sub allTerms {
  my @intialAnnotations = @{$_[0]};
  my %children = %{$_[1]};

  my @allAnnotations = @intialAnnotations;
  foreach my $term (@intialAnnotations) {
    if (exists $children{$term}) {
      push @allAnnotations, @{$children{$term}};
    }
  }
  my @uniqueAnnotations = uniq(@allAnnotations);

  return \@uniqueAnnotations;
}

sub mapIDs {
  my @ids = @{$_[0]};
  my %map = %{$_[1]};

  my @mapped;

  foreach my $id (@ids) {
    push @mapped, $map{$id}
  }

  return join '; ', @mapped;
}

sub parseLocalization {
  my $localization = $_[0];
  my $fileType = $_[1];
  my %ranks = %{$_[2]};
  my %domains = %{$_[3]};

  my @localizations;
  if ($fileType eq 'n') {
    @localizations = @{$ranks{'go'}[$localization]};
  } elsif ($fileType eq 's'){
    @localizations = @{$domains{'go'}[$localization]};
  } elsif ($fileType eq 'h'){
    my @goAnnotations = split /;/, $localization;
    foreach my $annotation (@goAnnotations) {
      my ($goID) = $annotation =~ /[^\(]+\((GO:\d+)?/;
      if ($goID) {
        push @localizations, $goID;
      }
    }
  }else {}

  return \@localizations;
}

sub readLocalizations {
  my $file = $_[0];
  my $fileType = $_[1];
  my %ranks = %{$_[2]};
  my %domains = %{$_[3]};

  # Get columns to read gene and localization from.
  my $geneColumn;
  my $localizationColumn;
  if ($fileType eq 'n') {
    $geneColumn = 0;
    $localizationColumn = 1;
  } elsif ($fileType eq 'h'){
    $geneColumn = 1;
    $localizationColumn = 10;
  } elsif ($fileType eq 's'){
    $geneColumn = 0;
    $localizationColumn = 2;
  } else {
    $geneColumn = 0;
    $localizationColumn = 1;
  }

  my %localizations;

  open my $fh, '<', $file or die "Could not open $file: $!";
  my $tsv = Text::CSV_XS->new(\%tsvParams);
  $tsv->getline($fh); # discard header
  while(my $row = $tsv->getline($fh)) {
    my $gene = @{$row}[$geneColumn];
    my $localization = @{$row}[$localizationColumn];
    $localizations{$gene} = parseLocalization($localization, $fileType, \%ranks, \%domains);
  }
  close $fh;

  return \%localizations;
}

# get child terms for each GO term
print STDOUT "Reading hierarchy\n";
open my $goChildrenFH, '<', $gfile or die "Could not open $gfile: $!";
my $goChildrenTSV = Text::CSV_XS->new(\%tsvParams);
my %children;
while(my $row = $goChildrenTSV->getline($goChildrenFH)) {
	@{$children{@{$row}[0]}} = split ',', @{$row}[1];
}
close $goChildrenFH;

# Read map of GO ID to term
print STDOUT "Reading GO map\n";
open my $goMapFH, '<', $mfile or die "Could not open $mfile: $!";
my $goMapTSV = Text::CSV_XS->new(\%tsvParams);
my %map;
while(my $row = $goMapTSV->getline($goMapFH)) {
  my $id = @{$row}[0];
  my $term = @{$row}[1];
	$map{$id} = $term;
}
close $goMapFH;

# reading known localizations
print STDERR "Reading known localizations\n";
open my $localizationFH, '<', $lfile or die "Could not open $lfile: $!";
my $localizationTSV = Text::CSV_XS->new(\%tsvParams);
my %geneLocalizations;
while(my $row = $localizationTSV->getline($localizationFH)) {
  my $namespace = @{$row}[8] ;
	if ($namespace eq 'C') {
		push @{$geneLocalizations{@{$row}[2]}}, @{$row}[4];
	}
}
close $localizationFH;

# read nmf summary terms
my %rankToGo;
if ($nfile) {
  print STDOUT "Reading NMF details\n";
  my ($noRanks, $nmfMapToGo) = parseSummary($nfile);
  %rankToGo = %{$nmfMapToGo};
}

# read SAFE summary terms
my %domainToGo;
if ($sfile) {
  print STDOUT "Reading SAFE details\n";
  my ($noDomains, $safeMapToGo) = parseSummary($sfile);
  %domainToGo = %{$safeMapToGo};
}

# read first dataset terms
print STDOUT "Reading first dataset\n";
my ($localizationsFile1) = readLocalizations($file1, $file1Type, \%rankToGo, \%domainToGo);
my %localizationsFirst = %{$localizationsFile1};

# read second dataset terms
print STDOUT "Reading second dataset\n";
my ($localizationsFile2) = readLocalizations($file2, $file2Type, \%rankToGo, \%domainToGo);
my %localizationsSecond = %{$localizationsFile2};

# Compare datasets
print STDOUT "Comparing datasets\n";
open my $outputFH, '>', $outPrefix . '.txt';
print $outputFH "gene\tdataset 1 term(s)\tdataset 2 term(s)\tdataset 1 known\tdataset 2 known\tmatch\n";
my $matches = 0;
my $totalGenes = 0;
my $unmatchedKnownFirst = 0;
my $unmatchedKnownSecond = 0;
foreach my $gene (sort {lc $a cmp lc $b} keys %localizationsFirst) {
  my $matched = 'FALSE';
  if (
    scalar @{$localizationsFirst{$gene}} > 0
    && exists $localizationsSecond{$gene}
    && scalar @{$localizationsSecond{$gene}} > 0
  ) {
    $totalGenes++;
    my @childAnnotationsFirst = @{allTerms(\@{$localizationsFirst{$gene}}, \%children)};
    my @childAnnotationsSecond = @{allTerms(\@{$localizationsSecond{$gene}}, \%children)};

    my @isectFirst = intersect(@{$localizationsFirst{$gene}}, @childAnnotationsSecond);
    my @isectSecond = intersect(@{$localizationsSecond{$gene}}, @childAnnotationsFirst);

    if (
      scalar @isectFirst > 0
      || scalar @isectSecond > 0
    ) {
      $matches++;
      $matched = 'TRUE'
    }

    # Check if localization is previously known.
    my $knownFirst = 'FALSE';
    my $knownSecond = 'FALSE';
    my $isectFirst = intersect(@childAnnotationsFirst, @{$geneLocalizations{$gene}});
    my $isectSecond = intersect(@childAnnotationsSecond, @{$geneLocalizations{$gene}});
    if ($isectFirst > 0) {
      $knownFirst = 'TRUE';
    }
    if ($isectSecond > 0) {
      $knownSecond = 'TRUE';
    }

    # Count unmatched hits that have known localizations
    if (
      $matched eq 'FALSE'
      && $isectFirst > 0
    ) {
      $unmatchedKnownFirst++;
    }
    if (
      $matched eq 'FALSE'
      && $isectSecond > 0
    ) {
      $unmatchedKnownSecond++;
    }

    # Map GO IDs to terms
    my $termsFirst = mapIDs(\@{$localizationsFirst{$gene}}, \%map);
    my $termsSecond = mapIDs(\@{$localizationsSecond{$gene}}, \%map);
    print $outputFH "$gene\t$termsFirst\t$termsSecond\t$knownFirst\t$knownSecond\t$matched\n";
  }
}
close $outputFH;

my $frac = sprintf "%.3f", $matches / $totalGenes;
print STDOUT "Matches: $matches / $totalGenes ($frac)\n";

my $unmatched = $totalGenes - $matches;
my $unknownFracFirst = sprintf "%.3f", $unmatchedKnownFirst / $unmatched;
my $unknownFracSecond = sprintf "%.3f", $unmatchedKnownSecond / $unmatched;
print STDOUT "Unmatched genes that have known localizations in dataset 1: $unmatchedKnownFirst / $unmatched ($unknownFracFirst)\n";
print STDOUT "Unmatched genes that have known localizations in dataset 2: $unmatchedKnownSecond / $unmatched ($unknownFracSecond)\n";