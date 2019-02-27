#!/usr/bin/perl

# 27/2/2019

use strict;
use warnings;

# libraries
use FindBin;
use lib "$FindBin::RealBin/../lib";

use Go::Children qw(readChildren);
use Go::Gaf qw(readGaf);
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# parameters
my $namespace = 'cc';

# command line parameters
my $cfile = ''; # list with child terms for each GO category
my $gfile = ''; # tsv file (with header), three columns with gene, its rank and score
my $lfile = ''; # GO localization file (goa_human.gaf without first 4 header lines)
my $rfile = ''; # file with rank details
my $threshold = 0;

if ($#ARGV==0){
	print "\nTakes a list of GO child terms, a list of official GO annotations for genes,\n";
	print "assigned localization form NMF and assesses the calculates the recovery of known\n";
	print "localization per NMF score.\n";
	print "\nusage:\n $0\n";
  print "-c [child terms for each GO term]\n";
	print "-g [gene localization file]\n";
	print "-l [GO annotations file (goa_human.gaf)]\n";
  print "-n [GO namespace, (default: cc, [options: bp, cc, mf])]\n";
	print "-r [rank/details]\n";
  print "-t [threshold]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-c'){
      $i++;
      $cfile = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-n'){
      $i++;
      $namespace = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-r'){
			$i++;
			$rfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
      $i++;
      $threshold = $ARGV[$i];
    } else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# set namespace
my $shortNameSpace;
my $longNameSpace;
if($namespace eq "bp") {
	$shortNameSpace = "P";
	$longNameSpace = "biological_process";
} elsif($namespace eq "cc") {
	$shortNameSpace = "C";
	$longNameSpace = "cellular_component";
} elsif($namespace eq "mf") {
	$shortNameSpace = "F";
	$longNameSpace = "molecular_function";
} else {
}

# read child terms for each GO term
my $childTerms = readChildren $cfile;
my %children = %{$childTerms};

# reading known localizations
my $gafAnnotations = readGaf($lfile, $shortNameSpace);
my %geneAnnotations = %{$gafAnnotations};

# read rank details
my %rankDetails;
print STDERR "Reading rank details\n";
my $rankTSV = Text::CSV_XS->new({
  binary => 1,
  sep_char => "\t",
  escape_char => undef,
  allow_loose_escapes => 1,
});
open my $rankfh, '<', $rfile or die "Could not open $rfile: $!";
$rankTSV->getline($rankfh); # discard header
while(my $row = $rankTSV->getline($rankfh)) {
  my $idString = @{$row}[3];
  $idString =~ s/[\[\]"]//g;
  my $termString = @{$row}[1];
  $termString =~ s/[\[\]"]//g;
  @{$rankDetails{'term'}[@{$row}[0]]} = split /, /, $termString;
  @{$rankDetails{'go'}[@{$row}[0]]} = split /, /, $idString;
}
close $rankfh;

# read assigned localizations and scores
print STDERR "Reading assigned gene localizations\n";
my %actualLocalizationName;
my %assignedLocalization;
my %indexMap;
my @localizationsPerGene;
open my $assignedlocalizationfh, '<', $gfile or die "Could not open $gfile: $!";
my $localizationTSV = Text::CSV_XS->new({ sep_char => "\t" });
$localizationTSV->getline($assignedlocalizationfh); # discard header
while(my $row = $localizationTSV->getline($assignedlocalizationfh)) {
  my $gene = @{$row}[0];
  $assignedLocalization{$gene}{'rank'} = @{$row}[1];
  $assignedLocalization{$gene}{'score'} = @{$row}[2];
  $assignedLocalization{$gene}{'known'} = 0;
}
close $assignedlocalizationfh;

# initialize recovery arrays
my @recoveryHundredth;
my @thresholdTotalHundredth;
for (my $i = 0; $i <= 100; $i++) {
  $recoveryHundredth[$i] = 0;
  $thresholdTotalHundredth[$i] = 0;
}

my @recoveryTwentieth;
my @thresholdTotalTwentieth;
for (my $i = 0; $i <= 20; $i++) {
  $recoveryTwentieth[$i] = 0;
  $thresholdTotalTwentieth[$i] = 0;
}

my @recoveryTenth;
my @thresholdTotalTenth;
for (my $i = 0; $i <= 10; $i++) {
  $recoveryTenth[$i] = 0;
  $thresholdTotalTenth[$i] = 0;
}

# assess localization
print STDERR "Assessing assigned localizations\n";
foreach my $gene (keys %assignedLocalization) {
  my $rank = $assignedLocalization{$gene}{'rank'};

  # get all possible GO terms for rank
	my @possibleTerms = @{$rankDetails{'go'}[$rank]};
	foreach my $term (@possibleTerms) {
		if (exists $children{$term}) {
			push @possibleTerms, @{$children{$term}};
		}
	}
	@possibleTerms = uniq @possibleTerms;

  # check if any of the possible rank terms is known for the gene
	if (scalar @possibleTerms > 0 && exists $geneAnnotations{$gene}) {
		foreach my $term (@{$geneAnnotations{$gene}}) {
			if (grep /^$term$/, @possibleTerms) {
				$assignedLocalization{$gene}{'known'} = 1;
				last;
			}
		}
	}

  # determine rank's array index
  my $indexHundredth = int(100 * $assignedLocalization{$gene}{'score'});
  my $indexTwentieth = int(20 * $assignedLocalization{$gene}{'score'});
  my $indexTenth = int(10 * $assignedLocalization{$gene}{'score'});
  if ($indexHundredth > 100) {
    $indexHundredth = 100;
  }
  if ($indexTwentieth > 20) {
    $indexTwentieth = 20;
  }
  if ($indexTenth > 10) {
    $indexTenth = 10;
  }
  $thresholdTotalHundredth[$indexHundredth]++;
  $thresholdTotalTwentieth[$indexTwentieth]++;
  $thresholdTotalTenth[$indexTenth]++;
  if ($assignedLocalization{$gene}{'known'}) {
    $recoveryHundredth[$indexHundredth]++;
    $recoveryTwentieth[$indexTwentieth]++;
    $recoveryTenth[$indexTenth]++;
  }
}

# output recovery
open my $outfile, '>', 'nmf-recovery.txt';
print $outfile "score\tfraction known\ttotal\trecovered\n";

# calculate recovery at one tenth binds
for (my $i = 10; $i >= 0; $i--) {
  my $total = $thresholdTotalTenth[$i];
  my $recovered = $recoveryTenth[$i];
  my $score = sprintf("%.2f", $i / 10);
  my $frac = 0;
  if ($total > 0) {
    $frac = sprintf("%.2f", $recovered / $total);
  }
  print $outfile "$score\t$frac\t$total\t$recovered\n";
}

print $outfile "\n";

# calculate recovery at one twentieth binds
for (my $i = 20; $i >= 0; $i--) {
  my $total = $thresholdTotalTwentieth[$i];
  my $recovered = $recoveryTwentieth[$i];
  my $score = sprintf("%.2f", $i / 20);
  my $frac = 0;
  if ($total > 0) {
    $frac = sprintf("%.2f", $recovered / $total);
  }
  print $outfile "$score\t$frac\t$total\t$recovered\n";
}

print $outfile "\n";

# calculate recovery at one hundreth binds
for (my $i = 100; $i >= 0; $i--) {
  my $total = $thresholdTotalHundredth[$i];
  my $recovered = $recoveryHundredth[$i];
  my $score = sprintf("%.2f", $i / 100);
  my $frac = 0;
  if ($total > 0) {
    $frac = sprintf("%.2f", $recovered / $total);
  }
  print $outfile "$score\t$frac\t$total\t$recovered\n";
}

close $outfile;

# output confidence and known status
open my $confidenceOutfile, '>', 'prey-confidence.txt';
print $confidenceOutfile "gene\trank\tscore\tknown\tconfident\n";
foreach my $gene (sort {lc $a cmp lc $b} keys %assignedLocalization) {
  my $known = $assignedLocalization{$gene}{'known'};
  my $rank = $assignedLocalization{$gene}{'rank'};
  my $score = $assignedLocalization{$gene}{'score'};
  my $confident = 0;
  if ($score >= $threshold) {
    $confident = 1;
  }
  print $confidenceOutfile "$gene\t$rank\t$score\t$known\t$confident\n";
}
close $confidenceOutfile;
