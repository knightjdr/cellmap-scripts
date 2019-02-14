#!/usr/bin/perl

# 31/5/2017

use strict;
use warnings;

# libraries
use FindBin;
use lib "$FindBin::RealBin/../lib"; 

use Algorithm::MedianSelect qw(median);
use Bait::Gene qw(parseMap);
use Interactions::Parse qw(readInteractions);
use Interactions::Saint qw(readSaintInteractions);
use Interactions::Summarize qw(summarizeTop);

#parameters
my $adjust = 1; # if true, spectral counts will be adjusted to length
my $controlSubtraction = 1; # subtract average control value
my $fileType = 'b';
my $score = 0.01; # fdr
my $topX = 25;

# command line parameters
my $bfile = ''; # BioGRID file
my $mfile = '';	# file with a map of Bait names to gene names "bait.dat" from ProHITS
my $sfile = '';	# SAINT file

if ($#ARGV==0){
  print "\nTakes a SAINT file and a BioGRID/IntAct file, and for each bait it calculates the number of\n";
  print "newly discovered preys, number of discovered previsouly known preys, and the number of\n";
  print "previously known preys that were missed. Also requires a file for mapping bait names to\n";
  print "gene names. Unlike bait-summary.pl, this is for only the top X preys (X = 25 default).\n\n";
	print "\nusage:\n $0\n";
  print "-a [adjust spectral counts to length? [1 (default): true, 0: false]]\n";
	print "-b [BioGRID/IntAct file]\n";
  print "-c [should control avg be subtracted from AvgSpec: 1 (default) = true, 0 = false]";
  print "-f [FDR cutoff]\n";
	print "-m [SAINT bait-gene map]\n";
	print "-s [SAINT file]\n";
	print "-t [file type: b = BioGRID (default), i = IntAct, m = merged\n";
  print "-x [Number of preys to do this for (25 default)]\n";
	die "\n";
} else{
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-a'){
			$i++;
			$adjust = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-c'){
			$i++;
			$controlSubtraction = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-f'){
			$i++;
			$score = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-m'){
			$i++;
			$mfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
			$i++;
			$fileType = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-x'){
			$i++;
			$topX = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read bait name to gene map
my ($baitGeneMap, $geneBaitMap) = parseMap $mfile;
my %baitMap = %{$baitGeneMap};
my %geneMap = %{$geneBaitMap};

# parse Interaction file
my $parsedInteractions = readInteractions($bfile, \%geneMap);
my %interactions = %{$parsedInteractions};

# parse SAINT Interactions
my $parsedSaintInteractions = readSaintInteractions($sfile, 0.01, 1);
my %saintInteractions = %{$parsedSaintInteractions};

# sort preys for each bait based on spectral counts, and only keep top X
my %topPreys;
foreach my $bait (keys %saintInteractions) {
  if ($adjust) {
    # adjust prey spectral count to length
    my @lengths;
    foreach my $prey (keys %{$saintInteractions{$bait}}) {
      push @lengths, $saintInteractions{$bait}{$prey}{'length'};
    }
    my $median = median @lengths;
    @{$topPreys{$bait}{'list'}} = sort {
      $saintInteractions{$bait}{$b}{'avgspec'} * $median / $saintInteractions{$bait}{$b}{'length'} <=>
      $saintInteractions{$bait}{$a}{'avgspec'} * $median / $saintInteractions{$bait}{$a}{'length'}
      or
      lc $a cmp lc $b
    } keys %{$saintInteractions{$bait}};
  } else {
    @{$topPreys{$bait}{'list'}} = sort {
      $saintInteractions{$bait}{$b}{'avgspec'} <=>
      $saintInteractions{$bait}{$a}{'avgspec'}
      or
      lc $a cmp lc $b
    } keys %{$saintInteractions{$bait}};
  }
  my $limit = $topX;
  if (scalar @{$topPreys{$bait}{'list'}} < $limit) {
    $limit = scalar @{$topPreys{$bait}{'list'}};
  }
  @{$topPreys{$bait}{'list'}} = @{$topPreys{$bait}{'list'}}[0..($limit - 1)];
}

# check SAINT pairs
my ($recoveredInteractions, $interactionsStats) = summarizeTop(\%topPreys, \%interactions, \%baitMap);
my %recovered = %{$recoveredInteractions};
my %stats = %{$interactionsStats};

# output stats
my $fileHandle = 'bait-level-recovery_top' . $topX;
if ($adjust) {
  $fileHandle .= '_lengthAdjusted'
}
$fileHandle .= '.txt';
open my $resultsfh, '>', $fileHandle;
print $resultsfh "bait\tunknown\tknown\tfraction\trecovered preys\ttop preys\n";
my $known = 0;
my $total = 0;
foreach my $bait (keys %topPreys) {
  my $recoveredString = '-';
  if (exists $recovered{$bait}) {
    $recoveredString = join ', ', sort @{$recovered{$bait}};
  }
  my $preyString = join ', ', sort @{$topPreys{$bait}{'list'}};
  my $fraction = sprintf "%.3f", ($stats{$bait}{'recovered'} / ($stats{$bait}{'recovered'} + $stats{$bait}{'discovered'}));
  print $resultsfh "$bait\t$stats{$bait}{'discovered'}\t$stats{$bait}{'recovered'}\t$fraction\t$recoveredString\t$preyString\n";
  $known += $stats{$bait}{'recovered'};
  $total += ($stats{$bait}{'discovered'} + $stats{$bait}{'recovered'});
}
close $resultsfh;

# print summary
my $percentage = 100 * $known / $total;
print STDOUT "known: $known\ntotal: $total\npercentarge: $percentage\n";