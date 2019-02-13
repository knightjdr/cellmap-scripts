#!/usr/bin/perl

# 31/5/2017

use strict;
use warnings;

# libraries
use FindBin;
use lib "$FindBin::RealBin/../lib"; 

use Bait::Gene qw(parseMap);
use Interactions::Parse qw(readInteractions);
use Interactions::Saint qw(readSaintInteractions);
use Interactions::Summarize qw(summarize);

#parameters
my $fileType = 'b';
my $score = 0.01; # FDR

# command line parameters
my $bfile = ''; # BioGRID file
my $mfile = '';	# file with a map of Bait names to gene names "bait.dat" from ProHITS
my $sfile = '';	# SAINT file

if ($#ARGV==0){
  print "\nTakes a SAINT file and a BioGRID/IntAct file, and for each bait it calculates the number of\n";
  print "newly discovered preys, number of discovered previsouly known preys, and the number of\n";
  print "previously known preys that were missed. Also requires a file for mapping bait names to\n";
  print "gene names. \n\n";
	print "\nusage:\n $0\n";
	print "-b [BioGRID/IntAct file]\n";
  print "-f [FDR cutoff]\n";
	print "-m [SAINT bait-gene map]\n";
	print "-s [SAINT file]\n";
	print "-t [file type: b = BioGRID (default), i = IntAct, m = merged\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
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
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read bait name to gene map
my ($baitGeneMap, $geneBaitMap) = parseMap $mfile;
my %baitMap = %{ $baitGeneMap };
my %geneMap = %{ $geneBaitMap };

# parse Interaction file
my $parsedInteractions = readInteractions($bfile, \%geneMap);
my %interactions = %{$parsedInteractions};

# parse SAINT Interactions
my $parsedSaintInteractions = readSaintInteractions($sfile, 0.01, 0);
my %saintInteractions = %{$parsedSaintInteractions};

# summarizing edges
print STDERR "Summarizing edges\n";
open my $edgefh, '>', 'edge-recovered.txt';
print $edgefh "name\tknown\n";
foreach my $bait (keys %saintInteractions) {
  my %knownPreys = map { $_ => 1 } @{$interactions{lc $baitMap{$bait}}};
  foreach my $prey (keys %{$saintInteractions{$bait}}) {
    if (exists $knownPreys{lc $prey}) {
      print $edgefh "$baitMap{$bait} (interacts with) $prey\t1\n";
    }
  }
}
close $edgefh;

# check SAINT pairs
my ($recoveredInteractions, $interactionsStats) = summarize(\%saintInteractions, \%interactions, \%baitMap);
my %recovered = %{$recoveredInteractions};
my %stats = %{$interactionsStats};

# output stats
open my $resultsfh, '>', 'bait-level-recovery.txt';
print $resultsfh "bait\tdiscovered\trecovered\tmissed\trecovered preys\n";
my $known = 0;
my $total = 0;
foreach my $bait (keys %saintInteractions) {
  my $recoveredString = '-';
  if (exists $recovered{$bait}) {
    $recoveredString = join ', ', sort @{$recovered{$bait}};
  }
  print $resultsfh "$bait\t$stats{$bait}{'discovered'}\t$stats{$bait}{'recovered'}\t$stats{$bait}{'missed'}\t$recoveredString\n";
  $known += $stats{$bait}{'recovered'};
  $total += ($stats{$bait}{'discovered'} + $stats{$bait}{'recovered'});
}
close $resultsfh;

# print summary
my $percentage = 100 * $known / $total;
print STDOUT "known: $known\ntotal: $total\n$percentage: (100*$known/$total)\n";
