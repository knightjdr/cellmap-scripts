#!/usr/bin/perl

# 31/5/2017

use strict;
use warnings;

# libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

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
print STDERR "Creating map\n";
my %baitMap;
my %baitActualMap;
my %baitLCtoUCMap;
my @baits;
my $tsv = Text::CSV_XS->new({ sep_char => "\t" });
open my $fh, '<', $mfile or die "Could not open $mfile: $!";
$tsv->getline($fh); #discard header
while(my $row = $tsv->getline($fh)) {
  $baitActualMap{@{$row}[0]} = @{$row}[1];
	my $bait = lc @{$row}[0];
	my $gene = lc @{$row}[1];
	$baitMap{$bait} = $gene;
  $baitLCtoUCMap{$bait} = @{$row}[0];
  push @baits, lc @{$row}[1];
}
close($fh);
@baits = uniq @baits;
my %baitsHash = map { $_ => 1 } @baits;

# parse BioGRID file
print STDERR "Getting list of interactors\n";
my %biogrid;
$tsv = Text::CSV_XS->new({ sep_char => "\t" });
open $fh, '<', $bfile or die "Could not open $bfile: $!";
$tsv->getline($fh); #discard header
if ($fileType eq 'i') {
	while(my $row = $tsv->getline($fh)) {
		my $source = lc @{$row}[0];
		my $target = lc @{$row}[1];
    if (exists $baitsHash{$source}) {
      my %currHash = map { $_ => 1 } @{$biogrid{$source}};
      if (!(exists $currHash{$target})) {
        push @{$biogrid{$source}}, $target;
      }
    }
    if (exists $baitsHash{$target}) {
      my %currHash = map { $_ => 1 } @{$biogrid{$target}};
      if (!(exists $currHash{$source})) {
        push @{$biogrid{$target}}, $source;
      }
    }
  }
} elsif ($fileType eq 'm') {
	while(my $row = $tsv->getline($fh)) {
		my $source = lc @{$row}[0];
		my $target = lc @{$row}[1];
    if (exists $baitsHash{$source}) {
      my %currHash = map { $_ => 1 } @{$biogrid{$source}};
      if (!(exists $currHash{$target})) {
        push @{$biogrid{$source}}, $target;
      }
    }
    if (exists $baitsHash{$target}) {
      my %currHash = map { $_ => 1 } @{$biogrid{$target}};
      if (!(exists $currHash{$source})) {
        push @{$biogrid{$target}}, $source;
      }
    }
  }
} else {
  while(my $row = $tsv->getline($fh)) {
    my $source = lc @{$row}[7];
    my $target = lc @{$row}[8];
    if (exists $baitsHash{$source}) {
      my %currHash = map { $_ => 1 } @{$biogrid{$source}};
      if (!(exists $currHash{$target})) {
        push @{$biogrid{$source}}, $target;
      }
    }
    if (exists $baitsHash{$target}) {
      my %currHash = map { $_ => 1 } @{$biogrid{$target}};
      if (!(exists $currHash{$source})) {
        push @{$biogrid{$target}}, $source;
      }
    }
  }
}
close($fh);

# check SAINT pairs
print STDERR "Reading SAINT file\n";
my %recovered;
my %saintInteractions;
$tsv = Text::CSV_XS->new({ sep_char => "\t" });
open $fh, '<', $sfile or die "Could not open $sfile: $!";
$tsv->getline($fh); #discard header
open my $edgefh, '>', 'edge-recovered.txt';
print $edgefh "name\tknown\n";
while(my $row = $tsv->getline($fh)) {
	my $bait = lc @{$row}[0];
  my $prey = lc @{$row}[2];
  my $fdr = @{$row}[15];
  if (!(exists $saintInteractions{$bait})) {
    $saintInteractions{$bait}{'discovered'} = 0;
    if (exists $biogrid{$baitMap{$bait}}) {
      $saintInteractions{$bait}{'missed'} = scalar @{$biogrid{$baitMap{$bait}}};
    } else {
      $saintInteractions{$bait}{'missed'} = 0;
    }
    $saintInteractions{$bait}{'recovered'} = 0;
  }
  if ($fdr <= $score) {
    my %currHash = map { $_ => 1 } @{$biogrid{$baitMap{$bait}}};
    if (exists $currHash{$prey}) {
      $saintInteractions{$bait}{'missed'}--;
      $saintInteractions{$bait}{'recovered'}++;
      print $edgefh "$baitActualMap{@{$row}[0]} (interacts with) @{$row}[2]\t1\n";
      push @{$recovered{$bait}}, @{$row}[2];
    } else {
      $saintInteractions{$bait}{'discovered'}++;
    }
  }
}
close $fh;

# print bait results
open my $resultsfh, '>', 'bait-level-recovery.txt';
print $resultsfh "bait\tdiscovered\trecovered\tmissed\trecovered preys\n";
foreach my $bait (keys %saintInteractions) {
  my $recoveredString = '-';
  if (exists $recovered{$bait}) {
    $recoveredString = join ', ', sort @{$recovered{$bait}};
  }
  print $resultsfh "$baitLCtoUCMap{$bait}\t$saintInteractions{$bait}{'discovered'}\t$saintInteractions{$bait}{'recovered'}\t$saintInteractions{$bait}{'missed'}\t$recoveredString\n";
}
close $resultsfh;
