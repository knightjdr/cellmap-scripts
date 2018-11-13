#!/usr/bin/perl

# 31/5/2017

use strict;
use warnings;

# libraries
use Algorithm::MedianSelect qw(median);
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(uniq);
use List::Util qw(sum);
use Text::CSV_XS;

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
print STDERR "Creating map\n";
my %baitMap;
my @baits;
my $tsv = Text::CSV_XS->new({
  sep_char => "\t"
});
open my $fh, '<', $mfile or die "Could not open $mfile: $!";
$tsv->getline($fh); #discard header
while(my $row = $tsv->getline($fh)) {
	my $bait = lc @{$row}[0];
	my $gene = lc @{$row}[1];
	$baitMap{$bait} = $gene;
  push @baits, lc @{$row}[1];
}
close($fh);
@baits = uniq @baits;
my %baitsHash = map { $_ => 1 } @baits;

# parse BioGRID file
print STDERR "Getting list of interactors\n";
my %biogrid;
my $biogridTSV = Text::CSV_XS->new({
  sep_char => "\t"
});
open my $biogridFH, '<', $bfile or die "Could not open $bfile: $!";
$biogridTSV->getline($biogridFH); #discard header
if ($fileType eq 'i') {
	while(my $row = $biogridTSV->getline($biogridFH)) {
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
	while(my $row = $tsv->getline($biogridFH)) {
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
  while(my $row = $biogridTSV->getline($biogridFH)) {
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
close($biogridFH);

# get SAINT data
print STDERR "Reading SAINT file\n";
my %saintInteractions;
my $saintTSV = Text::CSV_XS->new({
  sep_char => "\t"
});
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
$saintTSV->getline($saintFH); #discard header
while(my $row = $saintTSV->getline($saintFH)) {
  my $fdr = @{$row}[15];
  my $avgSpec = @{$row}[5];
  if ($controlSubtraction) {
    my @controlArray = split '\|', @{$row}[7];
    my $controlAvg = sum(@controlArray) / @controlArray;
    $avgSpec -= $controlAvg;
  }
	my $bait = @{$row}[0];
  my $length = @{$row}[20];
  my $prey = @{$row}[2];
  if ($fdr <= $score) {
    $saintInteractions{$bait}{$prey}{'avgspec'} = $avgSpec;
    $saintInteractions{$bait}{$prey}{'length'} = $length;
  }
}
close $saintFH;

# sort baits based on prey spectral counts, and only keep top X
my %topPreys;
foreach my $bait (keys %saintInteractions) {
  # adjust prey spectral count to length
  my @lengths;
  foreach my $prey (keys %{$saintInteractions{$bait}}) {
    push @lengths, $saintInteractions{$bait}{$prey}{'length'};
  }
  my $median = median @lengths;
  if ($adjust) {
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
  $topPreys{$bait}{'discovered'} = 0;
  $topPreys{$bait}{'recovered'} = 0;
}

# check SAINT pairs
print STDERR "Checking interactions\n";
my %recovered;
foreach my $bait (keys %topPreys) {
  foreach my $prey (@{$topPreys{$bait}{'list'}}) {
    my %currHash = map { $_ => 1 } @{$biogrid{$baitMap{lc $bait}}};
    if (exists $currHash{lc $prey}) {
      $topPreys{$bait}{'recovered'}++;
      push @{$recovered{$bait}}, $prey;
    } else {
      $topPreys{$bait}{'discovered'}++;
    }
  }
}

# print bait results
my $fileHandle = 'bait-level-recovery_top' . $topX;
if ($adjust) {
  $fileHandle .= '_lengthAdjusted'
}
$fileHandle .= '.txt';
open my $resultsfh, '>', $fileHandle;
print $resultsfh "bait\tunknown\tknown\tfraction\trecovered preys\ttop preys\n";
foreach my $bait (keys %topPreys) {
  my $recoveredString = '-';
  if (exists $recovered{$bait}) {
    $recoveredString = join ', ', sort @{$recovered{$bait}};
  }
  my $preyString = join ', ', sort @{$topPreys{$bait}{'list'}};
  my $fraction = sprintf "%.3f", ($topPreys{$bait}{'recovered'} / ($topPreys{$bait}{'recovered'} + $topPreys{$bait}{'discovered'}));
  print $resultsfh "$bait\t$topPreys{$bait}{'discovered'}\t$topPreys{$bait}{'recovered'}\t$fraction\t$recoveredString\t$preyString\n";
}
close $resultsfh;
