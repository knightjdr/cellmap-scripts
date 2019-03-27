#!/usr/bin/perl

# 27/3/2019

use strict;
use warnings;

# libraries
use Array::Utils qw(:all);
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # bioplex file
my $fdr = 0.01;
my $minInteraction = 1;
my $sfile = ''; # saint file

if ($#ARGV==0) {
	print "\nComputes the overlap between SAINT file and BioPlex interaction list\n";
	print "\nusage:\n $0\n";
	print "-b [BioPlex file]\n";
  print "-f [SAINT prey FDR]\n";
	print "-m [minimum number of bioplex baits]\n";
	print "-s [saint file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-f'){
			$i++;
			$fdr = $ARGV[$i];
		} elsif($ARGV[$i] eq '-m'){
			$i++;
			$minInteraction = $ARGV[$i];
		} elsif($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else{
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

# read BioPlex file
print STDERR "Reading BioPlex file\n";
open my $bioplexFH, '<', $bfile or die "Could not open $bfile: $!";
my $bioplexTSV = Text::CSV_XS->new(\%tsvParams);
my %bioplex;
$bioplexTSV->getline($bioplexFH); # discard header
while(my $row = $bioplexTSV->getline($bioplexFH)) {
  my $source = @{$row}[4];
  my $target = @{$row}[5];
  if (exists $bioplex{$source}) {
    push @{$bioplex{$source}}, $target;
  } else {
    @{$bioplex{$source}} = ($target);
  }
  if (exists $bioplex{$target}) {
    push @{$bioplex{$target}}, $source;
  } else {
    @{$bioplex{$target}} = ($source);
  }
}
close $bioplexFH;

# Filter BioPlex preys that pass interaction cutoff
my @bioplexPreys;
foreach my $prey (keys %bioplex) {
  if (scalar @{$bioplex{$prey}} >= $minInteraction) {
    push @bioplexPreys, $prey;
  }
}
my $len = scalar @bioplexPreys;

# read SAINT file
print STDERR "Reading SAINT file\n";
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
my $saintTSV = Text::CSV_XS->new(\%tsvParams);
my @saintPreys;
$saintTSV->getline($saintFH); # discard header
while(my $row = $saintTSV->getline($saintFH)) {
  my $currFDR = @{$row}[15];
	if ($currFDR <= $fdr) {
    my $currPrey = @{$row}[2];
		push @saintPreys, $currPrey;
	}
}
close $saintFH;
@saintPreys = uniq @saintPreys;

my $bioplexPreyNum = scalar @bioplexPreys;
my $saintPreyNum = scalar @saintPreys;
my $intersection = intersect(@bioplexPreys, @saintPreys);
my $intersectionPercentage = sprintf "%.2f", ($intersection / $saintPreyNum) * 100;

print STDOUT "Filters: FDR - $fdr, min. interactions: $minInteraction\n";
print STDOUT "BioPlex interactors: $bioplexPreyNum\n";
print STDOUT "SAINT interactors: $saintPreyNum\n";
print STDOUT "Overlap: $intersection ($intersectionPercentage%)\n";
