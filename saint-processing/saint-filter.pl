#!/usr/bin/perl

# 21/6/2017

use strict;
use warnings;

# libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;
use String::Util qw(trim);
use Text::CSV_XS;

# paramaters
my $minimumPreys = 5;
my $score = 0.01; # FDR

# command line parameters
my $cfile = ''; # contaminants file
my $rfile = ''; # optional list of baits to remove
my $sfile = '';	# SAINT file

if($#ARGV==0){
	print "\nTakes a SAINT file and filters out decoys (from list),\n";
  print "tags, etc, and removes any baits with less preys that a 1% FDR cutoff (default 5)\n\n";
	print "\nusage:\n $0\n";
  print "-c [contaminants file]\n";
	print "-f [score filter]\n";
	print "-m [minimum number of preys]\n";
	print "-r [optional list of baits to remove]\n";
	print "-s [SAINT file]\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
    if ($ARGV[$i] eq '-c'){
			$i++;
			$cfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-f'){
			$i++;
			$score = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-m'){
			$i++;
			$minimumPreys = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-r'){
			$i++;
			$rfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read contaminants file
print STDERR "Reading contaminants file\n";
my @contaminants;
my $tsv = Text::CSV_XS->new({
  sep_char => "\t"
});
open my $contaminantsFH, '<', $cfile or die "Could not open $cfile: $!";
while(my $row = $tsv->getline($contaminantsFH)) {
  push @contaminants, @{$row}[0];
}
close $contaminantsFH;

# check for optional baits to remove
my @excludeBaits;
if ($rfile) {
	my $excludeTSV = Text::CSV_XS->new({
	  sep_char => "\t"
	});
	open my $excludeFH, '<', $rfile or die "Could not open $rfile: $!";
	while(my $row = $tsv->getline($excludeFH)) {
	  push @excludeBaits, @{$row}[0];
	}
	close $excludeFH;
}

# read SAINT file
print STDERR "Reading SAINT file\n";
my %baitsAboveCutoff;
my %preys;
$tsv = Text::CSV_XS->new({
  sep_char => "\t"
});
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
my $header = join "\t", @{$tsv->getline($saintFH)};
while(my $row = $tsv->getline($saintFH)) {
	if (!(exists $baitsAboveCutoff{@{$row}[0]})) {
		$baitsAboveCutoff{@{$row}[0]} = 0;
	}
	my $trimmedPrey = trim @{$row}[1];
  if (!(grep(/^$trimmedPrey$/, @contaminants)) &&
    $trimmedPrey !~ /^DECOY/ &&
		!grep(/^@{$row}[0]$/, @excludeBaits)
  ) {
    if (@{$row}[15] <= $score) {
      $baitsAboveCutoff{@{$row}[0]}++;
    }
    push @{$preys{@{$row}[0]}}, join "\t", @{$row};
  }
}
close $saintFH;

# output new file
print STDERR "Outputting file\n";
my ($newFileHandle) = $sfile =~ /(.*)\./;
open my $filteredFH, '>', 'filtered_baits.txt';
open my $outputFH, '>', $newFileHandle . '_filtered.txt';
print $outputFH "$header\n";
foreach my $bait (sort keys %baitsAboveCutoff) {
  if ($baitsAboveCutoff{$bait} >= $minimumPreys) {
    for my $prey (@{$preys{$bait}}) {
      print $outputFH "$prey\n";
    }
  } else {
		print $filteredFH "$bait, preys: $baitsAboveCutoff{$bait}\n";
	}
}
close $filteredFH;
close $outputFH;
