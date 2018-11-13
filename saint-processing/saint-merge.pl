#!/usr/bin/perl

# 31/5/2018
# This script will merge a group of SAINT files. It takes a list of bait names
# to keep and will only keep those from the SAINT files. It must be run in a folder
# with the SAINT files and it will run on any files with the .txt extension. It will
# also only keep the first 15 columns from the files.

use strict;
use warnings;

# Packages
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # file with Baits. Should have a header and the first column has the bait names
my $noColumns = 22;

if ($#ARGV==0) {
	print "\nRun in a folder with SAINT files (and only SAINT files) and will merge\n";
  print "the files keeping the first 16 columns. It will only keep baits in the\n";
  print "supplied bait list file\n";
	print "\nusage:\n $0\n";
	print "-b [bait list]\n";
  print "-n [number of columns to keep (default 22)\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-n'){
			$i++;
			$noColumns = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Read bait list.
print STDERR "Reading bait list\n";
my %baitList;
my $tsv = Text::CSV_XS->new({ sep_char => "\t" });
open my $fh, '<', $bfile or die "Could not open $bfile: $!";
$tsv->getline($fh); #discard header
while(my $row = $tsv->getline($fh)) {
	my $bait = lc @{$row}[0];
	$baitList{$bait} = 1;
}
close($fh);

# Header columns
my @header = ('Bait', 'Prey', 'PreyGene', 'Spec', 'SpecSum', 'AvgSpec', 'NumReplicates', 'ctrlCounts', 'AvgP', 'MaxP', 'TopoAvgP', 'TopoMaxP', 'SaintScore', 'logOddsScore', 'FoldChange', 'BFDR', 'boosted_by', 'UniqueSpec', 'UniqueSpecSum', 'UniqueAvgSpec', 'PreySequenceLength', 'UniProtID');
my $headerString = join "\t", splice(@header, 0, $noColumns);

# Open merge file for writing and add header.
open my $merge, '>', 'merged.txt';
print $merge "$headerString\n";

# keep track of baits found so I can output baits missing.
my @baitsFound;
my %baitsSkipped;

# Read SAINT files.
my $dir = '.';
foreach my $file (glob("$dir/*.txt")) {
  next if $file eq './merged.txt';
  printf "%s\n", $file;
  my $tsv = Text::CSV_XS->new({ sep_char => "\t" });
  open my $fh, '<', $file or die "Could not open $file: $!";
  my $header = $tsv->getline($fh);

  # Not all SAINT files have the logsOddsScore column. Check if they do. If not,
  # add 0 in place.
  my $hasLogOdds = 0;
  if (@{$header}[13] eq 'logOddsScore') {
    $hasLogOdds = 1;
  }
  while(my $row = $tsv->getline($fh)) {
  	my $bait = lc @{$row}[0];
  	if (exists $baitList{$bait}) {
      push @baitsFound, $bait;
      my @rowSubset;
      if ($hasLogOdds == 1) {
        @rowSubset = splice @{$row}, 0, $noColumns;
      } else {
        my @end = splice @{$row}, 13, $noColumns - 13 - 1;
        @rowSubset = splice @{$row}, 0, 13;
        push @rowSubset, 0;
        push @rowSubset, @end;
      }
      my $row = join "\t", @rowSubset;
      print $merge "$row\n";
    } elsif (!exists $baitsSkipped{$bait}) {
      $baitsSkipped{$bait} = 1;
    }
  }
  close $fh;
}

close $merge;

# Output baits missing.
my %found = map { $_ => 1 } @baitsFound;
open my $missing, '>', 'missing_merged.txt';
foreach my $bait (sort keys %baitList) {
  if (!exists $found{$bait}) {
    print $missing "$bait\n";
  }
}
close $missing;

# Output baits skipped.
open my $skipped, '>', 'skipped_merged.txt';
foreach my $bait (sort keys %baitsSkipped) {
  print $skipped "$bait\n";
}
close $skipped;
