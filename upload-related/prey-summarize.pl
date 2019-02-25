#!/usr/bin/perl

# 19/6/2018

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use List::Util qw(max min);
use POSIX qw(floor);
use Statistics::Basic qw(mean);
use Text::CSV_XS;

# parameters
my $fdr = 0.01;

# command line parameters
my $bfile = ''; # basis matrix
my $sfile = ''; # SAINT file.

if ($#ARGV==0){
	print "\nTakes a SAINT file and an NMF basis matrix and for each significant prey counts the number";
  print "of baits it was seen with, the average spectral count and the NMF ratio for its.";
  print "best localization relative to the max value for that rank\n";
	print "\nusage:\n $0\n";
  print "-b [basis.csv]\n";
	print "-s [SAINT file]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
      $i++;
      $bfile = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Read basis matrix. Grab each prey's top value and determine the maximum
# value for each rank.
print STDERR "Reading basis matrix\n";
open my $basisFH, '<', $bfile or die "Could not open $bfile: $!";
my $basisTsv = Text::CSV_XS->new({ sep_char => "," });
$basisTsv -> getline($basisFH); # Discard header
my @maxRank;
my %preyTopRank;
while(my $row = $basisTsv->getline($basisFH)) {
  my $gene = @{$row}[0];
  shift @{$row};
  my @rankValues = @{$row};
  $preyTopRank{$gene}{'rank'} = 0;
  $preyTopRank{$gene}{'value'} = 0;
  for (my $i = 0; $i < scalar @rankValues; $i++) {
    if ($rankValues[$i] > $maxRank[$i]) {
      $maxRank[$i] = $rankValues[$i];
    }
    if ($rankValues[$i] > $preyTopRank{$gene}{'value'}) {
      $preyTopRank{$gene}{'rank'} = $i;
      $preyTopRank{$gene}{'value'} = $rankValues[$i];
    }
  }
}
close $basisFH;

# Read SAINT file and for each significant prey, count number of baits store avgspec.
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
my $saintTsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
$saintTsv->getline($saintFH); # discard header
my %preys;
while(my $row = $saintTsv->getline($saintFH)) {
	my $avgspec = @{$row}[5];
	my $prey = @{$row}[2];
	if (@{$row}[15] <= $fdr) {
		if (exists $preys{$prey}) {
			$preys{$prey}{'baits'}++;
			push @{$preys{$prey}{'avgspec'}}, $avgspec;
		} else {
			$preys{$prey}{'baits'} = 1;
			@{$preys{$prey}{'avgspec'}} = [$avgspec];
		}
	}
}
close $saintFH;

# Calculate average spectral counts for each prey across the baits it was found with.
# Average of the AvgSpec.
foreach my $prey (keys %preys) {
	$preys{$prey}{'avg'} = mean @{$preys{$prey}{'avgspec'}};
}

# Output prey numbers with color representing level.
open my $outputFH, '>', 'prey-summary.txt';
print $outputFH "prey\tbaits\tavgspec\tratio\n";
foreach my $prey (sort {lc $a cmp lc $b} keys %preys) {
  my $avg = sprintf "%.3f", $preys{$prey}{'avg'};
  my $ratio = sprintf "%.3f", $preyTopRank{$prey}{'value'} / $maxRank[$preyTopRank{$prey}{'rank'}];
  print $outputFH "$prey\t$preys{$prey}{'baits'}\t$avg\t$ratio\n";
}
close $outputFH;
