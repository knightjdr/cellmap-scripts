#!/usr/bin/perl

# 18/6/2018

use strict;
use warnings;

# libraries
use List::Util qw(max);
use Statistics::Basic qw(mean);
use Text::CSV_XS;

# Command line parameters.
my $sfile = ''; # SAINT file.

if ($#ARGV==0){
	print "\nTakes a SAINT file and for each prey calculates its average spectral\n";
  print "count and the best FDR it was seen with. The spectral count is calculated\n";
  print "as is or as control subtracted.\n";
	print "\nusage:\n $0\n";
	print "-s [SAINT file]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Vars.
my %baits; # Use for storing unique baits.
my %preys;

# Read SAINT file and sum spectral counts and best FDR.
print STDERR "Reading SAINT file\n";
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!";
my $tsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
$tsv->getline($saintFH); # Discard header.
while(my $row = $tsv->getline($saintFH)) {
  my $avgspec = @{$row}[5];
  my $bait = @{$row}[0];
  my $ctrlavg = mean split '\|', @{$row}[7];
  my $fdr = @{$row}[15];
  my $prey = @{$row}[2];
  if (exists $preys{$prey}) {
    if ($fdr < $preys{$prey}{'fdr'}) {
      $preys{$prey}{'fdr'} = $fdr;
    }
    $preys{$prey}{'spec'} += $avgspec;
    $preys{$prey}{'specContSub'} += max($avgspec - $ctrlavg, 0);
  } else {
    $preys{$prey}{'fdr'} = $fdr;
    $preys{$prey}{'spec'} = $avgspec;
    $preys{$prey}{'specContSub'} = max($avgspec - $ctrlavg, 0);
  }
  if (!(exists $baits{$bait})) {
    $baits{$bait} = 1;
  }
}
close $saintFH;

# Calculate averages.
open my $averagesfh, '>', 'prey-averages.txt';
print $averagesfh "prey\tavgSpec\tcontrolSubAvgSpec\tfdr\n";
my $noBaits = scalar(keys %baits);
foreach my $prey (sort {lc $a cmp lc $b} keys %preys) {
  $preys{$prey}{'spec'} = sprintf "%.2f", $preys{$prey}{'spec'} / $noBaits;
  $preys{$prey}{'specContSub'} = sprintf "%.2f", $preys{$prey}{'specContSub'} / $noBaits;
  print $averagesfh "$prey\t$preys{$prey}{'spec'}\t$preys{$prey}{'specContSub'}\t$preys{$prey}{'fdr'}\n";
}
close $averagesfh;
