#!/usr/bin/perl

# 6/7/2017

use strict;
use warnings;

# libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # BioGRID file

if ($#ARGV==0){
  print "\nTakes a BioGRID file and reformats it so that edges can be displayed as known\n";
  print "Intractions in Cytoscape\n\n";
	print "\nusage:\n $0\n";
	print "-b [BioGRID file]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# parse BioGRID file
print STDERR "Parsing BioGrid\n";
my $biogridTSV = Text::CSV_XS->new({ sep_char => "\t" });
open my $biogridFH, '<', $bfile or die "Could not open $bfile: $!";
$biogridTSV->getline($biogridFH); #discard header
my @printed;
print STDOUT "name\tknown\n";
while(my $row = $biogridTSV->getline($biogridFH)) {
  my $source = @{$row}[7];
  my $target = @{$row}[8];
  my @pair = sort ($source, $target);
  my $mergedPair = $pair[0] . '*' . $pair[1];
  if (!grep(/^$mergedPair$/, @printed)) {
    push @printed, $mergedPair;
    print STDOUT "$source (interacts with) $target\t1\n";
    print STDOUT "$target (interacts with) $source\t1\n";
  }
}
close $biogridFH;
