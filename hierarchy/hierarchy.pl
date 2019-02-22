#!/usr/bin/perl

# 26/6/2017

use strict;
use warnings;

# libraries
use JSON::XS qw( );
use Text::CSV_XS;

# command line parameters
my $gfile = ''; # GO term-id map file
my $hfile = ''; # hierarchy in json format
my $nfile = ''; # rank details
my $sfile = ''; # domain details

if ($#ARGV==0) {
	print "\nAnnotates a hierarchy for use on the cell map.\n";
	print "\nusage:\n $0\n";
  print "-g [GO term-ID map]\n";
	print "-h [hierarchy in json format]\n";
	print "-n [NMF details file]\n";
	print "-s [SAFE details file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV) {
		if ($ARGV[$i] eq '-g') {
			$i++;
			$gfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-h') {
			$i++;
			$hfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-n') {
			$i++;
			$nfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-s') {
			$i++;
			$sfile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read in a map of GO terms to IDs if SAFE
print STDERR "Reading GO term-id map\n";
my %goMap;
open my $gomapFH, '<', $gfile or die "Could not open $gfile: $!";
my $gomapTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
while(my $row = $gomapTSV->getline($gomapFH)) {
	$goMap{@{$row}[1]} = @{$row}[0];
}
close $gomapFH;

# get NMF rank details
print STDERR "Reading NMF rank details\n";
my @usedNMFRanks;
my $nmfTSV = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	escape_char => undef,
	allow_loose_escapes => 1,
});
open my $nmfFH, '<', $nfile or die "Could not open $nfile: $!";
$nmfTSV->getline($nmfFH); # discard header
while(my $row = $nmfTSV->getline($nmfFH)) {
	my $termString = @{$row}[1];
	$termString =~ s/[\[\]"]//g;
	push @usedNMFRanks, split /, /, $termString;
}
close $nmfFH;

# get SAFE details
print STDERR "Reading SAFE domain details\n";
my @usedSAFERanks;
my $safeTSV = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	escape_char => undef,
	allow_loose_escapes => 1,
});
open my $safeFH, '<', $sfile or die "Could not open $sfile: $!";
$safeTSV->getline($safeFH); # discard header
while(my $row = $safeTSV->getline($safeFH)) {
	my $termString = @{$row}[1];
	$termString =~ s/[\[\]"]//g;
	push @usedSAFERanks, split /, /, $termString;
}
close $safeFH;

# read hierarchy
print STDERR "Reading and updating hierarchy\n";

sub annotateChildren {
	my @hierarchy = @{$_[0]};
	foreach my $child (@hierarchy) {
		if (exists $goMap{$child->{'name'}}) {
			$child->{'id'} = $goMap{$child->{'name'}};
		}
		if (grep(/^$child->{'name'}$/, @usedNMFRanks)) {
			$child->{'nmf'} = "true";
		} else {
			$child->{'nmf'} = "false";
		}
		if (grep(/^$child->{'name'}$/, @usedSAFERanks)) {
			$child->{'safe'} = "true";
		} else {
			$child->{'safe'} = "false";
		}
		if ($child->{'children'}) {
			$child->{'children'} = annotateChildren($child->{'children'});
		}
	}
	return \@hierarchy;
}

my $json_text = do {
  open my $json_fh, '<:encoding(UTF-8)', $hfile or die "Could not open $hfile: $!\n";
  local $/;
  <$json_fh>
};
my $json = JSON::XS->new;
my $hierarchy = annotateChildren($json->decode($json_text));
open my $outputFH, ">", "hierarchy_annotated.json";
print $outputFH $json->pretty(1)->encode($hierarchy);
close $outputFH;
