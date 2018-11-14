#!/usr/bin/perl

# 23/6/2017

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# parameters
my $namespace = 'A';

# command line parameters
my $hfile = ''; # GO heirarchy file (obo format)

if ($#ARGV==0){
	print "\nTakes GO hierarchy in obo format and generates a two column list of\n";
	print "GO terms to IDs\n";
	print "\nusage:\n $0\n";
  print "-h [.obo (hierarchy) file]\n";
	print "-n [namespace to use, 'A' = all (default), 'C', 'F', 'P']\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-h'){
			$i++;
			$hfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-n'){
			$i++;
			$namespace = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# set namespace
my $outputFilehandle = '';
if ($namespace eq 'C') {
	$namespace = 'cellular_component';
	$outputFilehandle= '_cc';
} elsif ($namespace eq 'F') {
	$namespace = 'molecular_function';
	$outputFilehandle= '_mf';
} elsif ($namespace eq 'P') {
	$namespace = 'biological_process';
	$outputFilehandle= '_bp';
} else {
	$namespace = 'all';
}

# get direct parent terms for each GO term
open my $hierarchyFH, '<', $hfile or die "Could not open $hfile: $!";
open my $gotermsfh, '>', 'go-map' . $outputFilehandle . '.txt';
my $currId;
my $currTerm;
while(<$hierarchyFH>) {
	if ($_ =~ /^id: (\S+)/) {
		$currId = $1;
	} elsif ($_ =~ /^name: (.+)/) {
    $currTerm = $1;
	} elsif ($_ =~ /^namespace: (.+)/) {
		if (
			$namespace eq 'all' ||
			$1 eq $namespace
		) {
			print $gotermsfh "$currId\t$currTerm\n";
		}
	}
}
close $hierarchyFH;
close $gotermsfh;
