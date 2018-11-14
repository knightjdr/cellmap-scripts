#!/usr/bin/perl

# 23/6/2016

use strict;
use warnings;

# command line options
my $lfile = ''; 	#preys_cytoscape.txt

if ($#ARGV==0) {
	print "\nGets a unique list of genes from a correlation file\n";
	print "\nusage:\n $0\n";
	print "-l [preys_cytoscape.txt]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		}
		else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

my %geneList = ();

open(LFILE, "<$lfile") || die "$lfile can't be opened: $!";
while(<LFILE>){
	if($_ =~ /^(\S+)/) {
		if($1 ne "gene1") {
			$geneList{$1} = 1;
		}
	}
}
close(LFILE);

open(OFILE, ">gene_list.txt");
foreach my $gene (keys %geneList) {
	print OFILE "$gene\n";
}
close(OFILE);
