#!/usr/bin/perl

# 23/6/2016

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# command line
my $gfile = ''; 	#go annotation file
my $ofile = ''; 	#go heirarchy file (obo format)
my $lfile = ''; 	#gene list
my $namespace = 'cc';

if ($#ARGV==0) {
	print "\nTakes GO association file, a GO hierarchy in obo format, a list of genes and the GO namespace to use and outputs\n";
	print "a tab delimited table of Genes versus terms\n";
	print "\nusage:\n $0\n";
	print "-g [.go annotation file]\n";
	print "-f [.obo file]\n";
	print "-l [gene list]\n";
	print "-n [namespace (must be one of 'all', 'bp', 'cc' or 'mf']\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		}
		elsif($ARGV[$i] eq '-f'){
			$i++;
			$ofile = $ARGV[$i];
		}
		elsif($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		}
		elsif($ARGV[$i] eq '-n'){
			$i++;
			$namespace = $ARGV[$i];
		}
		else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

my $shortterm;
my $longNameSpace;
if($namespace eq "bp") {
	$shortterm = "P";
	$longNameSpace = "biological_process";
}
elsif($namespace eq "cc") {
	$shortterm = "C";
	$longNameSpace = "cellular_component";
}
elsif($namespace eq "mf") {
	$shortterm = "F";
	$longNameSpace = "molecular_function";
}
else {
	$longNameSpace = "all";
}

#retrieve child parent relations for obo file
my %alt_id;
my $currId;
my $currName;
my $inNamespace;
my $obsolete;
my %idMap = ();
my %hierarchy = ();
open(OFILE, "<$ofile") || die "$ofile can't be opened: $!";
while(<OFILE>) {
	if ($_ =~ /^\[Term\]/) {
		$currId = '';
		$currName = '';
		$inNamespace = 0;
		$obsolete = 0;
	} elsif ($_ =~ /^id: (\S+)/) {
		$currId = $1;
	} elsif ($_ =~ /^alt_id: (\S+)/) {
		$alt_id{$1} = $currId;
	} elsif ($_ =~ /^name: (.+)/) {
		$idMap{$currId} = $1;
	} elsif ($_ =~ /^namespace: (\S+)/) {
		if($longNameSpace eq "all" || $longNameSpace eq $1) {
			$inNamespace = 1;
		}
	} elsif ($_ =~ /^is_obsolete: true/) {
		$obsolete = 1;
	} elsif ($_ =~ /^is_a: (\S+) ! .+/) {
		my $parent = $1;
		if($inNamespace && !$obsolete) {
			if(exists($hierarchy{$currId})) {
				push(@{$hierarchy{$currId}}, $parent);
			}
			else {
				$hierarchy{$currId} = [$parent];
			}
		}
	} elsif ($_ =~ /^relationship: part_of (\S+) ! .+/) {
		my $parent = $1;
		if($inNamespace && !$obsolete) {
			if(exists($hierarchy{$currId})) {
				push(@{$hierarchy{$currId}}, $parent);
			}
			else {
				$hierarchy{$currId} = [$parent];
			}
		}
	}
	else {
	}
}
close(OFILE);

#get gene list
my @geneList;
open(LFILE, "<$lfile") || die "$lfile can't be opened: $!";
while(<LFILE>){
	if($_ =~ /^(\S+)/) {
		push(@geneList, $1);
	}
}
my %geneListHash = map { $_ => 1 } @geneList;

# get Gene terms
my %genes;
my %terms;
open my $annotationsFH, '<', $gfile or die "Could not open $gfile: $!";
my $annotationsTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
while(my $row = $annotationsTSV->getline($annotationsFH)) {
	my $currGene = @{$row}[2];
	my $currTerm = @{$row}[4];
	if(
		@{$row}[8] eq $shortterm &&
		exists($geneListHash{$currGene})
	) {
		# change alt_ids to normal
		if (exists $alt_id{$currTerm}) {
			$currTerm = $alt_id{$currTerm};
		}
		# in case alt_ids and regular ids are both being listed for a gene
		if (!(grep(/^$currTerm$/, @{$genes{$currGene}}))) {
			push @{$genes{$currGene}}, $currTerm;
		}
		$terms{$currTerm} = 1;
	}
}
close $annotationsFH;

#add parent terms to genes GO terms
foreach my $gene (keys %genes) {
	for(my $i = 0; $i < scalar(@{$genes{$gene}}); $i++) {
		my $currTerm = $genes{$gene}[$i];
		if(exists($hierarchy{$currTerm})) {
			for(my $j = 0; $j < scalar(@{$hierarchy{$currTerm}}); $j++) {
				my $parentTerm = $hierarchy{$currTerm}[$j];
				if(!exists $terms{$parentTerm}) {
					$terms{$parentTerm} = 1;
				}
				my %geneTerms = map { $_ => 1 } @{$genes{$gene}};
				if(!exists($geneTerms{$parentTerm})) {
					push(@{$genes{$gene}}, $parentTerm);
				}
			}
		}
	}
}

#map all terms to genes
my %genesBinary = ();
my @terms = keys %terms;
my $numTerms = scalar(@terms);
foreach my $gene (keys %genes) {
	my %currTerms = map { $_ => 1 } @{$genes{$gene}};
	for(my $i = 0; $i < $numTerms; $i++) {
		if (exists($currTerms{$terms[$i]})) {
			$genesBinary{$gene}{$terms[$i]} = 1;
		} else {
			$genesBinary{$gene}{$terms[$i]} = 0;
		}
	}
}

#print table
my $outfileName = $namespace . '_annotations.txt';
open(OFILE, ">$outfileName");
print OFILE "\t";
for(my $i = 0; $i < $numTerms; $i++) {
	if(exists($idMap{$terms[$i]})) {
		print OFILE "$idMap{$terms[$i]}";
	}
	else {
		print OFILE "$terms[$i]";
	}
	if($i < $numTerms - 1) {
		print OFILE "\t";
	}
}
print OFILE "\n";
foreach my $gene (keys %genes) {
	print OFILE "$gene\t";
	for(my $i = 0; $i < $numTerms; $i++) {
		print OFILE "$genesBinary{$gene}{$terms[$i]}";
		if($i < $numTerms - 1) {
			print OFILE "\t";
		}
	}
	print OFILE "\n";
}
close(OFILE);
