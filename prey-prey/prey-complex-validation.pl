#!/usr/bin/perl

# 21/7/2017

use strict;
use warnings;

# libraries
use FindBin;
use lib "$FindBin::RealBin/../lib";

use List::MoreUtils qw(uniq);
use Prey::InteractorOutput qw(output);
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # BioGRID file
my $cfile = '';	# correlation file
my $gfile = ''; # list of all genes to check for interactors
my $fileType = 'c'; # input file type

if ($#ARGV==0) {
	print "\nTakes a ProHits-viz correlation file and a list of genes in protein complexes, and calculates the percentage\n";
	print "of recovered known complex pairs for each correlation cutoff. Need to include a file with a list of all genes to\n";
	print "check.\n\n";
	print "\nusage:\n $0\n";
	print "-b [protein complex file]\n";
	print "-c [correlation file]\n";
  print "-g [list of prey genes in correlation file]\n";
	print "-t [input complex file type, c = CORUM or h = huMap]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-c'){
			$i++;
			$cfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
			$i++;
			$fileType = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# get list of prey genes
my @geneFilter;
my $geneTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
open my $geneFH, '<', $gfile or die "Could not open $gfile: $!";
while(my $row = $geneTSV->getline($geneFH)) {
	push @geneFilter, lc @{$row}[0];
}
close($geneFH);
my %geneFilterHash = map { $_ => 1 } @geneFilter;

# parse protein complex file
print STDERR "Getting list of complexes\n";
my %complexes;
open my $complexFH, '<', $bfile or die "Could not open $bfile: $!";
if ($fileType eq 'c') {
	my $tsv = Text::CSV_XS->new({
		sep_char => "\t",
	});
	$tsv->getline($complexFH); #discard header
	while(my $row = $tsv->getline($complexFH)) {
		my @currComplex = split ';', @{$row}[12];
		foreach my $gene (@currComplex) {
			$gene = lc($gene);
		}
		foreach my $gene (@currComplex) {
			if (exists $complexes{$gene}) {
				my $arrLength = scalar @{$complexes{$gene}};
				@{$complexes{$gene}[$arrLength]} = @currComplex;
			} else {
				@{$complexes{$gene}[0]} = @currComplex;
			}
		}
	}
} else {
	while(<$complexFH>) {
		chomp $_;
		my @currComplex = split '\t', $_;
		foreach my $gene (@currComplex) {
			$gene = lc($gene);
		}
		foreach my $gene (@currComplex) {
			if (exists $complexes{$gene}) {
				my $arrLength = scalar @{$complexes{$gene}};
				@{$complexes{$gene}[$arrLength]} = @currComplex;
			} else {
				@{$complexes{$gene}[0]} = @currComplex;
			}
		}
	}
}
close($complexFH);

# check correlation pairs
print STDERR "Reading correlation file\n\nProgress:\n";
my @correlationInteractions;
my $countLines = 0;
my $corrTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
open my $corrFH, '<', $cfile or die "Could not open $cfile: $!";
$corrTSV->getline($corrFH); #discard header
while(my $row = $corrTSV->getline($corrFH)) {
  my $correlation = int(@{$row}[2] * 100);
  if ($correlation >= 0) {
    my $source = lc @{$row}[0];
  	my $target = lc @{$row}[1];
		if ($source ne $target) {
			my $isComplex = 0;
			foreach my $complex (@{$complexes{$source}}) {
				my @currComplex = @{$complex};
				if (grep(/^$target$/, @currComplex)) {
					$isComplex = 1;
					last;
				}
			}
			my @pair = ($source, $target);
	    @pair = sort @pair;
	    my $joinedPair = join '_', @pair;
			if ($isComplex) {
				push @{$correlationInteractions[$correlation]{'recovered'}}, $joinedPair;
			} else {
				push @{$correlationInteractions[$correlation]{'discovered'}}, $joinedPair;
			}
		}
  }
  $countLines++;
  if ($countLines % 100000 == 0) {
    print STDERR ".";
  }
  if ($countLines % 1000000 == 0) {
    print STDERR "$countLines\n";
  }
}
close($corrFH);
print STDERR "\nRemoving duplicate interaction\n";
for(my $i = 0; $i <=100; $i++) {
  @{$correlationInteractions[$i]{'discovered'}} = uniq @{$correlationInteractions[$i]{'discovered'}};
  @{$correlationInteractions[$i]{'recovered'}} = uniq @{$correlationInteractions[$i]{'recovered'}};
}

# print fraction known
my $filename = 'prey-complex-recovered-corum.txt';
if ($fileType eq 'h') {
  $filename = 'prey-complex-recovered-humap.txt';
}
output($filename, \@correlationInteractions, 0, 100);
