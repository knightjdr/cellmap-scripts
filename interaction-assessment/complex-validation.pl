#!/usr/bin/perl

# 31/5/2017

use strict;
use warnings;

# libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# paramaters
my $fileType = 'b';
my $requiredEvidence = 1;

# command line parameters
my $bfile = ''; # BioGRID file
my $ffile = ''; # list of genes to include in the analysis
my $mfile = '';	# file with a map of Bait names to gene names "bait.dat" from ProHITS
my $sfile = '';	# SAINT file

if($#ARGV==0){
	print "\nTakes a SAINT file and a BioGRID/IntAct file, and calculates the percentage of recovered\n";
	print "known interactions for each FDR probability. Also requires a file for mapping\n";
	print "Bait names to gene names. Also, optionally, a list of genes can be include to filter\n";
	print "the bait list by\n\n";
	print "\nusage:\n $0\n";
	print "-b [BioGRID/IntAct file]\n";
	print "-f [List of genes to include]\n";
	print "-m [SAINT bait-gene map]\n";
	print "-s [SAINT file]\n";
	print "-t [file type: b = BioGRID (default), i = IntAct, m = merged\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-f'){
			$i++;
			$ffile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-m'){
			$i++;
			$mfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
			$i++;
			$fileType = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read bait name to gene map
print STDERR "Creating map\n";
my %baitMap;
my $tsv = Text::CSV_XS->new({ sep_char => "\t" });
open my $fh, '<', $mfile or die "Could not open $mfile: $!";
$tsv->getline($fh); #discard header
while(my $row = $tsv->getline($fh)) {
	my $bait = @{$row}[0];
	my $gene = @{$row}[1];
	$baitMap{$bait} = $gene;
	if (!$gene) {
		print STDERR "Missing gene name in bait map file, for bait $bait\n";
		die;
	}
}
close($fh);

# get list of bait genes (from filter file if specified, otherwise all baits)
my @baitFilter;
if ($ffile) {
	$tsv = Text::CSV_XS->new({ sep_char => "\t" });
	open $fh, '<', $ffile or die "Could not open $ffile: $!";
	while(my $row = $tsv->getline($fh)) {
		push @baitFilter, lc @{$row}[0];
	}
	close($fh);
} else {
	@baitFilter = values %baitMap;
	@baitFilter = uniq @baitFilter;
	for my $v (@baitFilter) {
		$v = lc($v);
	}
}
my %baitFilterHash = map { $_ => 1 } @baitFilter;

# parse Interactions file
print STDERR "Retrieving interactors\n";
my %biogrid;
$tsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
open $fh, '<', $bfile or die "Could not open $bfile: $!";
$tsv->getline($fh); #discard header
if ($fileType eq 'i') {
	while(my $row = $tsv->getline($fh)) {
		my $source = lc @{$row}[0];
		my $target = lc @{$row}[1];
		my @pair = ($source, $target);
		@pair = sort @pair;
		my $joinedPair = join '_', @pair;
		if (exists $baitFilterHash{$source} || exists $baitFilterHash{$target}) {
			my $approach = @{$row}[2];
			push @{$biogrid{$joinedPair}}, $approach;
		}
	}
} elsif ($fileType eq 'm') {
	while(my $row = $tsv->getline($fh)) {
		my $source = lc @{$row}[0];
		my $target = lc @{$row}[1];
		my @pair = ($source, $target);
		@pair = sort @pair;
		my $joinedPair = join '_', @pair;
		if (exists $baitFilterHash{$source} || exists $baitFilterHash{$target}) {
			@{$biogrid{$joinedPair}} = split /;/, @{$row}[2];
		}
	}
} else {
	while(my $row = $tsv->getline($fh)) {
		my $approach = @{$row}[11];
		my $source = lc @{$row}[7];
		my $target = lc @{$row}[8];
		my @pair = ($source, $target);
		@pair = sort @pair;
		my $joinedPair = join '_', @pair;

		# if the evidence types need not be unique, uncomment the next line
		if (exists $baitFilterHash{$source} || exists $baitFilterHash{$target}) {
			push @{$biogrid{$joinedPair}}, $approach;
		}

		# if the evidence types must be unique, uncomment next if/else block
		#if (exists $biogrid{$joinedPair}) {
		#	if ( !(grep /^$approach$/, @{$biogrid{$joinedPair}}) ) {
		#		push @{$biogrid{$joinedPair}}, $approach;
		#	}
		#} else {
		#	push @{$biogrid{$joinedPair}}, $approach;
		#}
	}
}
close($fh);

# remove all Interactions pairs that do not have at least $requiredEvidence pieces of evidence
print STDERR "Removing evidence-type pairs\n";
my @biogridPairs;
while (my ($key, $value) = each %biogrid) {
	if (scalar @{$value} >= $requiredEvidence) {
		push @biogridPairs, $key
	}
}
my %biogridPairsHash = map { $_ => 1 } @biogridPairs;

# check SAINT pairs
print STDERR "Reading SAINT file\n";
my @saintInteractions;
for(my $i = 0; $i <=100; $i++) {
	$saintInteractions[$i]{'known'} = 0;
	$saintInteractions[$i]{'total'} = 0;
}
$tsv = Text::CSV_XS->new({ sep_char => "\t" });
open $fh, '<', $sfile or die "Could not open $sfile: $!";
$tsv->getline($fh); #discard header
while(my $row = $tsv->getline($fh)) {
	if (!(exists $baitMap{@{$row}[0]})) {
		print STDERR "This bait is not in the map: @{$row}[0]\n";
		die;
	}
	my $bait = lc $baitMap{@{$row}[0]};
	if (exists $baitFilterHash{$bait}) {
		my $prey = lc @{$row}[2];
		my $fdr = @{$row}[15] * 100;
		my @pair = ($bait, $prey);
		@pair = sort @pair;
		my $joinedPair = join '_', @pair;
		if (exists $biogridPairsHash{$joinedPair}) {
			$saintInteractions[$fdr]{'known'}++;
		}
		$saintInteractions[$fdr]{'total'}++;
	}
}
close($fh);

# print fraction known
open my $fractionfh, '>', 'fraction-recovered.txt';
print $fractionfh "cutoff (FDR)\tfraction\tknown\ttotal\n";
for(my $i = 0; $i <= 100; $i++) {
	my $known = 0;
	my $total = 0;
	for(my $j = 0; $j <= $i; $j++) {
		$known += $saintInteractions[$j]{'known'};
		$total += $saintInteractions[$j]{'total'};
	}
	my $saintCutoff = $i / 100;
	if($total != 0) {
		my $fraction =  sprintf "%.5f", $known / $total;
		print $fractionfh "$saintCutoff\t$fraction\t$known\t$total\n"
	} else {
		print $fractionfh "$saintCutoff\tNA\t0\t0\n"
	}
}
close $fractionfh;
