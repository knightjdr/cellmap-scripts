#!/usr/bin/perl

# 10/7/2017

use strict;
use warnings;

# libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;
use IO::Handle;
use List::MoreUtils qw(uniq);
use String::Util qw(trim);
use Text::CSV_XS;

# redirect STDERR to file
#open ERROR,  '>', "error.txt"  or die $!;
#STDERR->fdopen( \*ERROR,  'w' ) or die $!;

# parameters
my $namespace = 'C';

# command line parameters
my $ifile = ''; # file with information content
my $gfile = ''; # goa_human file without header lines
my $mfile = ''; # map of bait names to desired gene names to retrieve information for (or use for analysis)
my $ufile = ''; # Uniprot file for getting gene synonyms.

if ($#ARGV==0) {
	print "\nFor a list of bait names, output synonyms and the best GO terms (based\n";
	print "on IC content score).\n";
	print "\nusage:\n $0\n";
	print "-i [IC content files]\n";
  print "-g [GO localization file (goa_human.gaf)]\n";
	print "-m [map of bait names to gene names]\n";
	print "-n [GO namespace, C (default), F or P]\n";
	print "-u [Uniprot database]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-i'){
			$i++;
			$ifile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-m'){
			$i++;
			$mfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-u'){
			$i++;
			$ufile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-n'){
			$i++;
			$namespace = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

sub removeValue {
	my $str = $_[0];
	my @arr = @{$_[1]};
	my $index = 0;
	if (scalar @arr) {
		$index++ until $arr[$index] eq $str;
		splice(@arr, $index, 1);
	}
	return @arr;
}

# get bait gene names
print STDERR "Reading bait names\n";
open my $baitFH, '<', $mfile or die "Could not open $mfile: $!";
my $baitTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
my %baits;
my @genes;
$baitTSV->getline($baitFH); # discard header
while(my $row = $baitTSV->getline($baitFH)) {
	$baits{@{$row}[0]} = @{$row}[1];
	push @genes, @{$row}[1];
}
close $baitFH;
@genes = uniq @genes;

# get IC
print STDERR "Reading information content\n";
open my $icFH, '<', $ifile or die "Could not open $ifile: $!";
my $icTSV = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
my %ic;
$icTSV->getline($icFH); # discard header
while(my $row = $icTSV->getline($icFH)) {
	$ic{@{$row}[0]}{'name'} = @{$row}[2];
	$ic{@{$row}[0]}{'score'} = @{$row}[5];
}
close $icFH;

# get GO terms
print STDERR "Retrieving GO terms\n";
open my $goFH, '<', $gfile or die "Could not open $gfile: $!";
my $goTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
my @goMatchedGenes;
my %termsPerGene;
while(my $row = $goTSV->getline($goFH)) {
	if (grep(/^@{$row}[2]$/, @genes)) {
		push @goMatchedGenes, @{$row}[2];
		if (@{$row}[8] eq $namespace) {
			push @{$termsPerGene{@{$row}[2]}}, @{$row}[4];
		}
	}
}
close $goFH;
@goMatchedGenes = uniq @goMatchedGenes;

# get gene synonyms
print STDERR "Parsing UniProt data\n";
my $currGene;
my $currLocString = '';
my @currNames;
my %geneLocalizations;
my %geneSynonyms;
my $gettingLoc = 0;
my $reviewed;
my $species;
open(UFILE, "<$ufile") || die "$ufile can't be opened: $!";
while(<UFILE>){
	if ($_ =~ /^ID\s+\S+\s+Reviewed;/) {
		$reviewed = 1;
	} elsif ($_ =~ /^OS   ([^(]+)/) {
		$species = trim $1;
	} elsif ($_ =~ /^GN/) {
		if ($_ =~ /Name=([^;{]+)/) {
			$currGene = trim $1;
			push @currNames, $currGene;
		}
		if ($_ =~ /Synonyms=([^;{]+)/) {
			push @currNames, split ', ', trim $1;
		}
		if ($_ =~ /OrderedLocusNames=([^;{]+)/) {
			push @currNames, split ', ', trim $1;
		}
		if ($_ =~ /ORFNames=([^;{]+)/) {
			push @currNames, split ', ', trim $1;
		}
	} elsif ($_ =~ /^CC   -!- SUBCELLULAR LOCATION:(.*)/) {
		$currLocString = $1;
		$gettingLoc = 1
	} elsif ($_ =~ /^CC   -!-/) {
		$gettingLoc = 0;
	} elsif (
		$_ =~ /^CC\s+(.+)/ &&
		$gettingLoc
	) {
		$currLocString .= ' ' . $1;
	} elsif ($_ =~ /^\/\//) {
		if (
			$reviewed &&
			$species eq 'Homo sapiens'
		) {
			if (scalar @currNames) {
				@currNames = uniq @currNames;
				@currNames = sort @currNames;
			}
			# create localization array
			my @currLoc = split '\.', $currLocString;
			if (scalar @currLoc) {
				for(my $i = 0 ; $i < scalar @currLoc; $i++) {
					$currLoc[$i] = trim $currLoc[$i];
					$currLoc[$i] =~ s/\s?{.*};?//g;
					$currLoc[$i] =~ s/[;,]/ -/g;
					if (!$currLoc[$i]) {
						splice @currLoc, $i, 1;
						$i--;
					} elsif ($currLoc[$i] =~ /^Note/) {
						splice @currLoc, $i, scalar @currLoc - $i;
					}
				}
				@currLoc = sort { lc($a) cmp lc($b) } @currLoc;
			}
			# for each offical gene name, assign values to it
			if ($currGene) {
				@{$geneLocalizations{$currGene}} = @currLoc;
				@{$geneSynonyms{$currGene}} = removeValue($currGene, \@currNames);
			}
			# for none official gene names, check if they exist first
			foreach my $gene (@currNames) {
				if(!(exists $geneSynonyms{$gene})) {
					@{$geneLocalizations{$gene}} = @currLoc;
					@{$geneSynonyms{$gene}} = removeValue($gene, \@currNames);
				}
			}
		}
		#reset some params
		$currGene = '';
		$currLocString = '';
		@currNames = ();
		$gettingLoc = 0;
		$reviewed = 0;
		$species = '';
	} else {
		$gettingLoc = 0;
	}
}
close(UFILE);

# output
print STDERR "Outputting information\n";
open my $outputFH, '>', 'bait_details.txt';
print $outputFH "bait name\tsymbol\taliases\tspecific GO term\tUniProt localization\n";
foreach my $bait (sort {lc $a cmp lc $b} keys %baits) {
	my $bestIC = 0;
	my $currGene = $baits{$bait};
	my @currTerms;
	if (exists $termsPerGene{$currGene}) {
		for(my $i = 0; $i < scalar @{$termsPerGene{$currGene}}; $i++) {
			my $currID = $termsPerGene{$currGene}[$i];
			if (exists $ic{$currID}) {
				if ($ic{$currID}{'score'} > $bestIC) {
					$bestIC = $ic{$currID}{'score'};
					@currTerms = ($ic{$currID}{'name'} . ' (' . $currID . ')');
				} elsif($ic{$currID}{'score'} == $bestIC) {
					if (scalar @currTerms > 0) {
						push @currTerms, $ic{$currID}{'name'} . ' (' . $currID . ')'
					} else {
						@currTerms = ($ic{$currID}{'name'} . ' (' . $currID . ')');
					}
				}
			}

		}
	}
	my $goString;
	if (scalar @currTerms) {
		@currTerms = sort { lc($a) cmp lc($b) } uniq @currTerms;
		$goString = join ', ', @currTerms;
	} else {
		if (grep(/^$currGene$/, @goMatchedGenes)) {
			$goString = '-'
		} else {
			$goString = 'WARNING: Gene not in GO';
		}
	}
	my $uniprotString;
	if (exists $geneLocalizations{$currGene}) {
		$uniprotString = join ', ', @{$geneLocalizations{$currGene}};
	} else {
		$uniprotString = 'WARNING: Gene not in UniProt';
	}
	my $aliasString;
	if (exists $geneSynonyms{$currGene}) {
		$aliasString = join ', ', @{$geneSynonyms{$currGene}};
	} else {
		$aliasString = 'WARNING: Gene not in UniProt';
	}
	print $outputFH "$bait\t$currGene\t$aliasString\t$goString\t$uniprotString\n";
}
close $outputFH;
