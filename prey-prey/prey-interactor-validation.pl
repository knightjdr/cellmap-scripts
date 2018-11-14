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
my $cfile = '';	# correlation file
my $gfile = ''; # list of all genes to check for interactors

if ($#ARGV==0) {
	print "\nTakes a ProHits-viz correlation file and a BioGRID, IntAct or merged file, and calculates\n";
	print "the percentage of recovered known interactions for each correlation cutoff. Need to include a\n";
	print "file with a list of all genes to check.\n\n";
	print "\nusage:\n $0\n";
	print "-b [interaction file]\n";
	print "-c [correlation file]\n";
  print "-g [list of prey genes in correlation file]\n";
  print "-t [file type: b = BioGRID (default), i = IntAct, m = merged\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV) {
		if ($ARGV[$i] eq '-b') {
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-c') {
			$i++;
			$cfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-g') {
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t') {
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
my $tsv = Text::CSV_XS->new({ sep_char => "\t" });
open my $fh, '<', $gfile or die "Could not open $gfile: $!";
while(my $row = $tsv->getline($fh)) {
	push @geneFilter, lc @{$row}[0];
}
close($fh);
my %geneFilterHash = map { $_ => 1 } @geneFilter;

# parse BioGRID/IntAct file
print STDERR "Getting list of interactors\n";
my %biogrid;
my $totalInteractions = 0;
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
    my $sourceSpeciesCell = @{$row}[9];
  	my $targetSpeciesCell = @{$row}[10];
  	my ($sourceSpecies) = $sourceSpeciesCell =~ /^taxid:([-\d]+)\(/;
  	my ($targetSpecies) = $targetSpeciesCell =~ /^taxid:([-\d]+)\(/;
  	if ($sourceSpecies &&
  		$targetSpecies &&
  		(
  			($sourceSpecies == 9606 && $targetSpecies > 0) ||
  			($targetSpecies == 9606 && $sourceSpecies > 0)
  		)
  	) {
  		my $sourceCell = @{$row}[4];
  		my $targetCell = @{$row}[5];
  		my ($source) = $sourceCell =~ /uniprotkb:([^\(]+)\(gene name\)/;
  		my ($target) = $targetCell =~ /uniprotkb:([^\(]+)\(gene name\)/;
  		if($source && $target) {
        $source = lc $source;
    		$target = lc $target;
        my @pair = ($source, $target);
    		@pair = sort @pair;
    		my $joinedPair = join '_', @pair;
    		if (exists $geneFilterHash{$source} ||
          exists $geneFilterHash{$target} &&
          !(exists $biogrid{$joinedPair})
        ) {
    			$biogrid{$joinedPair} = 1;
          $totalInteractions++
    		}
      }
    }
  }
} elsif ($fileType eq 'm') {
	while(my $row = $tsv->getline($fh)) {
		my $source = lc @{$row}[0];
		my $target = lc @{$row}[1];
		my @pair = ($source, $target);
		@pair = sort @pair;
		my $joinedPair = join '_', @pair;
		if (
			exists $geneFilterHash{$source} ||
			exists $geneFilterHash{$target} &&
			!(exists $biogrid{$joinedPair})
		) {
			$biogrid{$joinedPair} = 1;
			$totalInteractions++;
		}
	}
} else {
  while(my $row = $tsv->getline($fh)) {
    my $source = lc @{$row}[7];
		my $target = lc @{$row}[8];
		my @pair = ($source, $target);
		@pair = sort @pair;
		my $joinedPair = join '_', @pair;
		if (exists $geneFilterHash{$source} ||
      exists $geneFilterHash{$target} &&
      !(exists $biogrid{$joinedPair})
    ) {
			$biogrid{$joinedPair} = 1;
      $totalInteractions++
		}
  }
}
close($fh);

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
		#if ($source ne $target) {
			my @pair = ($source, $target);
	    @pair = sort @pair;
	    my $joinedPair = join '_', @pair;
	    if (exists $biogrid{$joinedPair}) {
	      push @{$correlationInteractions[$correlation]{'recovered'}}, $joinedPair;
	  	} else {
	      push @{$correlationInteractions[$correlation]{'discovered'}}, $joinedPair;
	    }
		#}
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
print STDERR "\nRemoving duplicate interactions\n";
for(my $i = 0; $i <=100; $i++) {
  @{$correlationInteractions[$i]{'discovered'}} = uniq @{$correlationInteractions[$i]{'discovered'}};
  @{$correlationInteractions[$i]{'recovered'}} = uniq @{$correlationInteractions[$i]{'recovered'}};
}

# print fraction known
print STDERR "Formatting output\n";
open my $fractionfh, '>', 'prey-interactors-recovered.txt';
print $fractionfh "cutoff\tdiscovered\trecovered\tmissed\tfraction\n";
for(my $i = 0; $i <= 100; $i++) {
  my @discoveredArray;
  my @recoveredArray;
  for(my $j = $i; $j <= 100; $j++) {
    if ($correlationInteractions[$j]) {
      if (exists $correlationInteractions[$j]{'discovered'}) {
        push @discoveredArray, @{$correlationInteractions[$j]{'discovered'}};
      }
      if (exists $correlationInteractions[$j]{'recovered'}) {
        push @recoveredArray, @{$correlationInteractions[$j]{'recovered'}};
      }
    }
  }
  my $discovered = scalar uniq @discoveredArray;
	my $recovered = scalar uniq @recoveredArray;
  my $missed = $totalInteractions - $recovered;
	my $fraction = sprintf "%.3f", $recovered / ($discovered + $recovered);
	my $correlationCutoff = sprintf "%.2f", $i / 100;
	print $fractionfh "$correlationCutoff\t$discovered\t$recovered\t$missed\t$fraction\n"
}
close $fractionfh;
