#!/usr/bin/perl

# 31/5/2017

use strict;
use warnings;

# libraries
use FindBin;
use lib "$FindBin::RealBin/../lib"; 

use List::MoreUtils qw(uniq);
use Prey::InteractorOutput qw(output);
use Text::CSV_XS;

# paramaters
my $requiredEvidence = 1;

# command line parameters
my $bfile = ''; # interaction file
my $cfile = '';	# correlation file
my $gfile = ''; # list of all genes to check for interactors

if ($#ARGV==0) {
	print "\nTakes a ProHits-viz correlation file and an interactions file, and calculates\n";
	print "the percentage of recovered known interactions for each correlation cutoff. Need to include a\n";
	print "file with a list of all genes to check.\n\n";
	print "\nusage:\n $0\n";
	print "-b [interaction file]\n";
	print "-c [correlation file]\n";
  print "-g [list of prey genes in correlation file]\n";
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

# parse interactor file
print STDERR "Getting list of interactors\n";
my %interactions;
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
while(my $row = $tsv->getline($fh)) {
  my $source = lc @{$row}[0];
  my $target = lc @{$row}[1];
  my @pair = ($source, $target);
  @pair = sort @pair;
  my $joinedPair = join '_', @pair;
  if (
    exists $geneFilterHash{$source} ||
    exists $geneFilterHash{$target} &&
    !(exists $interactions{$joinedPair})
  ) {
    $interactions{$joinedPair} = 1;
    $totalInteractions++;
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
    # next "if" is commented out to include self interaction
		# if ($source ne $target) {
		my @pair = ($source, $target);
	  @pair = sort @pair;
	  my $joinedPair = join '_', @pair;
	  if (exists $interactions{$joinedPair}) {
	    push @{$correlationInteractions[$correlation]{'recovered'}}, $joinedPair;
	  } else {
	    push @{$correlationInteractions[$correlation]{'discovered'}}, $joinedPair;
	  }
		# }
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
my $filename = 'prey-interactors-recovered.txt';
output($filename, \@correlationInteractions, $totalInteractions, 100);
