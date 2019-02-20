#!/usr/bin/perl

# 4/7/2017

use strict;
use warnings;

# libraries
use Algorithm::MedianSelect qw(median);
use Array::Utils qw(:all);
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(indexes);
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# parameters
my $dirAttributes = './attributes';
my $dirNodes = './nodes';
my $namespace = 'C';

# command line parameters
my $gfile = ''; # file with all child terms for each GO term
my $lfile = ''; # goa_human.gaf
my $mfile = ''; # map of GO terms to ID's

if ($#ARGV==0){
	print "\nTakes a group of NMF results and assesses them, outputting the fraction\n";
  print "of nodes assigned a known term. All 'attribute' files must be put in a folder\n";
  print "called 'attributes', and all 'node' files must be put in a folder called\n";
  print "'nodes'.\n";
	print "\nusage:\n $0\n";
	print "-g [GO file with child terms for each GO term]\n";
  print "-l [GO localization file (goa_human.gaf)]\n";
  print "-m [map of GO terms to names]\n";
  print "-n [GO namespace, should be one of 'C' (default), 'F' or 'P']\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-m'){
			$i++;
			$mfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-n'){
			$i++;
			$namespace = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# get child terms for each GO term
print STDERR "Reading hierarchy file\n";
open my $hierarchyFH, '<', $gfile or die "Could not open $gfile: $!";
my $tsv = Text::CSV_XS->new({
	sep_char => "\t",
});
my %children;
while(my $row = $tsv->getline($hierarchyFH)) {
	@{$children{@{$row}[0]}} = split ',', @{$row}[1];
}
close $hierarchyFH;

# get map of GO terms to ids
print STDERR "Reading GO map\n";
open my $mapFH, '<', $mfile or die "Could not open $mfile: $!";
my $mapTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
my %goMap;
while(my $row = $mapTSV->getline($mapFH)) {
	$goMap{@{$row}[1]} = @{$row}[0];
}
close $mapFH;

# reading known localizations
print STDERR "Reading known localizations\n";
open my $localizationfh, '<', $lfile or die "Could not open $lfile: $!";
my $locTSV = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
my %geneLocalizations;
while(my $row = $locTSV->getline($localizationfh)) {
	if (@{$row}[8] eq $namespace) {
		push @{$geneLocalizations{@{$row}[2]}}, @{$row}[4];
	}
}
close $localizationfh;

# read first attribute file (and associated nodes)
print STDERR "Getting attributes and node assigments\n";
opendir my ($attrDir), $dirAttributes or die "Couldn't open dir '$dirAttributes': $!";
my @attrFiles = readdir $attrDir;
closedir $attrDir;
opendir my ($nodeDir), $dirNodes or die "Couldn't open dir '$dirNodes': $!";
my @nodeFiles = readdir $nodeDir;
closedir $nodeDir;
print STDOUT "conditions\tnumber of nodes\tknown assignments\tfraction\tknown assignments (top 5)\tfraction (top 5)\tknown per domain\tnumber of NMF ranks\tworst Jaccard\tbest Jaccard\n";
foreach my $attrFile (@attrFiles) {
  if ($attrFile =~ /(.+)-terms_perrank/) {
    # get terms (and children) per domain
    my %attr;
		my %attrAll;
		my %attrTop5;
		my %termsAdded;
    my $attrHandle = $1;
    $attrFile = $dirAttributes . '/' . $attrFile;
    open my $attrFH, '<', $attrFile or die "Could not open $attrFile: $!";
    my $attrTSV = Text::CSV_XS->new({
    	binary => 1,
    	sep_char => "\t",
    	quote_char => undef,
    	escape_char => undef,
    	allow_loose_quotes => 1,
    	allow_loose_escapes => 1,
    });
    $attrTSV->getline($attrFH); # discard header
    while(my $row = $attrTSV->getline($attrFH)) {
			if (exists $termsAdded{@{$row}[0]}) {
				$termsAdded{@{$row}[0]}++;
			} else {
				$termsAdded{@{$row}[0]} = 1;
			}
			push @{$attr{@{$row}[0]}}, @{$row}[4];
      push @{$attrAll{@{$row}[0]}}, @{$row}[4];
			if (exists $children{@{$row}[4]}) {
        push @{$attrAll{@{$row}[0]}}, @{$children{@{$row}[4]}};
      }
			if ($termsAdded{@{$row}[0]} <= 5) {
				push @{$attrTop5{@{$row}[0]}}, @{$row}[4];
				if (exists $children{@{$row}[4]}) {
	        push @{$attrTop5{@{$row}[0]}}, @{$children{@{$row}[4]}};
	      }
			}
    }
    close $attrFH;
    my $noDomains = keys %attrAll;
    foreach my $domain (keys %attrAll) {
      @{$attrAll{$domain}} = uniq @{$attrAll{$domain}};
			@{$attrTop5{$domain}} = uniq @{$attrTop5{$domain}};
    }
		# calculate Jaccard distance between ranks
		my $jaccardBest = 0;
		my $jaccardWorst = 1;
		for(my $i = 1; $i <= $noDomains; $i++) {
			for(my $j = $i + 1; $j <= $noDomains; $j++) {
				my @intersect = intersect(@{$attr{$i}}, @{$attr{$j}});
				my @union = uniq (@{$attr{$i}}, @{$attr{$j}});
				my $currJaccard = 1 - (scalar @intersect / scalar @union);
				if ($currJaccard < $jaccardWorst) {
					$jaccardWorst = $currJaccard;
				}
				if ($currJaccard > $jaccardBest) {
					$jaccardBest = $currJaccard;
				}
			}
		}
		$jaccardBest = sprintf "%.3f", $jaccardBest;
		$jaccardWorst = sprintf "%.3f", $jaccardWorst;
    # check if node is assigned to known term
    my @indexes = indexes { $_ =~ /^$attrHandle/ } @nodeFiles;
    my $nodeFile = $dirNodes . '/' . $nodeFiles[$indexes[0]];
    open my $nodeFH, '<', $nodeFile or die "Could not open $nodeFile: $!";
    my $nodeTSV = Text::CSV_XS->new({
    	binary => 1,
    	sep_char => "\t",
    	quote_char => undef,
    	escape_char => undef,
    	allow_loose_quotes => 1,
    	allow_loose_escapes => 1,
    });
    $nodeTSV->getline($nodeFH); # discard header
    my $knownNodes = 0;
		my $knownNodesTop5 = 0;
		my @perDomainKnown;
		my @perDomainNodes;
		my $totalNodes = 0;
    while(my $row = $nodeTSV->getline($nodeFH)) {
      $totalNodes++;
      my $assignedDomain = @{$row}[1];
      if ($assignedDomain > 0) {
				if (!($perDomainKnown[$assignedDomain - 1])) {
					$perDomainKnown[$assignedDomain - 1] = 0;
					$perDomainNodes[$assignedDomain - 1] = 0;
				}
				$perDomainNodes[$assignedDomain - 1]++;
        my @knownTerms;
        if (exists $geneLocalizations{@{$row}[0]}) {
          @knownTerms = @{$geneLocalizations{@{$row}[0]}};
        }
        my @availableTerms = @{$attrAll{$assignedDomain}};
				my @availableTermsTop5 = @{$attrTop5{$assignedDomain}};
        my @isect = intersect(@availableTerms, @knownTerms);
				my @isectTop5 = intersect(@availableTermsTop5, @knownTerms);
        if (scalar @isect > 0) {
          $knownNodes++;
					$perDomainKnown[$assignedDomain - 1]++;
        }
				if (scalar @isectTop5 > 0) {
          $knownNodesTop5++;
        }
      }
    }
		my @perDomain;
		for(my $i = 0; $i < scalar @perDomainKnown; $i++){
			$perDomain[$i] = 0;
			if ($perDomainNodes[$i]) {
				$perDomain[$i] = sprintf "%.3f", $perDomainKnown[$i] / $perDomainNodes[$i];
			}
		}
		my $perDomainString = join ',', sort @perDomain;
		my $median = median @perDomain;
		my $fraction = sprintf "%.3f", $knownNodes / $totalNodes;
		my $fractionTop5 = sprintf "%.3f", $knownNodesTop5 / $totalNodes;
    print STDOUT "$attrHandle\t$totalNodes\t$knownNodes\t$fraction\t$knownNodesTop5\t$fractionTop5\t$perDomainString\t$noDomains\t$jaccardWorst\t$jaccardBest\n";
    close $nodeFH;
  }
}
