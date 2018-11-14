#!/usr/bin/perl

# 4/7/2017

use strict;
use warnings;

# libraries
use Array::Utils qw(:all);
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(indexes);
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# parameters
my $dirAttributes = 'attributes';
my $dirNodes = 'nodes';
my $namespace = 'C';

# command line parameters
my $gfile = ''; # file with all child terms for each GO term
my $lfile = ''; # goa_human.gaf
my $mfile = ''; # map of GO terms to ID's

if ($#ARGV==0){
	print "\nTakes a group of SAFE results and assesses them, outputting the fraction\n";
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
print STDOUT "conditions\tnumber of nodes\tknown assignments\tfraction\tnumber of SAFE domains\n";
foreach my $attrFile (@attrFiles) {
  if ($attrFile =~ /(.+)-attribute_properties_annotation/) {
    # get terms (and children) per domain
    my %attr;
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
      push @{$attr{@{$row}[2]}}, $goMap{@{$row}[1]};
      if (exists $children{$goMap{@{$row}[1]}}) {
        push @{$attr{@{$row}[2]}}, @{$children{$goMap{@{$row}[1]}}};
      }
    }
    close $attrFH;
    my $noDomains = keys %attr;
    $noDomains++;
    foreach my $domain (keys %attr) {
      @{$attr{$domain}} = uniq @{$attr{$domain}};
    }
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
    my $totalNodes = 0;
    my $knownNodes = 0;
    while(my $row = $nodeTSV->getline($nodeFH)) {
      $totalNodes++;
      my $assignedDomain = @{$row}[2];
      if ($assignedDomain > 1) {
        my @knownTerms;
        if (exists $geneLocalizations{@{$row}[0]}) {
          @knownTerms = @{$geneLocalizations{@{$row}[0]}};
        }
        my @availableTerms = @{$attr{$assignedDomain}};
        my @isect = intersect(@availableTerms, @knownTerms);
        if (scalar @isect > 0) {
          $knownNodes++;
        }
      }
    }
		my $fraction = sprintf "%.3f", $knownNodes / $totalNodes;
    print STDOUT "$attrHandle\t$totalNodes\t$knownNodes\t$fraction\t$noDomains\n";
    close $nodeFH;
  }
}
