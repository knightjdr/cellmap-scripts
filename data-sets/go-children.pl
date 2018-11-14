#!/usr/bin/perl

# 5/6/2017

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# parameters
my $namespace = 'cc';

# command line parameters
my $hfile = ''; # GO heirarchy file (obo format)

if ($#ARGV==0){
	print "\nTakes GO hierarchy in obo format and for each term it gets its children\n";
	print "and outputs them to a file for use in downstream scripts.\n";
	print "\nusage:\n $0\n";
  print "-h [.obo (hierarchy) file]\n";
	print "-n [a = all, bp = biological process, cc = cellular component (default), mf = molecular function]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-h'){
			$i++;
			$hfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-n'){
			$i++;
			$namespace = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# subs
sub getChildren {
  my @terms = @{$_[0]};
	my %parents = %{$_[1]};
  my @children = @{$_[2]};
	my @newChildren;
  my $print = 0;
  my @newTerms;
  for(my $i = 0; $i < scalar @terms; $i++) {
    foreach my $child (keys %parents) {
      if (grep /^$terms[$i]$/, @{$parents{$child}}) {
        push @newTerms, $child;
        push @newChildren, $child;
      }
    }
  }
	# because graph is acyclic, need to remove any new terms already checked
	if (scalar @children > 0) {
		my %childTest = map { $_ => 1 } @children;
		for(my $i = scalar @newTerms - 1; $i >= 0; $i--) {
			if (exists $childTest{$newTerms[$i]}) {
				splice @newTerms, $i, 1;
			}
		}
	}
	push @children, @newChildren;
  if (scalar @newTerms != 0) {
    @children = getChildren(\@newTerms, \%parents, \@children);
  }
  return uniq @children;
}

sub reportProgress {
  my ($current, $interval, $total) = @_;
  if ($current % $interval == 0) {
    my $noIntervals = $current / $interval;
    print STDERR "[";
    for (my $j = 0; $j < $noIntervals; $j++) {
      print STDERR ".";
    }
    for (my $j = 0; $j < int($total / $interval) + 1 - $noIntervals; $j++) {
      print STDERR " ";
    }
    print STDERR "]\n";
  }
}

# set namespace
my $longNameSpace;
if ($namespace eq "bp") {
	$longNameSpace = "biological_process";
} elsif ($namespace eq "cc") {
	$longNameSpace = "cellular_component";
} elsif ($namespace eq "mf") {
	$longNameSpace = "molecular_function";
} else {
	$longNameSpace = "all";
}

# get direct parent terms for each GO term
open my $hierarchyFH, '<', $hfile or die "Could not open $hfile: $!";
my $currId;
my $currName;
my $inNamespace;
my %idMap = ();
my @goTerms;
my $obsolete = 0;
my %parents = ();
print STDERR "Reading GO hierarchy\n";
open my $gotermsfh, '>', 'go-terms.txt';
while(<$hierarchyFH>) {
	if ($_ =~ /^\[Term\]/) {
		$inNamespace = 0;
		$obsolete = 0;
	} elsif ($_ =~ /^id: (\S+)/) {
		$currId = $1;
	} elsif ($_ =~ /^name: obsolete/) {
		$obsolete = 1;
	} elsif ($_ =~ /^name: (.+)/) {
		$currName = $1;
	} elsif ($_ =~ /^namespace: (\S+)/) {
		if (
			($longNameSpace eq 'all' || $longNameSpace eq $1) &&
			$obsolete == 0
		) {
			$inNamespace = 1;
			$idMap{$currId} = $currName;
			push @goTerms, $currId;
			print $gotermsfh "$currId\n";
		}
	} elsif ($_ =~ /^is_a: (\S+) ! .+/) {
		my $parent = $1;
		if ($inNamespace && $obsolete == 0) {
			push @{$parents{$currId}}, $parent;
		}
	} elsif ($_ =~ /^relationship: part_of (\S+) ! .+/) {
		my $parent = $1;
		if ($inNamespace && $obsolete == 0) {
			push @{$parents{$currId}}, $parent;
		}
	}
	else {
	}
}
close $hierarchyFH;
close $gotermsfh;
@goTerms = uniq @goTerms;

# get child terms for each GO ID
my %children;
my $noGoTerms = scalar @goTerms;
my $progress = 0;

print STDERR "Getting child terms for $noGoTerms total terms\n";
for(my $i = 0; $i < $noGoTerms; $i++) {
	my @currTerm;
	push @currTerm, $goTerms[$i];
	@{$children{$goTerms[$i]}} = getChildren(\@currTerm, \%parents, \@{$children{$goTerms[$i]}});
	$progress++;
	reportProgress($progress, 100, $noGoTerms);
}

# output children
open my $outputfh, '>', 'go-children_' . $namespace . '.txt';
foreach my $term (keys %children ){
	print $outputfh "$term\t";
	if (exists $children{$term}) {
		my $joinedTerms = join ',', @{$children{$term}};
		print $outputfh "$joinedTerms\n";
	} else {
		print $outputfh "\n";
	}
}
close $outputfh;
