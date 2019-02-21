#!/usr/bin/perl

# 25/7/2017

use strict;
use warnings;

# load libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;
use JSON::Parse 'parse_json';
use JSON::Parse 'valid_json';
use List::MoreUtils qw(uniq);
use LWP::UserAgent;
use String::Util qw(trim);

# parameters
my $species = 'Homo sapiens';
open my $log, '>', 'log.txt'; # error log

# command line parameters
my $ufile = ''; # uniprot file

if ($#ARGV==0){
	print "\nThis script will get all Pfam motifs for all human uniprot IDs, and\n";
  print "will out those motifs for each gene\n\n";
	print "\nusage:\n $0\n";
	print "-u [Uniprot file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-u'){
			$i++;
			$ufile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# read in uniprot and get uniprot ID and gene namespace
print STDERR "Parsing UniProt data\n";
my $currGene = '';
my $currID = '';
my $isReviewed = 0;
my $isSpecies = 0;
my %uniprotMap;
open my $uniprotFH, "<", $ufile || die "$ufile can't be opened: $!";
while(<$uniprotFH>){
	if ($_ =~ /^ID\s+\S+\s+Reviewed;/) {
		$isReviewed = 1;
	} elsif(
    $_ =~ /^OS   $species/ &&
    $isReviewed &&
    !$isSpecies
  ) {
		$isSpecies = 1;
	} elsif(
    $_ =~ /^GN   Name=([^;{]+)/ &&
    $isReviewed &&
    !$currGene
  ) {
		$currGene = trim($1);
	} elsif(
    $_ =~ /^GN   OrderedLocusNames=([^;{]+)/ &&
    $isReviewed &&
    !$currGene
  ) {
		$currGene = trim($1);
	} elsif(
    $_ =~ /^GN   ORFNames=([^;{]+)/ &&
    $isReviewed &&
    !$currGene
  ) {
		$currGene = trim($1);
	} elsif(
    $_ =~ /^AC\s+([A-Z0-9]+);/ &&
    $isReviewed &&
    !$currID
  ) {
		$currID = $1;
	} elsif($_ =~ /^\/\//) {
    if(
      $isReviewed &&
      $isSpecies
    ) {
			$uniprotMap{$currID} = $currGene;
		}
		#reset some params
		$currGene = '';
		$currID = '';
		$isReviewed = 0;
		$isSpecies = 0;
	}
}
close $uniprotFH;

# get Pfam information
print STDERR "Retrieving Pfam information\n";
#create user agent
my $ua = LWP::UserAgent->new;
$ua->env_proxy;
# init parameters for retreival
my $countQueries = 0;
my $numberQueries = keys %uniprotMap;
my %regions;
# output
open my $outputFH, '>', 'pfam_motifs.txt';
open my $outputCompleteFH, '>', 'pfam_motifs_complete.txt';
print $outputFH "gene\tmotif\n";
print $outputCompleteFH "gene\tmotif\tstart\tend\tmotif length\tprotein length\n";
foreach my $id (keys %uniprotMap) {
  # get JSON and parse
  my $res = $ua->get( 'http://pfam.xfam.org/protein/' . $id . '/graphic' );
  if ($res->is_success) {
    my $response = $res->content;
    if (valid_json $response) {
			my $currGene = $uniprotMap{$id};
      my $json = parse_json $response;
			my $proteinLength = $json->[0]{'length'};
      my @motifs = @{$json->[0]{'motifs'}};
      # add motifs
      foreach my $motif (@motifs) {
        push @{$regions{$id}}, $motif->{'type'};
				my $motifLength = $motif->{'end'} - $motif->{'start'} + 1;
				print $outputCompleteFH "$currGene\t$motif->{'type'}\t$motif->{'start'}\t$motif->{'end'}\t$motifLength\t$proteinLength\n";
      }
      # only report unique motifs (remove this line if duplicates allowed)
      @{$regions{$id}} = uniq @{$regions{$id}};
      # output current results
      foreach my $motif (@{$regions{$id}}) {
        print $outputFH "$currGene\t$motif\n";
      }
    } else {
      print $log "Invalid JSON for $id\n";
    }
  } else {
    print $log "error: failed to retrieve JSON for $id: " . $res->status_line . "\n";
  }
  $countQueries++;
  if ($countQueries % 100 == 0) {
    print STDERR "$countQueries of $numberQueries complete\n";
  }
}
close $outputFH;
close $outputCompleteFH;
