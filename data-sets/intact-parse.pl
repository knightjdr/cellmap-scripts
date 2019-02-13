#!/usr/bin/perl

# 1/6/2017

use strict;
use warnings;

# libraries
use Data::Dumper; # use like this to print an array print Dumper \@array;

use Text::CSV_XS;

# paramaters
my $fileType = 'b';
my $organism = 9606;
my $requiredEvidence = 1;

# command line parameters
my $ifile = ''; # IntAct file

if($#ARGV==0){
	print "Takes an IntAct database file and parses it for human genes\n\n";
	print "\nusage:\n $0\n";
	print "-i [IntAct file]\n";
  print "-o [organism, default: 9606]\n\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-i'){
			$i++;
			$ifile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

my $tsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
open my $fh, '<', $ifile or die "Could not open $ifile: $!";
open my $parsefh, '>', 'human-intact.txt';
print $parsefh "source\ttarget\tapproach\n";
$tsv->getline($fh); #discard header
while(my $row = $tsv->getline($fh)) {
	my $sourceSpeciesCell = @{$row}[9];
	my $targetSpeciesCell = @{$row}[10];
	my ($sourceSpecies) = $sourceSpeciesCell =~ /^taxid:([-\d]+)\(/;
	my ($targetSpecies) = $targetSpeciesCell =~ /^taxid:([-\d]+)\(/;
	if ($sourceSpecies &&
		$targetSpecies &&
		(
			($sourceSpecies == $organism && $targetSpecies > 0) ||
			($targetSpecies == $organism && $sourceSpecies > 0)
		)
	) {
		my $sourceCell = @{$row}[4];
		my $targetCell = @{$row}[5];
		my ($source) = $sourceCell =~ /uniprotkb:([^\(]+)\(gene name\)/;
		my ($target) = $targetCell =~ /uniprotkb:([^\(]+)\(gene name\)/;
		if($source && $target) {
			my $approachCell = @{$row}[6];
			my ($approach) = $approachCell =~ /\(([^\(^\)]+)\)/;
			if(!$approach) {
				$approach = 'missing approach';
				print STDERR "missing approach\n";
			}
			print $parsefh "$source\t$target\t$approach\n";
		}
	}
}
close $fh;
close $parsefh;
