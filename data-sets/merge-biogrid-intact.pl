#!/usr/bin/perl

# 21/6/2018

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Text::CSV_XS;

# command line parameters
my $bfile = ''; # BioGRID file
my $ifile = ''; # Intact file

if($#ARGV==0){
	print "\nMerges BioGrid and Intact interactor lists into a single file.\n";
  print "The intact file must first be run through intact-parse.pl.\n";
	print "\nusage:\n $0\n";
	print "-b [BioGRID]\n";
	print "-i [IntAct file]\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-b'){
			$i++;
			$bfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-i'){
			$i++;
			$ifile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Vars
my %interactions;

# parse BioGRID file
print STDERR "Retrieving Biogrid interactors\n";
open my $biogridFH, '<', $bfile or die "Could not open $bfile: $!";
my $biogridTsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
$biogridTsv->getline($biogridFH); #discard header
while(my $row = $biogridTsv->getline($biogridFH)) {
  my $approach = @{$row}[11];
  my $source = lc @{$row}[7];
  my $target = lc @{$row}[8];
  my @pair = ($source, $target);
  @pair = sort @pair;
  push @{$interactions{$pair[0]}{$pair[1]}}, $approach;
}
close($biogridFH);

# parse Intact file
print STDERR "Retrieving Intact interactors\n";
open my $intactFH, '<', $ifile or die "Could not open $ifile: $!";
my $intactTsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
$intactTsv->getline($intactFH); #discard header
while(my $row = $intactTsv->getline($intactFH)) {
  my $approach = @{$row}[2];
  my $source = lc @{$row}[0];
  my $target = lc @{$row}[1];
  my @pair = ($source, $target);
  @pair = sort @pair;
  push @{$interactions{$pair[0]}{$pair[1]}}, $approach;
}
close($intactFH);

# Output list of interactions.
open my $outputFH, '>', 'interactions.txt';
print $outputFH "source\ttarget\tapproach\n";
foreach my $source (sort {lc $a cmp lc $b} keys %interactions) {
  foreach my $target (sort {lc $a cmp lc $b} keys %{$interactions{$source}}) {
    my $approachString = join ";", uniq @{$interactions{$source}{$target}};
    print $outputFH "$source\t$target\t$approachString\n";
  }
}
