#!/usr/bin/perl

# Takes output from tSNE matlab script (x, y coordinates for genes) and
# generates a file to import to Cytoscape.

# 6/11/2018

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# command line parameters
my $colorFile = ''; # File with color and rank for each gene.
my $detailsFile = ''; # File with name for each NMF rank, called nmf_summary.txt.
my $tsneFile = ''; # File output from tSNE Matlab script, called tsne_nmf.txt

if ($#ARGV == 0) {
	print "\nTakes tSNE output file and formats it for Cytoscape, adding colors\n";
  print "\nand multiple the coordinates by 30, 50 and 100.\n";
	print "\nusage:\n $0\n";
	print "-c [color file]\n";
  print "-d [nmf_summary.txt]";
	print "-t [tSNE file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-c'){
			$i++;
			$colorFile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-d'){
			$i++;
			$detailsFile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
			$i++;
			$tsneFile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Vars.
my %genes;
my %ranks;

# Read tSNE file.
print STDERR "Reading tSNE\n";
open my $tsneFH, '<', $tsneFile or die "Could not open $tsneFile: $!";
my $tsneTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$tsneTSV->getline($tsneFH); # Discard header.
while(my $row = $tsneTSV->getline($tsneFH)) {
  $genes{@{$row}[0]}{'x'} = @{$row}[1];
  $genes{@{$row}[0]}{'y'} = @{$row}[2];
}
close $tsneFH;

# Read details file.
print STDERR "Reading NMF details\n";
open my $detailsFH, '<', $detailsFile or die "Could not open $detailsFile: $!";
my $detailsTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$detailsTSV->getline($detailsFH); # Discard header.
while(my $row = $detailsTSV->getline($detailsFH)) {
  my ($term) = @{$row}[1] =~ /\[(.+)\]/;
  $ranks{@{$row}[0]} = $term;
}
close $detailsFH;

# Read color file and add to gene hash.
print STDERR "Reading colors\n";
open my $colorFH, '<', $colorFile or die "Could not open $colorFile: $!";
my $colorTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$colorTSV->getline($colorFH); # Discard header.
while(my $row = $colorTSV->getline($colorFH)) {
  my $currGene = @{$row}[0];
  if (exists $genes{$currGene}) {
    my $currRank = @{$row}[1];
    $genes{$currGene}{'color'} = @{$row}[2];
    $genes{$currGene}{'rank'} = $currRank;
    $genes{$currGene}{'term'} = $ranks{$currRank};
    $genes{$currGene}{'trans'} = @{$row}[3];
  }
}
close $colorFH;

# Print out file for Cytoscape.
open my $cytoscapeFH, '>', 'tsne_nmf_details.txt';
print $cytoscapeFH "gene\tx\ty\trank\trank name\tcolor\ttransparency\tx30\tx50\tx100\ty30\ty50\ty100\n";
foreach my $gene (sort {lc $a cmp lc $b} keys %genes) {
  my $x30 = $genes{$gene}{'x'} * 30;
  my $x50 = $genes{$gene}{'x'} * 50;
  my $x100 = $genes{$gene}{'x'} * 100;
  my $y30 = $genes{$gene}{'y'} * 30;
  my $y50 = $genes{$gene}{'y'} * 50;
  my $y100 = $genes{$gene}{'y'} * 100;
  print $cytoscapeFH "$gene\t$genes{$gene}{'x'}\t$genes{$gene}{'y'}\t$genes{$gene}{'rank'}\t$genes{$gene}{'term'}\t$genes{$gene}{'color'}\t$genes{$gene}{'trans'}\t$x30\t$x50\t$x100\t$y30\t$y50\t$y100\n";
}
close $cytoscapeFH;
