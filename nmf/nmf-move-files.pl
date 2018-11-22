#!/usr/bin/perl

# 4/7/2017

use strict;
use warnings;

# libraries
use File::Copy qw(copy);

# parameters
my $attributesFile = 'terms_perrank.txt';
my $attributeTargetDir = './nmf_assessment/attributes/';
my $nodeFile = 'gene-localizations.txt';
my $nodeTargetDir = './nmf_assessment/nodes/';
my $resultsFolder ='./Results/';

# create subfolders if needed
mkdir 'nmf_assessment' unless -d 'nmf_assessment';
mkdir 'nmf_assessment/attributes' unless -d 'nmf_assessment/attributes';
mkdir 'nmf_assessment/nodes' unless -d 'nmf_assessment/node';

# get result subfolders and move files
opendir my ($directories), $resultsFolder or die "Couldn't open dir '$resultsFolder': $!";
while (my $entry = readdir $directories) {
  next unless -d $resultsFolder . '/' . $entry;
  next if $entry eq '.' or $entry eq '..';
  # move attribute file
  my $attrFile = $resultsFolder . $entry . '/' . $attributesFile;
  my $newAttrFileHandle = $entry . '-' . $attributesFile;
  my $attrDest = $attributeTargetDir . $newAttrFileHandle;
  copy $attrFile, $attrDest or die "Copy failed: $!";
  my $cmd = "sed -i .bak 1,4d $attrDest";
  system($cmd) == 0 or die "Couldn't launch $cmd: $!/$?";
  unlink $attrDest . '.bak';
  # move node file
  my $ndFile = $resultsFolder . $entry . '/' . $nodeFile;
  my $newNodeFileHandle = $entry . '-' . $nodeFile;
  my $nodeDest = $nodeTargetDir . $newNodeFileHandle;
  copy $ndFile, $nodeDest or die "Copy failed: $!";
}
closedir $directories;
