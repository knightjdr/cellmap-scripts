#!/usr/bin/perl

# 4/7/2017

use strict;
use warnings;

# libraries
use File::Copy qw(copy);

# parameters
my $attributesFile = 'attribute_properties_annotation-highest.txt';
my $attributeTargetDir = './safe_assessment/attributes/';
my $nodeFile = 'node_properties_annotation-highest.txt';
my $nodeTargetDir = './safe_assessment/nodes/';
my $resultsFolder ='./Results/';

# create subfolders if needed
mkdir 'safe_assessment' unless -d 'safe_assessment';
mkdir 'safe_assessment/attributes' unless -d 'safe_assessment/attributes';
mkdir 'safe_assessment/nodes' unless -d 'safe_assessment/node';

# get result subfolders and move files
opendir my ($directories), $resultsFolder or die "Couldn't open dir '$resultsFolder': $!";
while (my $entry = readdir $directories) {
  next unless -d $resultsFolder . '/' . $entry;
  next if $entry eq '.' or $entry eq '..';
  # move attribute file
  my $attrFile = $resultsFolder . $entry . '/' . $attributesFile;
  my $newAttrFileHandle = $entry . '-' . $attributesFile;
  my $attrDest = $attributeTargetDir . $newAttrFileHandle;
  copy $attrFile, $attrDest or next;
  my $cmd = "sed -i .bak 1,4d $attrDest";
  system($cmd) == 0 or die "Couldn't launch $cmd: $!/$?";
  unlink $attrDest . '.bak';
  # move node file
  my $ndFile = $resultsFolder . $entry . '/' . $nodeFile;
  my $newNodeFileHandle = $entry . '-' . $nodeFile;
  my $nodeDest = $nodeTargetDir . $newNodeFileHandle;
  copy $ndFile, $nodeDest or die "Copy failed: $!";
  $cmd = "sed -i .bak 1,4d $nodeDest";
  system($cmd) == 0 or die "Couldn't launch $cmd: $!/$?";
  unlink $nodeDest . '.bak';
}
closedir $directories;
