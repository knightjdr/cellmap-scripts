#!/usr/bin/perl

# 2/26/2019

use strict;
use warnings;

opendir my $dir, './' or die "Cannot open directory: $!";
my @files = grep {/\.json$/ and not /\.min\.json$/} readdir $dir;
closedir $dir;

foreach my $file (@files) {
  my $minFile = $file =~ s/json/min.json/r;
  print STDOUT "$minFile\n";
  my $cmd = "jq -c . < $file > $minFile";
  system($cmd) == 0 or die "Couldn't launch $cmd: $!/$?";
}
