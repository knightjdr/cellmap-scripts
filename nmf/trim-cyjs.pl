#!/usr/bin/perl

# 19/6/2018

use strict;
use warnings;

# libraries

# Command line parameters.
my $cfile = ''; # cyjs file.

if ($#ARGV==0){
	print "\nTakes a .cyjs file and removes coordinates attributes from data object\n";
	print "\nusage:\n $0\n";
	print "-c [.cyjs file]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-c'){
			$i++;
			$cfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Vars.
my $data = 0; # Boolean to see if the current line is in the data object.

# Read .cyjs file
print STDERR "Reading .cyjs file\n";
open my $cyjsFH, '<', $cfile or die "Could not open $cfile: $!";
open my $outputFH, '>', 'trimmed.cyjs';
while(my $row = <$cyjsFH>) {
  my $print = 1; # Should a line be printed?
  if ($row =~ /"data" : {/) {
    $data = 1;
  } elsif ($row =~ /"position" : {/) {
    $data = 0;
  } elsif (
    $data == 1 &&
    (
      $row =~ /"x\d*" :/ ||
      $row =~ /"y\d*" :/
    )
  ) {
    $print = 0;
  }
  if ($print == 1) {
    print $outputFH $row;
  }
}
close $cyjsFH;
close $outputFH;
