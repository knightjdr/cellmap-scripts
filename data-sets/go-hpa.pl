#!/usr/bin/perl

# 21/3/2019

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# command line parameters
my $gfile = ''; # GO GAF file

if($#ARGV==0){
	print "Removes HPA annotations from GO .gaf file\n\n";
	print "\nusage:\n $0\n";
  print "-g [goa_human.gaf without header]\n\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
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
open my $gafFH, '<', $gfile or die "Could not open $gfile: $!";
open my $outputFH, '>', 'goa_human_nohpa.gaf';
while(my $row = $tsv->getline($gafFH)) {
  my $database = @{$row}[14];
  if ($database ne 'HPA') {
    my $line = join "\t", @{$row};
    print $outputFH "$line\n";
  }
}

close $gafFH;
close $outputFH;
