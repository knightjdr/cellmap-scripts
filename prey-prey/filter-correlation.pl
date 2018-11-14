#!/usr/bin/perl

# Takes a file with gene pairs and a correlation score (for cytoscape),
# removes duplicate pairs and filters based on the specified cutoff.

# 27/5/2018

# libraries
use Text::CSV_XS;

# paramaters
$cutoff = 0.5;

# command line parameters
$gfile = ''; # file with gene pairs

if ($#ARGV==0) {
	print "\nTakes a file with gene pairs and a correlation score (for cytoscape),\n";
	print "removes duplicate pairs and filters based on the specified cutoff.\n\n";
	print "\nusage:\n $0\n";
  print "-c [cutoff]\n";
	print "-f [file with gene pairs (should have header)]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV) {
		if ($ARGV[$i] eq '-c') {
			$i++;
			$cutoff = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-f') {
			$i++;
			$gfile = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Open file for reading.
my $tsv = Text::CSV_XS->new({ sep_char => "\t" });
open my $fh, '<', $gfile or die "Could not open $gfile: $!";
my $header = $tsv->getline($fh); #discard header

# Open file for writing.
(my $filename = $gfile) =~ s/\.[^.]+$//;
open my $outfh, '>', $filename . '_filtered_' . $cutoff . '.txt';

# Output header to new file.
my $headerLine = join "\t", @{$header};
print $outfh "$headerLine\n";

# Iterate over input file and keep rows passing cutoff that have not already
# been printed. Also ignore self pairs.
my %pairs;
while(my $row = $tsv->getline($fh)) {
  my $source = @{$row}[0];
  my $score = @{$row}[2];
  my $target = @{$row}[1];
  if (
    (exists $pairs{$source} && exists $pairs{$source}{$target}) ||
    (exists $pairs{$target} && exists $pairs{$target}{$source})
  ) {
  } elsif (
    $source ne $target &&
    $score >= $cutoff
  ) {
    print $outfh "$source\t$target\t$score\n";
    $pairs{$source}{$target} = 1;
    $pairs{$target}{$source} = 1;
  }
}
close($fh);
close($outfh);
