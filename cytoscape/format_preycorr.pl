#!/usr/bin/perl

#7/7/2016

use strict;
use warnings;

my $cfile = ''; 	#.corr file

if($#ARGV==0){
	print "\nTakes a prey-prey correlation file and removes self-links and duplicate links\n";
	print "\nusage:\n $0\n";
	print "-c [corr file]\n";
	die "\n";
}
else{
	my $i = 0;
	while($i<=$#ARGV){
		if($ARGV[$i] eq '-c'){
			$i++;
			$cfile = $ARGV[$i];
		}
		else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

my %pairs = ();

open(CFILE, "<$cfile") || die "$cfile can't be opened: $!";
while(<CFILE>){
	if ($_ =~ /^gene/) {
	} elsif($_ =~ /^(\S+)\t(\S+)\t(\S+)/) {
		my $source = $1;
		my $target = $2;
		my $edge = $3;
		$source =~ s/\*//g;
		$target =~ s/\*//g;
		#sort source and target (easy way to check for duplicates)
		($source, $target) = sort { lc($a) cmp lc($b) } ($source, $target);
		if ($source ne $target) {
			if (exists $pairs{$source}) {
				if (!exists $pairs{$source}{$target}) {
					$pairs{$source}{$target} = $edge;
				}
			} else {
				$pairs{$source} = ();
				$pairs{$source}{$target} = $edge;
			}
		}
	}
	else{
	}
}
close(CFILE);

open(OFILE, ">nodes.txt");
print OFILE "source\ttarget\tedge\n";
foreach my $source (keys %pairs) {
	foreach my $target (keys %{$pairs{$source}}) {
		print OFILE "$source\t$target\t$pairs{$source}{$target}\n";
	}
}
close(OFILE);
