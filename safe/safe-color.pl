#!/usr/bin/perl

# 4/6/2018

use strict;
use warnings;

# modules
use Math::Round;
use Text::CSV_XS;

# paremeters
my @colors = (
	"#FFFFFF", "#1CE6FF", "#FF4A46", "#008941", "#006FA6", "#8FB0FF", "#FFB500",
	"#809693", "#4FC601", "#3B5DFF", "#00C2A0", "#FFAA92", "#D16100", "#7B4F4B",
	"#A1C299", "#0AA6D8", "#00846F", "#997D87", "#A079BF", "#C0B9B2", "#C2FF99",
	"#00489C", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
	"#885578", "#FAD09F", "#FF8A9A", "#BEC459", "#456648", "#0086ED", "#886F4C",
	"#B4A8BD", "#00A6AA", "#636375", "#A3C8C9", "#FF913F", "#938A81", "#00FECF",
	"#B05B6F", "#8CD0FF", "#C8A1A1", "#A77500", "#6367A9", "#A05837", "#D790FF",
	"#9B9700", "#549E79", "#72418F", "#99ADC0", "#0089A3", "#CB7E98", "#324E72",
	"#83AB58", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#BF5650", "#66796D",
	"#8ADBB4", "#C895C5", "#FF6832", "#66E1D3", "#CFCDAC", "#D0AC94", "#7ED379"
);

# command line parameters
my $nfile = '';	# node properties file
my $sfile = ''; # summary file

if ($#ARGV == 0) {
	print "\nTakes SAFE output on nodes and generates a file to\n";
	print "use for colouring those nodes in cytoscape.\n";
	print "\nusage:\n $0\n";
	print "-n [node_properties_annotation-highest_noHeader.txt]\n";
	print "-s [summary file]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV){
		if ($ARGV[$i] eq '-n'){
			$i++;
			$nfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s'){
			$i++;
			$sfile = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# Read domain names.
my @domainNames;
push @domainNames, '';
open my $domainFH, '<', $sfile or die "Could not open $sfile: $!";
my $domainTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$domainTSV->getline($domainFH); # discard header
while(my $row = $domainTSV->getline($domainFH)) {
	my ($domain) = @{$row}[1] =~ /\[(.+)\]/;
	push  @domainNames, $domain;
}
close $domainFH;

# Read node file.
open my $colorFH, '>', 'domain-color.tsv';
print $colorFH "name\tdomain\tdomain name\tcolor\n";
open my $nodeFH, '<', $nfile or die "Could not open $nfile: $!";
my $nodeTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$nodeTSV->getline($nodeFH); # discard header
while(my $row = $nodeTSV->getline($nodeFH)) {
	my $domain = @{$row}[2];
	print $colorFH "@{$row}[0]\t$domain\t$domainNames[$domain]\t$colors[$domain]\n";
}
close $colorFH;
