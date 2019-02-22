#!/usr/bin/perl

# 23/6/2017

use strict;
use warnings;

# libraries
use Array::Utils qw(:all);
use lib qw(..);
use JSON::XS qw( );
use Text::CSV_XS;

# parameters
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
my %transparency = (
	"100" => "FF",
	"95" => "F2",
	"90" => "E6",
	"85" => "D9",
	"80" => "CC",
	"75" => "BF",
	"70" => "B3",
	"65" => "A6",
	"60" => "99",
	"55" => "8C",
	"50" => "80",
	"45" => "73",
	"40" => "66",
	"35" => "59",
	"30" => "4D",
	"25" => "40",
	"20" => "33",
	"15" => "26",
	"10" => "1A",
	"5" => "0D",
	"0" => "00",
);
my $fileType = 'n';

# command line parameters
my $afile = ''; # annotations matrix
my $gfile = ''; # rank assigned to each localization
my $nfile = ''; # network in .cyjs format
my $outputPrefix = 'corr'; # prefix for output files
my $rfile = ''; # ranks.json file
my $sfile = ''; # SAINT file

if($#ARGV==0) {
	print "\tCreates a network for the cell map.\n";
	print "\nusage:\n $0\n";
	print "-a [annotations matrix]\n";
  print "-g [rank associated with each node]\n";
	print "-n [network in .cyjs format]\n";
  print "-o [output prefix]\n";
	print "-r [ranks.json file]\n";
	print "-s [SAINT file]\n";
	print "-t [n for NMF (default) or s for SAFE]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV) {
		if($ARGV[$i] eq '-a') {
			$i++;
			$afile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-g') {
			$i++;
			$gfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-n') {
			$i++;
			$nfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-o') {
      $i++;
      $outputPrefix = $ARGV[$i];
    } elsif($ARGV[$i] eq '-r') {
			$i++;
			$rfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-s') {
			$i++;
			$sfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-t') {
			$i++;
			$fileType = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

#get baits
print STDERR "Getting baits\n";
my %baits;
open my $saintFH, '<', $sfile or die "Could not open $sfile: $!\n";
my $saintTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$saintTSV->getline($saintFH); # discard header
while(my $row = $saintTSV->getline($saintFH)) {
	$baits{@{$row}[0]} = 1;
}
close $saintFH;

#get node categories for each gene (the score is only for SAFE, and can be used for shading nodes)
print STDERR "Get node categories\n";
open my $nodeFH, '<', $gfile or die "Could not open $gfile: $!\n";
my $nodeTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$nodeTSV->getline($nodeFH); # discard header
my %categories;
while(my $row = $nodeTSV->getline($nodeFH)) {
  my $domain;
  my $node;
  my $score = 1;
  if ($fileType eq 'n') {
    $domain = @{$row}[1];
    $node = @{$row}[0];
		$score = @{$row}[2];
  } else {
    $domain = @{$row}[2];
    $node = @{$row}[1];
    $score = @{$row}[3];
  }
  $categories{$node}{'domain'} = $domain;
  $categories{$node}{'score'} = $score;
  $categories{$node}{'known_name'} = 'false';
  $categories{$node}{'known_term'} = 'false';
  $categories{$node}{'bait'} = 'false';
  if (exists $baits{$node} ) {
    $categories{$node}{'bait'} = 'true';
  }
}
close $nodeFH;

# get rank names
print STDERR "Getting rank info\n";
my @rankDiseases;
my @rankDomains;
my @rankMotifs;
my @rankNames;
my @rankListed;
my @rankTerms;
my $json_text = do {
  open my $json_fh, '<:encoding(UTF-8)', $rfile or die "Could not open $rfile: $!\n";
  local $/;
  <$json_fh>
};
my $json = JSON::XS->new;
my $rankInfo = $json->decode($json_text);
for(@{$rankInfo}) {
	my $currRank;
	if ($fileType eq 'n') {
		$currRank = $_->{rank};
	} else {
		$currRank = $_->{domain};
	}
	my @currTerms;
	for(@{$_->{terms}}) {
		push @currTerms, $_->{go_name};
	}
	my @currDiseases;
	for(@{$_->{diseases}}) {
		push @currDiseases, $_->{disease_name};
	}
	my @currDomains;
	for(@{$_->{domains}}) {
		push @currDomains, $_->{domain_name};
	}
	my @currMotifs;
	for(@{$_->{motifs}}) {
		push @currMotifs, $_->{motif_name};
	}
	@{$rankTerms[$currRank]} = @currTerms;
	@{$rankDiseases[$currRank]} = @currDiseases;
	@{$rankDomains[$currRank]} = @currDomains;
	@{$rankMotifs[$currRank]} = @currMotifs;
  @{$rankNames[$currRank]} = @{$_->{names}};
	@{$rankListed[$currRank]} = @{$_->{listedname}};
}

#read in annotations file and check if annotation was known
print STDERR "Reading annotations\n";
open my $annotationFH, '<', $afile or die "Could not open $afile: $!\n";
my @headers;
while(<$annotationFH>) {
	if ($_ =~ /^\t(.+)/) {
		@headers = split(/\t/, $1);
	} elsif ($_ =~ /^(\S+)\t(.+)/) {
		my $node = $1;
		my @attributes = split(/\t/, $2);
		if (exists $categories{$node}) {
			my @currDomains;
			for(my $i = 0; $i < scalar @attributes; $i++) {
				if($attributes[$i] == 1) {
					push @currDomains, $headers[$i];
				}
			}
			my $assignedDomain = $categories{$node}{'domain'};
			my @isect = intersect(@{$rankNames[$assignedDomain]}, @currDomains);
			if (scalar @isect > 0) {
				$categories{$node}{'known_name'} = 'true';
			}
      @isect = intersect(@{$rankTerms[$assignedDomain]}, @currDomains);
			if (scalar @isect > 0) {
				$categories{$node}{'known_term'} = 'true';
			}
		}
	}
}
close $annotationFH;

#get node coordinates
print STDERR "Read node coordinates\n";
my $cyjsText = do {
  open my $json_fh, '<:encoding(UTF-8)', $nfile or die "Could not open $nfile: $!\n";
  local $/;
  <$json_fh>
};
my $cyjs = JSON::XS->new->decode($cyjsText);
my @edgeArray;
my @nodeArray;
my %nodes;
foreach my $node (@{$cyjs->{elements}{nodes}}) {
  my %nodeHash = (
    'name' => $node->{data}{name},
    'x' => $node->{position}{x},
    'y' => $node->{position}{y},
  );
	$nodes{$node->{data}{name}}{'x'} = $nodeHash{x};
	$nodes{$node->{data}{name}}{'y'} = $nodeHash{y};
  push @nodeArray, \%nodeHash;
}
foreach my $edge (@{$cyjs->{elements}{edges}}) {
  my ($gene1, $gene2) = $edge->{data}{shared_name} =~ /(\S+) \(interacts with\) (\S+)/;
  my %edgeHash = (
    'x1' => $nodes{$gene1}{'x'},
    'y1' => $nodes{$gene1}{'y'},
    'x2' => $nodes{$gene2}{'x'},
    'y2' => $nodes{$gene2}{'y'},
  );
  push @edgeArray, \%edgeHash;
}

# output network
print STDERR "Output network\n";
open my $networkFH, '>', $outputPrefix . '-network.json';
print $networkFH "{\n\t\"nodes\": [\n";
for(my $i = 0; $i < scalar @nodeArray; $i++) {
  my %currNode = %{$nodeArray[$i]};
  my $x = sprintf "%.2f", $currNode{'x'};
  my $y = sprintf "%.2f", $currNode{'y'};
	print $networkFH "\t\t{\n\t\t\t";
	print $networkFH "\"name\": \"$currNode{'name'}\",\n\t\t\t";
	print $networkFH "\"x\": $x,\n\t\t\t";
	print $networkFH "\"y\": $y,\n\t\t\t";
	if (exists $categories{$currNode{'name'}}) {
		print $networkFH "\"domain\": $categories{$currNode{'name'}}{'domain'},\n\t\t\t";
		print $networkFH "\"color\": \"$colors[$categories{$currNode{'name'}}{'domain'}]\",\n\t\t\t";
		print $networkFH "\"score\": $categories{$currNode{'name'}}{'score'},\n\t\t\t";
		print $networkFH "\"bait\": $categories{$currNode{'name'}}{'bait'},\n\t\t\t";
    print $networkFH "\"known_name\": $categories{$currNode{'name'}}{'known_name'},\n\t\t\t";
		print $networkFH "\"known_term\": $categories{$currNode{'name'}}{'known_term'}\n\t\t}";
	} else {
		print $networkFH "\"domain\": 0,\n\t\t\t";
		print $networkFH "\"color\": \"$colors[0]\",\n\t\t\t";
		print $networkFH "\"score\": 0,\n\t\t\t";
		print $networkFH "\"bait\": \"prey\",\n\t\t\t";
		print $networkFH "\"known\": 1\n\t\t}";
	}
	if ($i < scalar @nodeArray - 1) {
		print $networkFH ",\n";
	} else {
		print $networkFH "\n";
	}
}
print $networkFH "\t],\n";


print $networkFH "\t\"edges\": [\n";
for(my $i = 0; $i < scalar @edgeArray; $i++) {
  my $x1 = sprintf "%.2f", $edgeArray[$i]{'x1'};
  my $y1 = sprintf "%.2f", $edgeArray[$i]{'y1'};
  my $x2 = sprintf "%.2f", $edgeArray[$i]{'x2'};
  my $y2 = sprintf "%.2f", $edgeArray[$i]{'y2'};
	print $networkFH "\t\t{\n\t\t\t";
	print $networkFH "\"x1\": $x1,\n\t\t\t";
	print $networkFH "\"y1\": $y1,\n\t\t\t";
	print $networkFH "\"x2\": $x2,\n\t\t\t";
	print $networkFH "\"y2\": $y2\n\t\t}";
	if ($i < scalar @edgeArray - 1) {
		print $networkFH ",\n";
	}	else {
		print $networkFH "\n";
	}
}
print $networkFH "\t],\n";
print $networkFH "\t\"domains\": [\n";
for(my $i = 1; $i < scalar @rankNames; $i++) {
	print $networkFH "\t\t{\n\t\t\t";
	print $networkFH "\"domain\": $i,\n\t\t\t";
	print $networkFH "\"color\": \"$colors[$i]\",\n\t\t\t";
	print $networkFH "\"terms\": [";
	for(my $j = 0; $j < scalar(@{$rankTerms[$i]}); $j++) {
		print $networkFH "\"$rankTerms[$i][$j]\"";
		if ($j < scalar(@{$rankTerms[$i]}) - 1) {
			print $networkFH ", ";
		}
	}
	print $networkFH "],\n";
	print $networkFH "\t\t\t\"domains\": [";
	for(my $j = 0; $j < scalar(@{$rankDomains[$i]}); $j++) {
		print $networkFH "\"$rankDomains[$i][$j]\"";
		if ($j < scalar(@{$rankDomains[$i]}) - 1) {
			print $networkFH ", ";
		}
	}
	print $networkFH "],\n";
	print $networkFH "\t\t\t\"motifs\": [";
	for(my $j = 0; $j < scalar(@{$rankMotifs[$i]}); $j++) {
		print $networkFH "\"$rankMotifs[$i][$j]\"";
		if ($j < scalar(@{$rankMotifs[$i]}) - 1) {
			print $networkFH ", ";
		}
	}
	print $networkFH "],\n";
	print $networkFH "\t\t\t\"diseases\": [";
	for(my $j = 0; $j < scalar(@{$rankDiseases[$i]}); $j++) {
		print $networkFH "\"$rankDiseases[$i][$j]\"";
		if ($j < scalar(@{$rankDiseases[$i]}) - 1) {
			print $networkFH ", ";
		}
	}
	if ($i < scalar(@rankNames) - 1) {
		print $networkFH "]\n\t\t},\n";
	} else {
		print $networkFH "]\n\t\t}\n";
	}
}
print $networkFH "\t]\n}\n";
close $networkFH ;


# output legend
print STDERR "Output legend\n";
open my $networkLegendFH, '>', $outputPrefix . '-network-legend.json';
print $networkLegendFH "{\n\t\"domains\": [\n";
for(my $i = 1; $i < scalar @rankNames; $i++) {
	print $networkLegendFH "\t\t{\n\t\t\t";
	print $networkLegendFH "\"domain\": $i,\n\t\t\t";
	print $networkLegendFH "\"color\": \"$colors[$i]\",\n\t\t\t";
	print $networkLegendFH "\"names\": [";
	for(my $j = 0; $j < scalar @{$rankNames[$i]}; $j++) {
		print $networkLegendFH "\"$rankNames[$i][$j]\"";
		if($j < scalar @{$rankNames[$i]} - 1) {
			print $networkLegendFH ", ";
		}
	}
	print $networkLegendFH "],\n\t\t\t";
  print $networkLegendFH "\"listedname\": [";
	for(my $j = 0; $j < scalar @{$rankListed[$i]}; $j++) {
		print $networkLegendFH "\"$rankListed[$i][$j]\"";
		if($j < scalar @{$rankListed[$i]} - 1) {
			print $networkLegendFH ", ";
		}
	}
	print $networkLegendFH "],\n\t\t\t";
	print $networkLegendFH "\"terms\": [";
	for(my $j = 0; $j < scalar @{$rankTerms[$i]}; $j++) {
		print $networkLegendFH "\"$rankTerms[$i][$j]\"";
		if($j < scalar @{$rankTerms[$i]} - 1) {
			print $networkLegendFH ", ";
		}
	}
	print $networkLegendFH "],\n\t\t\t";
	print $networkLegendFH "\"diseases\": [";
	for(my $j = 0; $j < scalar @{$rankDiseases[$i]}; $j++) {
		print $networkLegendFH "\"$rankDiseases[$i][$j]\"";
		if($j < scalar @{$rankDiseases[$i]} - 1) {
			print $networkLegendFH ", ";
		}
	}
	print $networkLegendFH "],\n\t\t\t";
	print $networkLegendFH "\"domains\": [";
	for(my $j = 0; $j < scalar @{$rankDomains[$i]}; $j++) {
		print $networkLegendFH "\"$rankDomains[$i][$j]\"";
		if($j < scalar @{$rankDomains[$i]} - 1) {
			print $networkLegendFH ", ";
		}
	}
	print $networkLegendFH "],\n\t\t\t";
	print $networkLegendFH "\"motifs\": [";
	for(my $j = 0; $j < scalar @{$rankMotifs[$i]}; $j++) {
		print $networkLegendFH "\"$rankMotifs[$i][$j]\"";
		if($j < scalar @{$rankMotifs[$i]} - 1) {
			print $networkLegendFH ", ";
		}
	}
	if($i < scalar(@rankNames) - 1) {
		print $networkLegendFH "]\n\t\t},\n";
	}
	else {
		print $networkLegendFH "]\n\t\t}\n";
	}
}
print $networkLegendFH "\t]\n}\n";
close $networkLegendFH ;
