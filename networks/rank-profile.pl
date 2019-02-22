#!/usr/bin/perl

# 23/6/2017

use strict;
use warnings;

# libraries
use Text::CSV_XS;

# parameters
my $fileType = 'n';
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
my $dfile = ''; # domains associated with each rank
my $disfile = ''; # diseases associated with each rank
my $gfile = ''; # GO terms associated with each rank
my $motiffile = ''; # motifs associated with each rank
my $ofile = ''; # GO file mapping names to IDs
my $rfile = ''; # for NMF and SAFE, file with rank/domain details
my $sfile = ''; # order for NMF ranks

if ($#ARGV==0) {
	print "\nTakes information about each NMF rank or SAFE domain and creates a json\n";
	print "for use with the Cellmap.\n";
	print "\nusage:\n $0\n";
	print "-domains [domains associated with each rank]\n";
	print "-diseases [diseases associated with each rank (optional)]\n";
	print "-motifs [motifs associated with each rank (optional)]\n";
  print "-g [GO terms associated with each rank]\n";
	print "-o [GO file for mapping names to ids, need for SAFE data]\n";
	print "-r [NMF or SAFE file with rank/domain details]\n";
	print "-s [file with the NMF rank order]\n";
	print "-t [n for NMF (default) or s for SAFE]\n";
	die "\n";
} else {
	my $i = 0;
	while($i <= $#ARGV) {
		if ($ARGV[$i] eq '-domains') {
			$i++;
			$dfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-diseases') {
			$i++;
			$disfile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-motifs') {
			$i++;
			$motiffile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-g') {
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-o') {
			$i++;
			$ofile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-r') {
			$i++;
			$rfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-s') {
			$i++;
			$sfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t') {
			$i++;
			$fileType = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# get domains associated with each rank
print STDERR "Reading domains\n";
open my $domainFH, '<', $dfile or die "Could not open $dfile: $!";
my $domainTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$domainTSV->getline($domainFH); # discard header
my @domains;
while(my $row = $domainTSV->getline($domainFH)) {
	my %currHash = (
		'genes' => @{$row}[8],
		'fold' => @{$row}[4],
		'matched' => @{$row}[2],
		'pvalue' => @{$row}[6],
		'term' => @{$row}[1],
	);
	push @{$domains[@{$row}[0]]}, \%currHash;
}
close $domainFH;

# reading disease information
my @diseases;
if ($disfile) {
	print STDERR "Reading diseases\n";
	open my $diseaseFH, '<', $disfile or die "Could not open $disfile: $!";
	my $diseaseTSV = Text::CSV_XS->new({
		sep_char => "\t",
	});
	$diseaseTSV->getline($diseaseFH); # discard header
	while(my $row = $diseaseTSV->getline($diseaseFH)) {
		my %currHash = (
			'genes' => @{$row}[8],
			'fold' => @{$row}[4],
			'matched' => @{$row}[2],
			'pvalue' => @{$row}[6],
			'term' => @{$row}[1],
		);
		push @{$diseases[@{$row}[0]]}, \%currHash;
	}
	close $diseaseFH;
}

# reading motif information
my @motifs;
if ($motiffile) {
	print STDERR "Reading motifs\n";
	open my $motifFH, '<', $motiffile or die "Could not open $motiffile: $!";
	my $motifTSV = Text::CSV_XS->new({
		sep_char => "\t",
	});
	$motifTSV->getline($motifFH); # discard header
	while(my $row = $motifTSV->getline($motifFH)) {
		my %currHash = (
			'genes' => @{$row}[8],
			'fold' => @{$row}[4],
			'matched' => @{$row}[2],
			'pvalue' => @{$row}[6],
			'term' => @{$row}[1],
		);
		push @{$motifs[@{$row}[0]]}, \%currHash;
	}
	close $motifFH;
}

# read in a map of GO terms to IDs if SAFE
my %goMap;
if ($fileType eq 's') {
	print STDERR "Creating GO map\n";
	open my $gomapFH, '<', $ofile or die "Could not open $ofile: $!";
	my $gomapTSV = Text::CSV_XS->new({
		sep_char => "\t",
	});
	while(my $row = $gomapTSV->getline($gomapFH)) {
		$goMap{@{$row}[1]} = @{$row}[0];
	}
	close $gomapFH;
}

# get GO terms associated with each rank
print STDERR "Reading GO terms\n";
open my $goFH, '<', $gfile or die "Could not open $gfile: $!";
my $goTSV = Text::CSV_XS->new({
	sep_char => "\t",
});
$goTSV->getline($goFH); # discard header
my @goTerms;
while(my $row = $goTSV->getline($goFH)) {
	my %currHash;
	if ($fileType eq 'n') {
		%currHash = (
			'genes' => @{$row}[5],
			'id' => @{$row}[4],
			'matched' => @{$row}[2],
			'name' => @{$row}[1],
			'pvalue' => @{$row}[3],
		);
		push @{$goTerms[@{$row}[0]]}, \%currHash;
	} else {
		%currHash = (
			'id' => $goMap{@{$row}[1]},
			'name' => @{$row}[1],
		);
		push @{$goTerms[@{$row}[2]]}, \%currHash;
	}
}
close $goFH;

# read rank information
print STDERR "Reading rank information\n";
open my $rankFH, '<', $rfile or die "Could not open $rfile: $!";
my $rankTSV = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
$rankTSV->getline($rankFH); # discard header
my @rankInfo;
while(my $row = $rankTSV->getline($rankFH)) {
	my $idString = @{$row}[3];
	$idString =~ s/[\[\]"]//g;
	my $listedString = @{$row}[2];
	$listedString =~ s/[\[\]"]//g;
	my $termString = @{$row}[1];
	$termString =~ s/[\[\]"]//g;
	@{$rankInfo[@{$row}[0]]{'id'}} = split /, /, $idString;
	@{$rankInfo[@{$row}[0]]{'listedname'}} = split /, /, $listedString;
	@{$rankInfo[@{$row}[0]]{'term'}} = split /, /, $termString;
}
close $rankFH;

# read NMF rank order
my @rankOrder;
if ($fileType eq 'n') {
	print STDERR "Reading rank order\n";
	open my $orderFH, '<', $sfile or die "Could not open $sfile: s!";
	while(<$orderFH>) {
		my @rankArray = split /, /, $_;
		for(my $i = 0; $i < scalar @rankArray; $i++) {
			$rankOrder[$rankArray[$i]] = $i + 1;
		}
	}
	close $orderFH;
}

# output
my $categoryType = 'rank';
my $outfile = 'ranks.json';
if ($fileType eq 's') {
	$categoryType = 'domain';
	$outfile = 'domains.json';
}
open my $outputFH, '>', $outfile;
print $outputFH "[\n";
for(my $i = 1; $i < scalar @rankInfo; $i++) {
	my $rank = $i;
	print $outputFH "\t{\n";
	print $outputFH "\t\t\"$categoryType\": $rank,\n";
	print $outputFH "\t\t\"color\": \"$colors[$i]\",\n";
	print $outputFH "\t\t\"names\": [";
	for (my $j = 0; $j < scalar @{$rankInfo[$i]{'term'}}; $j++) {
		print $outputFH "\"$rankInfo[$i]{'term'}[$j]\"";
		if ($j < scalar @{$rankInfo[$i]{'term'}} - 1) {
			print $outputFH ", ";
		}
	}
	print $outputFH "],\n";
	print $outputFH "\t\t\"listedname\": [";
	for (my $j = 0; $j < scalar @{$rankInfo[$i]{'listedname'}}; $j++) {
		print $outputFH "\"$rankInfo[$i]{'listedname'}[$j]\"";
		if ($j < scalar @{$rankInfo[$i]{'listedname'}} - 1) {
			print $outputFH ", ";
		}
	}
	print $outputFH "],\n";
	print $outputFH "\t\t\"go_id\": [";
	for (my $j = 0; $j < scalar @{$rankInfo[$i]{'id'}}; $j++) {
		print $outputFH "\"$rankInfo[$i]{'id'}[$j]\"";
		if ($j < scalar @{$rankInfo[$i]{'id'}} - 1) {
			print $outputFH ", ";
		}
	}
	print $outputFH "],\n";
	if ($fileType eq 'n') {
		print $outputFH "\t\t\"order\": $rankOrder[$i],\n";
	}
	if ($disfile) {
		print $outputFH "\t\t\"diseases\": [\n";
		if ($diseases[$rank]) {
			for (my $j = 0; $j < scalar @{$diseases[$rank]}; $j++) {
				my %currHash = %{$diseases[$rank][$j]};
				print $outputFH "\t\t\t{\n";
				print $outputFH "\t\t\t\t\"disease_name\": \"$currHash{'term'}\",\n";
				print $outputFH "\t\t\t\t\"disease_foldenrich\": $currHash{'fold'},\n";
				print $outputFH "\t\t\t\t\"disease_pvalue\": $currHash{'pvalue'},\n";
				print $outputFH "\t\t\t\t\"disease_matched\": $currHash{'matched'},\n";
				print $outputFH "\t\t\t\t\"disease_genes\": \"$currHash{'genes'}\"\n";
				print $outputFH "\t\t\t}";
				if ($j < scalar @{$diseases[$rank]} - 1) {
					print $outputFH ",";
				}
				print $outputFH "\n";
			}
		}
		print $outputFH "\t\t],\n";
	}
	if ($dfile) {
		print $outputFH "\t\t\"domains\": [\n";
		if ($domains[$rank]) {
			for (my $j = 0; $j < scalar @{$domains[$rank]}; $j++) {
				my %currHash = %{$domains[$rank][$j]};
				print $outputFH "\t\t\t{\n";
				print $outputFH "\t\t\t\t\"domain_name\": \"$currHash{'term'}\",\n";
				print $outputFH "\t\t\t\t\"domain_foldenrich\": $currHash{'fold'},\n";
				print $outputFH "\t\t\t\t\"domain_pvalue\": $currHash{'pvalue'},\n";
				print $outputFH "\t\t\t\t\"domain_matched\": $currHash{'matched'},\n";
				print $outputFH "\t\t\t\t\"domain_genes\": \"$currHash{'genes'}\"\n";
				print $outputFH "\t\t\t}";
				if ($j < scalar @{$domains[$rank]} - 1) {
					print $outputFH ",";
				}
				print $outputFH "\n";
			}
		}
		print $outputFH "\t\t],\n";
	}
	if ($motiffile) {
		print $outputFH "\t\t\"motifs\": [\n";
		if ($motifs[$rank]) {
			for (my $j = 0; $j < scalar @{$motifs[$rank]}; $j++) {
				my %currHash = %{$motifs[$rank][$j]};
				print $outputFH "\t\t\t{\n";
				print $outputFH "\t\t\t\t\"motif_name\": \"$currHash{'term'}\",\n";
				print $outputFH "\t\t\t\t\"motif_foldenrich\": $currHash{'fold'},\n";
				print $outputFH "\t\t\t\t\"motif_pvalue\": $currHash{'pvalue'},\n";
				print $outputFH "\t\t\t\t\"motif_matched\": $currHash{'matched'},\n";
				print $outputFH "\t\t\t\t\"motif_genes\": \"$currHash{'genes'}\"\n";
				print $outputFH "\t\t\t}";
				if ($j < scalar @{$motifs[$rank]} - 1) {
					print $outputFH ",";
				}
				print $outputFH "\n";
			}
		}
		print $outputFH "\t\t],\n";
	}
	print $outputFH "\t\t\"terms\": [\n";
	if ($goTerms[$rank]) {
		for (my $j = 0; $j < scalar @{$goTerms[$rank]}; $j++) {
			my %currHash = %{$goTerms[$rank][$j]};
			print $outputFH "\t\t\t{\n";
			print $outputFH "\t\t\t\t\"go_name\": \"$currHash{'name'}\",\n";
			print $outputFH "\t\t\t\t\"go_id\": \"$currHash{'id'}\"";
			if ($fileType eq 'n') {
				print $outputFH ",\n";
				print $outputFH "\t\t\t\t\"go_pvalue\": $currHash{'pvalue'},\n";
				print $outputFH "\t\t\t\t\"go_matched\": $currHash{'matched'},\n";
				print $outputFH "\t\t\t\t\"go_genes\": \"$currHash{'genes'}\"\n";
			} else {
				print $outputFH "\n";
			}
			print $outputFH "\t\t\t}";
			if ($j < scalar @{$goTerms[$rank]} - 1) {
				print $outputFH ",";
			}
			print $outputFH "\n";
		}
	}
	print $outputFH "\t\t]\n";
	print $outputFH "\t}";
	if ($i < scalar @rankInfo - 1) {
		print $outputFH ",";
	}
	print $outputFH "\n";
 }
 print $outputFH "]";
 close $outputFH;
