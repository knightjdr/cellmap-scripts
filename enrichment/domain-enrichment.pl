#!/usr/bin/perl

# 13/6/2017

use strict;
use warnings;

# modules
use Data::Dumper; # use like this to print an array print Dumper \@array;
use List::MoreUtils qw(uniq);
use Spreadsheet::WriteExcel;
use String::Util qw(trim);
use Text::CSV_XS;
use Text::NSP::Measures::2D::Fisher::right;

# parameters
my $fdr = 0.01;
my $fileType = 'n';

# command line parameters
my $gfile = '';	# file with gene names and their associated domain or rank
my $nfile = ''; # file with NMF/SAFE details
my $ufile = '';	# file with list of domains for each gene

if ($#ARGV == 0) {
	print "\nTakes a file with domain information for each gene, and a list of NMF or SAFE\n";
	print "ranks/domains with all the genes in each and calculates domain enrichment.\n";
	print "\nusage:\n $0\n";
	print "-g [list of genes with ranks]\n";
	print "-n [file with NMF/SAFE details, for output. Optional. If absent ranks/domains just get numbered]\n";
	print "-u [list of domains for each gene]\n";
	print "-t [type of file analyzed. n = NMF (default), s = SAFE]\n";
	die "\n";
} else {
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-n'){
			$i++;
			$nfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-u'){
			$i++;
			$ufile = $ARGV[$i];
		} elsif($ARGV[$i] eq '-t'){
			$i++;
			$fileType = $ARGV[$i];
		} else{
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

# subs
sub bhCorrection {
	my @sortedPValues;
	my @sortedTerms;
	my %unsortedHash = %{$_[0]};
	my $fdr = $_[1];
	# Order p-values from smallest to largest.
	foreach my $term (sort { $unsortedHash{$a} <=> $unsortedHash{$b} } keys %unsortedHash) {
		push @sortedPValues, $unsortedHash{$term};
		push @sortedTerms, $term;
	}
	my $noTests = scalar @sortedPValues;
	# Set rank for each pvalue. Equivalent p-values get the same rank. The rank
	# only increments when the p-value changes.
	my $lastPvalue;
	my $rank = 1;
	my @ranks;
	for(my $i = 0; $i < $noTests; $i++) {
		if ($i == 0) {
			$ranks[0] = 1;
		} else {
			if ($sortedPValues[$i] > $lastPvalue) {
				$ranks[$i] = ++$rank;
			} else {
				$ranks[$i] = $rank;
			}
		}
		$lastPvalue = $sortedPValues[$i];
	}
	# I'm performing B-H correction two ways: 1) calculating the adjusted p-value
	# as p-value * (number of tests) / (p-value rank); and we want this to be less
	# than the specified FDR. 2) Using the B-H critical value which is the p-value's
	# rank * desired FDR / (number of tests); if the critical value >= p-value, that
	# is a false positive.
	my %adjustedPValues;
	my %correctedPValues;
	my $lastAdjusted = 1;
	for(my $i = $noTests - 1; $i >= 0; $i--) {
		my $adjustedP = $unsortedHash{$sortedTerms[$i]} * $noTests / $ranks[$i];
		if ($adjustedP > 1) {
			$adjustedP = 1;
		}
		if ($lastAdjusted < $adjustedP) {
			$adjustedP = $lastAdjusted;
		}
		$lastAdjusted = $adjustedP;
		$adjustedPValues{$sortedTerms[$i]} = $adjustedP;
		$correctedPValues{$sortedTerms[$i]} = $fdr * $ranks[$i] / $noTests;
	}
	return (\%adjustedPValues, \%correctedPValues);
}

sub fishers {
	# $a = genes that have the domain
	# $ab = genes tested in this set
	# $ac = genes with this domain in background
	# $n = total number of background genes
	my( $a, $ab, $ac, $n ) = @_;
	my $p = calculateStatistic(
		n11 => $a,
		n1p => $ab,
		np1 => $ac,
		npp => $n
	);
	return $p;
}

#get list of genes
print STDERR "Reading gene file\n";
open my $geneFH, '<', $gfile or die "Could not open $gfile: $!";
my $tsv = Text::CSV_XS->new({
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
});
$tsv->getline($geneFH); # discard header
my @geneList;
my %genesPerRank;
if ($fileType eq 'n') {
	while(my $row = $tsv->getline($geneFH)) {
		push @geneList, @{$row}[0];
		push @{$genesPerRank{@{$row}[1]}}, @{$row}[0];
	}
} else {
	while(my $row = $tsv->getline($geneFH)) {
		push @geneList, @{$row}[0];
		push @{$genesPerRank{@{$row}[1]}}, @{$row}[0];
	}
}
close $geneFH;

#get domain information
print STDERR "Reading domain info \n";
open my $domainFH, '<', $ufile or die "Could not open $ufile: $!";
$tsv = Text::CSV_XS->new({
	sep_char => "\t",
});
$tsv->getline($domainFH); # discard header
my %completeDomainsPerGene;
while(my $row = $tsv->getline($domainFH)) {
	my @currDomains = split ';', @{$row}[2];
	@{$completeDomainsPerGene{@{$row}[0]}{'domains'}} = uniq @currDomains;
	$completeDomainsPerGene{@{$row}[0]}{'uniprot'} = @{$row}[1];
}
close $domainFH;

# remove unneeded genes
my %domainList;
my %domainsPerGene;
my $totalGeneNumber = 0;
for my $gene (@geneList) {
	if (exists $completeDomainsPerGene{$gene}) {
		$totalGeneNumber++;
		@{$domainsPerGene{$gene}{'domains'}} = @{$completeDomainsPerGene{$gene}{'domains'}};
		# $domainsPerGene{$gene}{'transmembrane'} = $completeDomainsPerGene{$gene}{'transmembrane'};
		$domainsPerGene{$gene}{'uniprot'} = $completeDomainsPerGene{$gene}{'uniprot'};
		for my $domain (@{$domainsPerGene{$gene}{'domains'}}) {
			push @{$domainList{$domain}}, $gene;
		}
		# if ($domainsPerGene{$gene}{'transmembrane'}) {
		#	push @{$domainList{'transmembrane'}}, $gene;
		# }
	}
}

# read rank/domain details
my @rankMap;
if ($nfile) {
	print STDERR "Reading rank info \n";
	open my $rankFH, '<', $nfile or die "Could not open $nfile: $!";
	$tsv = Text::CSV_XS->new({
		sep_char => "\t",
	});
	$tsv->getline($rankFH); # discard header
	while(my $row = $tsv->getline($rankFH)) {
		my $termString = @{$row}[1];
		$termString =~ s/[\[\]"]//g;
		$rankMap[@{$row}[0]] = $termString;
	}
	close $rankFH;
}

# count domains per rank
open my $outputFH, '>', 'domain.txt';
print $outputFH "rank\tterm\tmatched\tbackground_size\tfold enrichment\tpvalue\tadj. pvalue\tbhfdr\tgenes\n";
my $workbook = Spreadsheet::WriteExcel->new('domain.xls');
my $worksheet = $workbook->add_worksheet('details');
$worksheet->write_string(0, 0, 'Input list size');
$worksheet->write_number(0, 1, scalar @geneList);
$worksheet->write_string(1, 0, 'Genes with information available');
$worksheet->write_number(1, 1, $totalGeneNumber);
for my $rank ( sort {$a <=> $b} keys %genesPerRank) {
	my %rankDetails;
	my %genesInTerm;
	my $noRankGenes = 0;
	for my $gene (@{$genesPerRank{$rank}}) {
		if (exists $domainsPerGene{$gene}) {
			$noRankGenes++; # this is here because we don't won't to count genes with no information
			# if ($domainsPerGene{$gene}{'transmembrane'}) {
			#	if (exists $rankDetails{'transmembrane'}) {
			#		$rankDetails{'transmembrane'}++;
			#	} else {
			#		$rankDetails{'transmembrane'} = 1;
			#	}
			#	push @{$genesInTerm{'transmembrane'}}, $gene;
			# }
			for my $geneDomain (@{$domainsPerGene{$gene}{'domains'}}) {
				if (exists $rankDetails{$geneDomain}) {
					$rankDetails{$geneDomain}++;
				} else {
					$rankDetails{$geneDomain} = 1;
				}
				push @{$genesInTerm{$geneDomain}}, $gene;
			}
		}
	}
	my %correctionArray;
	my %rankOutput;
	foreach my $term (keys %rankDetails) {
		$rankOutput{$term}{'genes'} = join ',', @{$genesInTerm{$term}};
		$rankOutput{$term}{'matched'} = $rankDetails{$term};
		$rankOutput{$term}{'totalGenes'} = $noRankGenes;
		my $totalSet = scalar @{$domainList{$term}};
		$rankOutput{$term}{'backgroundSet'} = scalar $totalSet;
		$rankOutput{$term}{'fEnrichment'} = sprintf '%.3f', ($rankDetails{$term} / $noRankGenes) / ($totalSet / $totalGeneNumber);
		$rankOutput{$term}{'pValue'} = fishers($rankDetails{$term}, $noRankGenes, $totalSet, $totalGeneNumber);
		$correctionArray{$term} = $rankOutput{$term}{'pValue'};
	}
	my ($adjustedPValueRef, $correctedFDRRef) = bhCorrection(\%correctionArray, $fdr);
	my %adjustedPValue = %{$adjustedPValueRef};
	my %correctedFDR = %{$correctedFDRRef};
	my @sortedRankOutput;
	my $sortedOrder = 0;
	foreach my $term (sort { $correctedFDR{$a} <=> $correctedFDR{$b} } keys %correctedFDR) {
		$sortedRankOutput[$sortedOrder++] = $term;
	}
	$worksheet = $workbook->add_worksheet($rank);
	$worksheet->set_column( 'A:A' , 25 );
	$worksheet->set_column( 'B:G' , 15 );
	$worksheet->write_string(0, 0, 'term');
	$worksheet->write_string(0, 1, 'matched');
	$worksheet->write_string(0, 2, 'genes in rank (with annotations)');
	$worksheet->write_string(0, 3, 'background term size');
	$worksheet->write_string(0, 4, 'fold enrichment');
	$worksheet->write_string(0, 5, 'p-value');
	$worksheet->write_string(0, 6, 'adj. p-value');
	$worksheet->write_string(0, 7, 'B-H FDR');
	$worksheet->write_string(0, 8, 'genes');
	if ($nfile) {
		$worksheet->write_string(0, 9, 'rank name: ' . $rankMap[$rank]);
	} else {
		$worksheet->write_string(0, 9, 'rank name: ' . $rank);
	}
	my $row = 1;
	for(my $i = 0; $i < scalar @sortedRankOutput; $i++) {
		if ($rankOutput{$sortedRankOutput[$i]}{'pValue'} == 0 ||
			$rankOutput{$sortedRankOutput[$i]}{'pValue'} < $correctedFDR{$sortedRankOutput[$i]}) {
			print $outputFH "$rank\t";
			print $outputFH "$sortedRankOutput[$i]\t";
			print $outputFH "$rankOutput{$sortedRankOutput[$i]}{'matched'}\t";
			print $outputFH "$rankOutput{$sortedRankOutput[$i]}{'backgroundSet'}\t";
			print $outputFH "$rankOutput{$sortedRankOutput[$i]}{'fEnrichment'}\t";
			print $outputFH "$rankOutput{$sortedRankOutput[$i]}{'pValue'}\t";
			print $outputFH "$adjustedPValue{$sortedRankOutput[$i]}\t";
			print $outputFH "$correctedFDR{$sortedRankOutput[$i]}\t";
			print $outputFH "$rankOutput{$sortedRankOutput[$i]}{'genes'}\n";
			$worksheet->write_string($row, 0, $sortedRankOutput[$i]);
			$worksheet->write_number($row, 1, $rankOutput{$sortedRankOutput[$i]}{'matched'});
			$worksheet->write_number($row, 2, $rankOutput{$sortedRankOutput[$i]}{'totalGenes'});
			$worksheet->write_number($row, 3, $rankOutput{$sortedRankOutput[$i]}{'backgroundSet'});
			$worksheet->write_number($row, 4, $rankOutput{$sortedRankOutput[$i]}{'fEnrichment'});
			$worksheet->write_number($row, 5, $rankOutput{$sortedRankOutput[$i]}{'pValue'});
			$worksheet->write_number($row, 6, $adjustedPValue{$sortedRankOutput[$i]});
			$worksheet->write_number($row, 7, $correctedFDR{$sortedRankOutput[$i]});
			$worksheet->write_string($row, 8, $rankOutput{$sortedRankOutput[$i]}{'genes'});
			$row++;
		}
	}
}
$workbook->close();
