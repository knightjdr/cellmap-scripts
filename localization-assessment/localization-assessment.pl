#!/usr/bin/perl

# 2/6/2017

use strict;
use warnings;

# libraries
use List::MoreUtils qw(uniq);
use Statistics::Basic qw(mean median);
use String::Util qw(trim);
use Text::CSV_XS;

# parameters
my $fileType = 'n';
my $namespace = 'cc';
my @icCutoffs = (2, 1.7, 1, 0.6, 0);
my %tsvParams = (
	binary => 1,
	sep_char => "\t",
	quote_char => undef,
	escape_char => undef,
	allow_loose_quotes => 1,
	allow_loose_escapes => 1,
);

# command line parameters
my $hfile = ''; # hierarchy file
my $gfile = ''; # tsv file (with header), two columns with gene and it's rank/domain/GO term
my $ifile = ''; # tsv file with information content
my $lfile = ''; # GO localization file (goa_human.gaf)
my $oPrefix = '';
my $rfile = ''; # for NMF and SAFE, file with rank/domain details

if ($#ARGV==0){
	print "\nTakes GO hierarchy from get_children.pl, a list of official localizations for genes,\n";
	print "assigned localization form NMF, SAFE or an external source, and assesses the\n";
	print "assigned localizations.\n";
	print "\nusage:\n $0\n";
	print "-g [gene localization file]\n";
  print "-h [hierarcy file procduced by get_children.pl]\n";
	print "-i [file with information content for each GO term]\n";
	print "-l [GO localization file (goa_human.gaf)]\n";
  print "-o [output-file prefix]\n";
	print "-r [when assessing NMF or SAFE, file with rank/domain details]\n";
	print "-t [type of assessment analyzed. n = NMF (default), h = HPA, o = other, s = SAFE]\n";
	die "\n";
} else{
	my $i = 0;
	while($i<=$#ARGV){
		if ($ARGV[$i] eq '-g'){
			$i++;
			$gfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-h'){
			$i++;
			$hfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-i'){
			$i++;
			$ifile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-l'){
			$i++;
			$lfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-o'){
      $i++;
      $oPrefix = $ARGV[$i];
    } elsif ($ARGV[$i] eq '-r'){
			$i++;
			$rfile = $ARGV[$i];
		} elsif ($ARGV[$i] eq '-t'){
			$i++;
			$fileType = $ARGV[$i];
		} else {
			die "\nIncorrect program usage\n\n";
		}
		$i++;
	}
}

sub assessLocalization {
  my @possibleTerms = @{$_[0]};
  my %ic = %{$_[1]};
  my %children = %{$_[2]};
  my @geneLocalizations = @{$_[3]};

  my $icWorst;
  my $match = 0;
  my $noInitialTerms = scalar @possibleTerms;
  foreach my $term (@possibleTerms) {
    if (exists $ic{$term}) {
      my $currIC = $ic{$term};
      if (!$icWorst) {
        $icWorst = $currIC;
      } elsif ($currIC < $icWorst) {
        $icWorst = $currIC;
      }
    }
    if (exists $children{$term}) {
      push @possibleTerms, @{$children{$term}};
    }
  }
  @possibleTerms = uniq @possibleTerms;
  if (scalar @possibleTerms > 0 && scalar @geneLocalizations > 0) {
    foreach my $term (@geneLocalizations) {
      if (grep /^$term$/, @possibleTerms) {
        $match = 1;
        last;
      }
    }
  }
  return ($icWorst, $match, $noInitialTerms);
}

sub assessTier {
  my $icWorst = $_[0];
  my $icCutoffs = @{$_[1]};

  my $tier;

  if ($icWorst && $icWorst >= $icCutoffs[0]) {
    $tier = 'tier1';
  } elsif ($icWorst && $icWorst >= $icCutoffs[1]) {
    $tier = 'tier2';
  } elsif ($icWorst && $icWorst >= $icCutoffs[2]) {
    $tier = 'tier3';
  } elsif ($icWorst && $icWorst >= $icCutoffs[3]) {
    $tier = 'tier4';
  } elsif ($icWorst) {
    $tier = 'tier5';
  } else {
    $tier = 'unknown';
  }
  return $tier;
}

sub outputSummary {
  my $file = $_[0];
  my %correctAssignments = %{$_[1]};
  my %totalGenes = %{$_[2]};
  my $noLocalizations = $_[3];
  my %localizationAssessment = %{$_[4]};
  my %rankToGo = %{$_[5]};
  my $type = $_[6];

  
  my $tier1Frac = 0;
  my $tier2Frac = 0;
  my $tier3Frac = 0;
  my $tier4Frac = 0;
  my $tier5Frac = 0;
  if ($totalGenes{'tier1'}) {
    $tier1Frac = sprintf "%.3f", $correctAssignments{'tier1'} / $totalGenes{'tier1'};
  }
  if ($totalGenes{'tier2'}) {
    $tier2Frac = sprintf "%.3f", $correctAssignments{'tier2'} / $totalGenes{'tier2'};
  }
  if ($totalGenes{'tier3'}) {
    $tier3Frac = sprintf "%.3f", $correctAssignments{'tier3'} / $totalGenes{'tier3'};
  }
  if ($totalGenes{'tier4'}) {
    $tier4Frac = sprintf "%.3f", $correctAssignments{'tier4'} / $totalGenes{'tier4'};
  }
  if ($totalGenes{'tier5'}) {
    $tier5Frac = sprintf "%.3f", $correctAssignments{'tier5'} / $totalGenes{'tier5'};
  }

  if ($type ne 's') {
    print $file "Summary\n";
    my $totalFrac = sprintf "%.3f", $correctAssignments{'total'} / $totalGenes{'total'};
    print $file "total\t$totalFrac\t$totalGenes{'total'}\t$correctAssignments{'total'}\n\n";
  }

  print $file "Tiers\n";
  print $file "tier1\t$tier1Frac\t$totalGenes{'tier1'}\t$correctAssignments{'tier1'}\n";
  print $file "tier2\t$tier2Frac\t$totalGenes{'tier2'}\t$correctAssignments{'tier2'}\n";
  print $file "tier3\t$tier3Frac\t$totalGenes{'tier3'}\t$correctAssignments{'tier3'}\n";
  print $file "tier4\t$tier4Frac\t$totalGenes{'tier4'}\t$correctAssignments{'tier4'}\n";
  print $file "tier5\t$tier5Frac\t$totalGenes{'tier5'}\t$correctAssignments{'tier5'}\n";
  print $file "unknown\t0\t$totalGenes{'unknown'}\t0\n\n";

  print $file "Per localization\n";
  for(my $i = 1; $i <= $noLocalizations; $i++) {
    my $joinedTerms;
    if ($type eq 'n') {
      $joinedTerms = join ', ', @{$rankToGo{'displayname'}[$i]};
    } else {
      $joinedTerms = $rankToGo{'index'}{$i};
    }
    my $locFrac = sprintf "%.3f", $localizationAssessment{'known'}[$i] / $localizationAssessment{'total'}[$i];
    print $file "rank $i, $joinedTerms\t$locFrac\t$localizationAssessment{'total'}[$i]\t$localizationAssessment{'known'}[$i]\n";
  }
}

# set namespace
my $shortNameSpace;
my $longNameSpace;
if($namespace eq "bp") {
	$shortNameSpace = "P";
	$longNameSpace = "biological_process";
} elsif($namespace eq "cc") {
	$shortNameSpace = "C";
	$longNameSpace = "cellular_component";
} elsif($namespace eq "mf") {
	$shortNameSpace = "F";
	$longNameSpace = "molecular_function";
} else {
}

# get child terms for each GO term
print STDERR "Reading hierarchy file\n";
open my $hierarchyFH, '<', $hfile or die "Could not open $hfile: $!";
my $tsv = Text::CSV_XS->new(\%tsvParams);
my %children;
while(my $row = $tsv->getline($hierarchyFH)) {
	@{$children{@{$row}[0]}} = split ',', @{$row}[1];
}
close $hierarchyFH;

# read information content
$tsv = Text::CSV_XS->new(\%tsvParams);
open my $icfh, '<', $ifile or die "Could not open $ifile: $!";
$tsv->getline($icfh); # discard header
my %informationContent;
while(my $row = $tsv->getline($icfh)) {
	$informationContent{@{$row}[0]} = @{$row}[2];
}
close $icfh;

# reading known localizations
print STDERR "Reading known localizations\n";
$tsv = Text::CSV_XS->new(\%tsvParams);
open my $localizationfh, '<', $lfile or die "Could not open $lfile: $!";
my %geneLocalizations;
while(my $row = $tsv->getline($localizationfh)) {
	if (@{$row}[8] eq $shortNameSpace) {
		push @{$geneLocalizations{@{$row}[2]}}, @{$row}[4];
	}
}
close $localizationfh;

# get rank details
my $noRanks = 0;
my %rankToGo;
if ($fileType eq 'n' || $fileType eq 's') {
	print STDERR "Reading rank details\n";
	my $tsv = Text::CSV_XS->new({
		binary => 1,
		sep_char => "\t",
		escape_char => undef,
		allow_loose_escapes => 1,
	});
	open my $rankfh, '<', $rfile or die "Could not open $rfile: $!";
	$tsv->getline($rankfh); # discard header
	while(my $row = $tsv->getline($rankfh)) {
		$noRanks++;
		my $idString = @{$row}[3];
		$idString =~ s/[\[\]"]//g;
		my $listedString = @{$row}[2];
		$listedString =~ s/[\[\]"]//g;
		my $termString = @{$row}[1];
		$termString =~ s/[\[\]"]//g;
		@{$rankToGo{'term'}[@{$row}[0]]} = split /, /, $termString;
		@{$rankToGo{'displayname'}[@{$row}[0]]} = split /, /, $listedString;
		@{$rankToGo{'go'}[@{$row}[0]]} = split /, /, $idString;
	}
	close $rankfh;
}

# get gene localizations
print STDERR "Reading assigned gene localizations\n";
my %actualLocalizationName;
my %assignedLocalization;
my %indexMap;
open my $assignedlocalizationfh, '<', $gfile or die "Could not open $gfile: $!";
if ($fileType eq 'n') {
	$tsv = Text::CSV_XS->new(\%tsvParams);
	$tsv->getline($assignedlocalizationfh); # discard header
	while(my $row = $tsv->getline($assignedlocalizationfh)) {
		$assignedLocalization{@{$row}[0]} = @{$row}[1];
	}
} elsif ($fileType eq 's') {
	$tsv = Text::CSV_XS->new(\%tsvParams);
	$tsv->getline($assignedlocalizationfh); # discard header
	while(my $row = $tsv->getline($assignedlocalizationfh)) {
		$assignedLocalization{@{$row}[1]} = @{$row}[2];
	}
} elsif ($fileType eq 'h') {
	$tsv = Text::CSV_XS->new({
		binary => 1,
		empty_is_undef => 1,
		escape_char => undef,
		allow_loose_escapes => 1,
		sep_char => "\t",
	});
	$tsv->getline($assignedlocalizationfh); # discard header
	while(my $row = $tsv->getline($assignedlocalizationfh)) {
		my @currGoList = split /;/, @{$row}[10];
		my %currTermMap;
		for(my $i = 0; $i < scalar @currGoList; $i++) {
			my ($currTerm, $currGoTerm) = $currGoList[$i] =~ /([^\(]+)\((GO:\d+)?/;
			$currTerm = trim $currTerm;
			if ($currGoTerm) {
				$currTermMap{$currTerm} = $currGoTerm;
			} else {
				$currTermMap{$currTerm} = 'GO:XXXXXXX';
			}
		}
		if (@{$row}[3]) {
			my @currTermArray = split /;/, @{$row}[3];
			for(my $i = 0; $i < scalar @currTermArray; $i++) {
				push @{$assignedLocalization{@{$row}[1]}{'enhanced'}}, $currTermMap{$currTermArray[$i]};
				push @{$assignedLocalization{@{$row}[1]}{'supported'}}, $currTermMap{$currTermArray[$i]};
				push @{$assignedLocalization{@{$row}[1]}{'approved'}}, $currTermMap{$currTermArray[$i]};
				push @{$assignedLocalization{@{$row}[1]}{'uncertain'}}, $currTermMap{$currTermArray[$i]};
			}
		}
		if (@{$row}[4]) {
			my @currTermArray = split /;/, @{$row}[4];
			for(my $i = 0; $i < scalar @currTermArray; $i++) {
				push @{$assignedLocalization{@{$row}[1]}{'supported'}}, $currTermMap{$currTermArray[$i]};
				push @{$assignedLocalization{@{$row}[1]}{'approved'}}, $currTermMap{$currTermArray[$i]};
				push @{$assignedLocalization{@{$row}[1]}{'uncertain'}}, $currTermMap{$currTermArray[$i]};
			}
		}
		if (@{$row}[5]) {
			my @currTermArray = split /;/, @{$row}[5];
			for(my $i = 0; $i < scalar @currTermArray; $i++) {
				push @{$assignedLocalization{@{$row}[1]}{'approved'}}, $currTermMap{$currTermArray[$i]};
				push @{$assignedLocalization{@{$row}[1]}{'uncertain'}}, $currTermMap{$currTermArray[$i]};
			}
		}
		if (@{$row}[6]) {
			my @currTermArray = split /;/, @{$row}[6];
			for(my $i = 0; $i < scalar @currTermArray; $i++) {
				push @{$assignedLocalization{@{$row}[1]}{'uncertain'}}, $currTermMap{$currTermArray[$i]};
			}
		}
	}
} elsif ($fileType eq 'o') {
	$tsv = Text::CSV_XS->new(\%tsvParams);
	$tsv->getline($assignedlocalizationfh); # discard header
	while(my $row = $tsv->getline($assignedlocalizationfh)) {
		$assignedLocalization{@{$row}[0]} = @{$row}[2];
		if (!exists $indexMap{'id'}) {
			++$noRanks;
			$indexMap{'id'}{@{$row}[2]} = $noRanks;
			$indexMap{'index'}{$noRanks} = @{$row}[1];
		} elsif (!exists $indexMap{'id'}{@{$row}[2]}) {
			++$noRanks;
			$indexMap{'id'}{@{$row}[2]} = $noRanks;
			$indexMap{'index'}{$noRanks} = @{$row}[1];
		}
	}
}
close $assignedlocalizationfh;

# assess localization
print STDERR "Assessing assigned localizations\n";
my %correctAssignments = (
	'total' => 0,
	'tier1' => 0,
	'tier2' => 0,
	'tier3' => 0,
	'tier4' => 0,
  'tier5' => 0,
);
my @ic;
my %localizationAssessment;
for(my $i = 1; $i <= $noRanks; $i++) {
	$localizationAssessment{'known'}[$i] = 0;
	$localizationAssessment{'total'}[$i] = 0;
}
my @localizationsPerGene;
my %totalGenes = (
	'total' => 0,
	'tier1' => 0,
	'tier2' => 0,
	'tier3' => 0,
	'tier4' => 0,
  'tier5' => 0,
	'unknown' => 0,
);
open my $outfile, '>', $oPrefix . 'localization-assessment.txt';
if ($fileType eq 'n' || $fileType eq 's') {
	open my $nodeAttributeFH, '>', 'node_attributes_cytoscape.txt';
	print $nodeAttributeFH "name\trank (known)\n";
	while( my( $gene, $rank ) = each %assignedLocalization ) {
		$localizationAssessment{'total'}[$rank]++;
		$totalGenes{'total'}++;
    my @possibleTerms = @{$rankToGo{'go'}[$rank]};
    my ($icWorst, $match, $noInitialTerms) = assessLocalization(\@possibleTerms, \%informationContent, \%children, \@{$geneLocalizations{$gene}});
		$correctAssignments{'total'} += $match;
		$localizationAssessment{'known'}[$rank] += $match;
    if ($noInitialTerms > 0) {
      push @localizationsPerGene, $noInitialTerms;
    }

		# print out a file for cytoscape, listing gene and rank (if known, otherwise 0)
		if ($match == 1) {
			print $nodeAttributeFH "$gene\t$rank\n";
		} else {
			print $nodeAttributeFH "$gene\t0\n";
		}
    my $tier = assessTier($icWorst, \@icCutoffs);
		$totalGenes{$tier}++;
		if ($match == 1){
			$correctAssignments{$tier}++;
		}
	}
  close $nodeAttributeFH;
  outputSummary($outfile, \%correctAssignments, \%totalGenes, $noRanks, \%localizationAssessment, \%rankToGo, 'n');
} elsif ($fileType eq 'h') {
	my %correct;
	my %total;
	foreach my $gene (keys %assignedLocalization) {
		if ($assignedLocalization{$gene}{'enhanced'}) {
      $total{'enhanced'}++;
      my ($icWorst, $match, $noInitialTerms) = assessLocalization(\@{$assignedLocalization{$gene}{'enhanced'}}, \%informationContent, \%children, \@{$geneLocalizations{$gene}});
      $correct{'enhanced'} += $match;
		}
		if ($assignedLocalization{$gene}{'supported'}) {
      $total{'supported'}++;
      my ($icWorst, $match, $noInitialTerms) = assessLocalization(\@{$assignedLocalization{$gene}{'supported'}}, \%informationContent, \%children, \@{$geneLocalizations{$gene}});
      $correct{'supported'} += $match;
		}
		if ($assignedLocalization{$gene}{'approved'}) {
      $total{'approved'}++;
      my ($icWorst, $match, $noInitialTerms) = assessLocalization(\@{$assignedLocalization{$gene}{'approved'}}, \%informationContent, \%children, \@{$geneLocalizations{$gene}});
      $correct{'approved'} += $match;
		}
		if ($assignedLocalization{$gene}{'uncertain'}) {
      $total{'uncertain'}++;
      my @possibleTerms = @{$assignedLocalization{$gene}{'uncertain'}};
      my ($icWorst, $match, $noInitialTerms) = assessLocalization(\@possibleTerms, \%informationContent, \%children, \@{$geneLocalizations{$gene}});
      $correct{'uncertain'} += $match;
      if ($noInitialTerms > 0) {
        push @localizationsPerGene, $noInitialTerms;
      }

			my $tier = assessTier($icWorst, \@icCutoffs);
		  $totalGenes{$tier}++;
		  if ($match == 1){
			  $correctAssignments{$tier}++;
		  }
		}
	}
  my %emptyHash = ();
  outputSummary($outfile, \%correctAssignments, \%totalGenes, 0, \%emptyHash, \%emptyHash, 's');

  my $enhancedFrac = sprintf "%.3f", $correct{'enhanced'} / $total{'enhanced'};
  my $supportedFrac = sprintf "%.3f", $correct{'supported'} / $total{'supported'};
  my $approvedFrac = sprintf "%.3f", $correct{'approved'} / $total{'approved'};
  my $uncertainFrac = sprintf "%.3f", $correct{'uncertain'} / $total{'uncertain'};
	print $outfile "enhanced\t$enhancedFrac\t$total{'enhanced'}\t$correct{'enhanced'}\n";
	print $outfile "supported\t$supportedFrac\t$total{'supported'}\t$correct{'supported'}\n";
	print $outfile "approved\t$approvedFrac\t$total{'approved'}\t$correct{'approved'}\n";
	print $outfile "uncertain\t$uncertainFrac\t$total{'uncertain'}\t$correct{'uncertain'}\n";
} elsif ($fileType eq 'o') {
	while( my( $gene, $go ) = each %assignedLocalization ){
		$localizationAssessment{'total'}[$indexMap{'id'}{$go}]++;
		$totalGenes{'total'}++;
    my @possibleTerms = split /;/, $go;
    my ($icWorst, $match, $noInitialTerms) = assessLocalization(\@possibleTerms, \%informationContent, \%children, \@{$geneLocalizations{$gene}});
    $correctAssignments{'total'} += $match;
    $localizationAssessment{'known'}[$indexMap{'id'}{$go}] += $match;
    if ($noInitialTerms > 0) {
      push @localizationsPerGene, $noInitialTerms;
    }

		my $tier = assessTier($icWorst, \@icCutoffs);
    $totalGenes{$tier}++;
    if ($match == 1){
      $correctAssignments{$tier}++;
    }
	}
  outputSummary($outfile, \%correctAssignments, \%totalGenes, $noRanks, \%localizationAssessment, \%indexMap, 'o');
}

my $meanAssignedLocalizations = mean @localizationsPerGene;
my $medianAssignedLocalizations = median @localizationsPerGene;
print $outfile "\nMean localizations per gene\t$meanAssignedLocalizations\n";
print $outfile "Median localizations per gene\t$medianAssignedLocalizations\n";

close $outfile;
