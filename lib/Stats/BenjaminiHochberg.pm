package Stats::BenjaminiHochberg;

use strict;
use warnings;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(bhCorrection);

# Benjamini-Hochberg correction
sub bhCorrection {
	my @sortedPValues;
	my @sortedTerms;
	my %unsortedHash = %{$_[0]};
	my $fdr = $_[1];
	# order p-values
	foreach my $term (sort { $unsortedHash{$a} <=> $unsortedHash{$b} } keys %unsortedHash) {
		push @sortedPValues, $unsortedHash{$term};
		push @sortedTerms, $term;
	}
	my $noTests = scalar @sortedPValues;
	# set rank for each
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
 