package Stats::Fishers;

use strict;
use warnings;
 
use Exporter qw(import);
use Text::NSP::Measures::2D::Fisher::right;
 
our @EXPORT_OK = qw(fishers);

# Fisher's exact test
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
 