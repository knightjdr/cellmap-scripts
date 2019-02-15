package Prey::InteractorOutput;

use strict;
use warnings;
 
use Exporter qw(import);
use List::MoreUtils qw(uniq);
 
our @EXPORT_OK = qw(output);

# Outputs summary statistics for prey-prey correlation
sub output {
  my $filename = $_[0];
  my @correlation = @{$_[1]};
  my $total = $_[2];
  my $limit = $_[3];

  print STDERR "Formatting output\n";
  open my $fh, '>', $filename;
  if ($total > 0) {
    print $fh "cutoff\tdiscovered\trecovered\tfraction\n";
  } else {
    print $fh "cutoff\tdiscovered\trecovered\tmissed\tfraction\n";
  }
  for(my $i = 0; $i <= $limit; $i++) {
    my @discoveredArray;
    my @recoveredArray;
    for(my $j = $i; $j <= $limit; $j++) {
      if ($correlation[$j]) {
        if (exists $correlation[$j]{'discovered'}) {
          push @discoveredArray, @{$correlation[$j]{'discovered'}};
        }
        if (exists $correlation[$j]{'recovered'}) {
          push @recoveredArray, @{$correlation[$j]{'recovered'}};
        }
      }
    }
    my $discovered = scalar uniq @discoveredArray;
    my $recovered = scalar uniq @recoveredArray;
    my $fraction = sprintf "%.3f", $recovered / ($discovered + $recovered);
    my $correlationCutoff = sprintf "%.2f", $i / 100;
    if ($total > 0) {
      my $missed = $total - $recovered;
      print $fh "$correlationCutoff\t$discovered\t$recovered\t$missed\t$fraction\n"
    } else {
      print $fh "$correlationCutoff\t$discovered\t$recovered\t$fraction\n"
    }
  }
  close $fh;
}
