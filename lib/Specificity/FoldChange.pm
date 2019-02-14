package Specificity::FoldChange;

use strict;
use warnings;
 
use Exporter qw(import);
 
our @EXPORT_OK = qw(specfc);

# Calculates specificity scores for all preys in an interation list.
sub specfc {
  my %saintInteractions = %{$_[0]};

  print STDERR "Calculating specificty scores\n";
  my $maxSpecificity = 0;
  my %specificity;

  my $numBaits = scalar keys %saintInteractions;
  foreach my $bait (keys %saintInteractions) {
    foreach my $prey (keys %{$saintInteractions{$bait}}) {
      # get average spectral count across all other baits
      my $avg = 0;
      foreach my $otherBait (keys %saintInteractions) {
        if (
          $otherBait ne $bait
          && exists $saintInteractions{$otherBait}{$prey}
        ) {
          $avg = $avg + $saintInteractions{$otherBait}{$prey}{'avgspec'};
        }
      }
      $avg = $avg / ($numBaits - 1);
      if ($avg == 0) {
        $specificity{$bait}{$prey} = "inf";
      } else {
        my $score = $saintInteractions{$bait}{$prey}{'avgspec'} / $avg;
        $specificity{$bait}{$prey} = $score;
        if ($score > $maxSpecificity) {
          $maxSpecificity = $score;
        }
      }
    }
  }

  # Replace all "inf" score with max specificity seen
  foreach my $bait (keys %specificity) {
    foreach my $prey (keys %{$specificity{$bait}}) {
      if ($specificity{$bait}{$prey} eq "inf") {
        $specificity{$bait}{$prey} = $maxSpecificity;
      }
    }
  }

  return \%specificity;
}
 