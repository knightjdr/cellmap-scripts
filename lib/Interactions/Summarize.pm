package Interactions::Summarize;

use strict;
use warnings;
 
use Exporter qw(import);
use Text::CSV_XS;
 
our @EXPORT_OK = qw(summarize summarizeTop);

# Takes detected SAINT interactions and summarizes how many were previoulsy known.
sub summarize {
  my %saintInteractions = %{$_[0]};
  my %knownInteractions = %{$_[1]};
  my %baitMap = %{$_[2]};

  print STDERR "Getting list of interactors\n";
  my %recovered;
  my %stats;
  foreach my $bait (keys %saintInteractions) {
    $stats{$bait}{'discovered'} = 0;
    $stats{$bait}{'recovered'} = 0;
    if (exists $knownInteractions{lc $baitMap{$bait}}) {
      $stats{$bait}{'missed'} = scalar @{$knownInteractions{lc $baitMap{$bait}}};
    } else {
      $stats{$bait}{'missed'} = 0;
    }
    my %knownPreys = map { $_ => 1 } @{$knownInteractions{lc $baitMap{$bait}}};
    foreach my $prey (keys %{$saintInteractions{$bait}}) {
      if (exists $knownPreys{lc $prey}) {
        $stats{$bait}{'missed'}--;
        $stats{$bait}{'recovered'}++;
        push @{$recovered{$bait}}, $prey;
      } else {
        $stats{$bait}{'discovered'}++;
      }
    }
  }

  return (\%recovered, \%stats);
}

# Takes detected SAINT interactions and summarizes how many were previoulsy known.
sub summarizeTop {
  my %saintInteractions = %{$_[0]};
  my %knownInteractions = %{$_[1]};
  my %baitMap = %{$_[2]};

  my %recovered;
  my %stats;
  foreach my $bait (keys %saintInteractions) {
    $stats{$bait}{'discovered'} = 0;
    $stats{$bait}{'recovered'} = 0;
    my %knownPreys = map { $_ => 1 } @{$knownInteractions{lc $baitMap{$bait}}};
    foreach my $prey (@{$saintInteractions{$bait}{'list'}}) {
      if (exists $knownPreys{lc $prey}) {
        $stats{$bait}{'recovered'}++;
        push @{$recovered{$bait}}, $prey;
      } else {
        $stats{$bait}{'discovered'}++;
      }
    }
  }

  return (\%recovered, \%stats);
}
 