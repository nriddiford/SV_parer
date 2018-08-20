#!/usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use autodie;

use Data::Dumper;

my @files = @ARGV;

open my $all, '>', 'all_samples.txt';

open my $file1, '<', $ARGV[0];
my $header = <$file1>;
close $file1;

# my $header = join("\t", qw/sample event source type chromosome1 bp1 chromosome2 bp2 split_reads pe_reads id length(Kb) position consensus|type microhomology configuration allele_frequency mechanism|log2(cnv) bp1_locus bp2_locus affected_genes T|F notes/);
print $all "sample" . "\t" . $header . "\n";

my %seen;
my @lines;
for (@files) {
  open my $in, '<', $_;
  my ($sample) = split(/\_/);
    while(<$in>){
      chomp;
      next if /^event/;
      push @lines, "$sample\t$_" unless $seen{$_}++;
    }
}

print $all "$_\n" foreach @lines;
