#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use feature qw/ say /;

use Data::Dumper;

my (%samples, %chroms, %genes, %gene_data, %gene_sample, %features, %types);

my $in_file = $ARGV[0];
open my $in, '<', $in_file;

my @omit = qw/ A373R1 A512R17 A373R7 /;

print "\n";

say "Omitting sample $_" for @omit;

print "\n";

while(<$in>){
  chomp;
  my ($sample, $chrom, $bp, $gene, $feature, $type) = (split);
  next if grep /$sample/, @omit;

  $feature =~ s/_\d+//g;

  $types{$type}++;
  $samples{$sample}++;
  $chroms{$chrom}++;
  $gene_data{$gene}{$sample}++ unless $gene eq 'intergenic';
  $genes{$gene}++ unless $gene eq 'intergenic';
  
  push @{$gene_sample{$gene}} , $sample;
  $features{$feature}++ unless $feature eq 'intergenic';
}

my $top_samples = 0;
print "Top 10 samples for structural variant count:\n";
for ( sort { $samples{$b} <=> $samples{$a} } keys %samples ){
  say "$_: " . ($samples{$_}/2);
  $top_samples++;
  last if $top_samples == 10;
}

print "\n";

print "Structural variant calls per chromosome for all samples:\n";
say "$_: " . $chroms{$_} for sort { $chroms{$b} <=> $chroms{$a} } keys %chroms;
print "\n";


my %genes_by_sample;
for my $gene ( keys %gene_data ){
  my $sample_count = 0;
  my $found;
  for my $sample ( sort { $gene_data{$gene}{$b} <=> $gene_data{$gene}{$a} } keys %{$gene_data{$gene}} ){
    $sample_count++;
    $genes_by_sample{$gene} = $sample_count;
  }
}

my $top_genes_by_sample = 0;
print "Genes hit in multiple samples\n";

for my $hit_genes (sort { $genes_by_sample{$b} <=> $genes_by_sample{$a} } keys %genes_by_sample ) {
  say "$hit_genes: " . "$genes_by_sample{$hit_genes} " . "[" . join(", ", keys %{$gene_data{$hit_genes}}) . "]";
  $top_genes_by_sample++;
  last if $top_genes_by_sample == 10;
}

print "\n";

my $top_genes = 0;
print "Most hit genes\n";
for ( sort { $genes{$b} <=> $genes{$a} } keys %genes ){
  say "$_: " . "$genes{$_} " . "[" . join(", ", @{$gene_sample{$_}}) . "]";
  $top_genes++;
  last if $top_genes == 10;
}

print "Types of event:\n";
for ( sort { $types{$b} <=> $types{$a} } keys %types ){
  say "$_: " . $types{$_};
}

print "\n";

my $top_features = 0;
print "10 most hit features accross all samples\n";
for ( sort { $features{$b} <=> $features{$a} } keys %features ){
  print join("\t", $_, $features{$_}) . "\n";
  $top_features++;
  last if $top_features == 10;
}

print "\n";