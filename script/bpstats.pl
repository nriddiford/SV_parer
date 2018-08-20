#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use feature qw/ say /;
use Data::Dumper;

use Getopt::Long qw/ GetOptions /;

my $notch = 0;

GetOptions( 'notch'  => \$notch );

my (%samples, %chroms, %genes, %gene_data, %gene_sample, %features, %types, %seen_event, %seen_event_bp);

# Read in 'all_bps_cleaned.txt'
my $in_file = $ARGV[0];
open my $in, '<', $in_file;

my @omit = qw/ A373R1 A512R17 A373R7 /;

print "\n";

say "Omitting sample $_" for @omit;

print "\n";

say "Exclduding bps around Notch" if $notch;

my $sv_count = 0;

while(<$in>){
  chomp;
  my ($event, $bp_id, $sample, $chrom, $bp, $gene, $feature, $type, $length) = (split);
  next if grep /$sample/, @omit;

  if ($notch){
    if ($chrom eq 'X' and $bp >= 3000000 and $bp <= 3300000){
      next;
    }
  }
  $feature =~ s/_\d+//g;

  unless ($seen_event{$sample}{$event}++){
    $sv_count++;
    $samples{$sample}++;
    $chroms{$chrom}++;
    $types{$type}++;
  }

  unless ( $seen_event_bp{$sample}{$event}{$bp_id}++ ){
    $genes{$gene}++ unless $gene eq 'intergenic';
    push @{$gene_sample{$gene}} , $sample;
    $gene_data{$gene}{$sample}++ unless $gene eq 'intergenic';
    $features{$feature}++;
    if ($notch){
    print join ("\t", $_ ) . "\n";
  }
  }

}

my $sample_count = keys %samples;
print "$sv_count structural variants detected across " . $sample_count . " samples\n";

print "\n";

my $top_samples = 0;
print "Top 10 samples for structural variant count:\n";
for ( sort { $samples{$b} <=> $samples{$a} } keys %samples ){
  say "$_: " . $samples{$_};
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

print "\n";

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
