#!/usr/bin/perl
use strict;
use warnings;
use feature qw/ say /;
use Data::Printer;
use Data::Dumper;
use autodie;

use File::Basename;
use FindBin qw($Bin);



my $cnvs = shift;
my $cnv_ref = extractVars($cnvs);


my @name_fields = split( /\_/, basename($cnvs) );
my $dir = "$Bin/../filtered/summary/";

open my $out, '>', "$dir" . $name_fields[0] . ".cnvseq.filtered.summary.txt";
print $out join("\t", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2", "split_reads", "disc_reads", 'genotype', "id", "length(Kb)", "position", "consensus|type", "microhomology", "configuration", "allele_frequency", "mechanism|log2(cnv)") . "\n";

my @lines = @{$cnv_ref};
print $out "$_\n" foreach @lines;


sub extractVars {
  my $in = shift;
  open my $annotated_cnvs, '<', $in;
  my @cnv;
  while(<$annotated_cnvs>){
    chomp;
    my ($chrom, $start, $stop, $type, $length, $fc) = split;

    push @cnv, join("\t", "CNV-Seq",               # source
                            $type,                 # type
                            $chrom,                # chrom1
                            $start,                # bp1
                            $chrom,                # chrom2
                            $stop,                 # bp2
                            '-',                   # split reads
                            '-',                   # paired reads
                            'somatic_tumour',      # genotype
                            '-',                   # id
                            $length,               # length
                            "$chrom:$start-$stop", # IGV
                            '-',                   # misc1 (type)
                            '-',                   # microhomology
                            '-',                   # configuration
                            '-',                   # allele frequency
                            $fc);                  # misc (cnv)

  }
  p(@cnv);
  return(\@cnv);
}
