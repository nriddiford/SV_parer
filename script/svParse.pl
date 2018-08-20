#!/usr/bin/perl

use strict;
use warnings;

use FindBin qw/ $Bin /;
use FindBin qw/ $Script /;
use File::Spec;
use lib File::Spec->catdir($FindBin::Bin, '..', 'bin/');

use svParser;

use feature qw/ say /;
use Data::Dumper;
use Data::Printer;
use File::Basename;
use File::Path qw/ make_path / ;
use File::Slurp;

use Getopt::Long qw/ GetOptions /;

my $vcf_file;
my $help;
my $id;
my $dump;
my $chromosome;
my $method = "guess";
my $print;
my $exclude = "";
my %filters;
my $PON_print = 5;
my %opt;
my $true_positives;

GetOptions( 'vcf=s'             =>    \$vcf_file,
            'method=s'          =>    \$method,
            'id=s'              =>    \$id,
            'dump'              =>    \$dump,
            'filter:s'          =>    \%filters,
            'print'             =>    \$print,
            'normalPrint=s'     =>    \$PON_print,
            'chromosome=s'      =>    \$chromosome,
            'exclude=s'         =>    \$exclude,
            'true_calls=s'      =>    \$true_positives,
            'help'              =>    \$help
) or die usage();

if ($help) { exit usage() }

if (not $vcf_file) {
   exit usage();
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

my ($filtered_out, $summary_out);

makeDirs() if $print;

my %legalFilters = ('su'  =>  1,
                    'dp'  =>  1,
                    'rdr' =>  1,
                    'sq'  =>  1,
                    'chr' =>  1,
                    'gr'  =>  1,
                    'gp'  =>  1,
                    'st'  =>  1,
                    'sn'  =>  1,
                    'e'   =>  1,
                    'a'   =>  1
);

my $filter_switch = 0;
my $filter_ref = \%filters;

if ( exists $filters{'a'} ){
  $filter_ref = allFilters(\%filters);
  %filters = %{ $filter_ref };
}

my $chroms;
$chroms = getChroms() if exists $filters{'chr'};

my $exclude_regions;
($filter_ref, $exclude_regions) =  readUnmappable($exclude, $filter_ref, $chroms) if $exclude;
($filter_switch, $filter_ref) = checkFilters($filter_ref, \%legalFilters, 0, $exclude, $chroms) if scalar keys %filters > 0;

print "\n";

# Retun SV and info hashes
my ( $SVs, $info, $filtered_vars, $call_lookup) = svParser::typer( $vcf_file, $method, $exclude_regions, $chroms, $filter_ref );
testCalls($true_positives, $SVs, $info, $filter_switch, $PON_print, $call_lookup) if $true_positives;

if ($method ne 'snp') {
  # Print summary to screen
  svParser::summarise_variants( $SVs, $filter_switch, $chromosome ) unless $id or $dump or $true_positives;
  # Print all info for specified id
  svParser::get_variant( $id, $SVs, $info, $filter_switch, $PON_print) if $id and not $true_positives;
  # Dump all variants to screen
  svParser::dump_variants( $SVs, $info, $filter_switch, $chromosome, $method, $PON_print ) if $dump and not $true_positives;
  # Write out variants passing filters
  # Write out some useful info to txt file
  svParser::print_variants( $SVs, $filtered_vars, $name, $filtered_out, 0 ) if $print;
  svParser::write_summary( $SVs, $name, $summary_out, $method, 0 ) if $print;
}


elsif ($method eq 'snp') {
  # Dump all variants to screen
  svParser::dump_variants( $SVs, $info, $filter_switch, $chromosome, $method) if $dump;
  if ($print or $id ){
    die "Print and get variants not supported for SNP data\n";
  }
}


sub makeDirs {
  $filtered_out = "$Bin/../filtered/";
  eval { make_path($filtered_out) };

  if ($@) {
    print "Couldn't create '$filtered_out': $@";
  }
  $summary_out = "$Bin/../filtered/summary/";
  eval { make_path($summary_out) };

  if ($@) {
    print "Couldn't create '$summary_out': $@";
  }
}

=item testCalls()

Test variants that pass filter against a
  list of true positives for that sample

=cut

sub testCalls {
  my ($truePositves, $SVs, $info, $filter_switch, $PON_print, $call_lookup) = @_;
  open my $tps, '<', $true_positives or die $!;
  my %truth;
  my %calls = %{ $call_lookup };
  my $found = 0;
  my $not_found = 0;
  my $total = 0;
  while(<$tps>){
    chomp;
    my ($sample, $caller, $type, $c1, $bp1, $c2, $bp2) = split;
    my $key = $c1 . "_" . $bp1 . "_" . $c2 . "_" . $bp2;
    $truth{$key} = 1;

    if (exists $calls{$key}){
      say "o True positve found in filtered data id: $calls{$key}";
      $found++;
    }
    else{
      say "!! True positve NOT found in filtered data!";
      say "  - $_";
      $not_found++;
    }
    $total++;
  }

  print "$found/$total true positves retained in filtered data\n";

  for my $k (keys %calls){
    if (not exists $truth{$k}){
      say "New variant found not in true positives: $calls{$k}";
    }
  }

}


sub allFilters {
  my $f = shift;
  my %filters = %{$f};
  $filters{'a'} = 1;
  $filters{'su'}  = 4;
  $filters{'dp'}  = 10;
  $filters{'rdr'} = 0.1;
  $filters{'sq'}  = 10;
  $filters{'chr'} = 1;
  return(\%filters)
}


sub checkFilters {
  my ($filter_given, $legalFilters, $filter_switch, $exclude, $chroms) = @_;
  my %legalFilters = %{$legalFilters};
  my %filters = %{$filter_given};
  say "Running in filter mode, using custom filters:";

  for my $k (keys %{ $filter_given } ){
    if (exists $legalFilters{$k}){
      explainFilters(\%legalFilters, $k, $exclude, $chroms, \%filters);
      $filter_switch = 1;
    }
    else {
      printIllegals($k);
    }
  }
  return($filter_switch, $filter_given);
}


sub explainFilters {
  my ($legals, $filter_given, $exclude, $chroms, $filters) = @_;
  my %legalFilters = %{ $legals };
  my %selectedFilters = %{ $filters };
      say " o Read support >= $selectedFilters{'su'}" if $filter_given eq 'su';
      say " o Read depth (in both tumour and normal) > $selectedFilters{'dp'}" if $filter_given eq 'dp';
      say " o Read support / depth > $selectedFilters{'rdr'}" if $filter_given eq 'rdr';
      say " o SQ quality > $selectedFilters{'sq'}" if $filter_given eq 'sq';
      say " o Chromosomes: " . join(' ', @{$chroms}) if $filter_given eq 'chr';
      say " o Running in germline PRIVATE mode" if $filter_given eq 'gp';
      say " o Running in germline RECURRANT mode" if $filter_given eq 'gr';
      say " o Running in somatic TUMOUR mode" if $filter_given eq 'st';
      say " o Running in somatic NORMAL mode" if $filter_given eq 'sn';
      say " o Excluding calls overlapping: $exclude" if $filter_given eq 'e';
}


sub printIllegals {
  my $illegal_option = shift;
  say "Illegal filter option used: '$illegal_option'. Please specify filters to run with (or use '-f or -f a' to run all defaults)";
  die usage();
}


sub getChroms {
  my @keys;
  if (-e 'chroms.txt'){
    say "\nReading chromosomes file from 'chroms.txt'";
    @keys = read_file('chroms.txt', chomp=>1);
  }
  else {
    say "\nFiltering for Drosophila chroms: 2L 2R 3L 3R 4 X Y ";
    @keys = qw / 2L 2R 3L 3R 4 X Y /;
  }
  return(\@keys);
}


sub readUnmappable {
  my ($x, $f, $c) = @_;
  my %filts = %{$f};
  my @keys = @{$c};
  my @exclude_regions;
  my @exclude = read_file($x, chomp=>1);
  if (exists $filts{'chr'}){
    my %chroms;
    @chroms{@keys} = ();
    # Should speed things up a bit by removing any chromosomes in exclude (.bed)
    # file that aren't in our list of chroms we want to look at (@keys)
    foreach(@exclude){
      my $chrom = (split)[0];
      next unless exists $chroms{$chrom};
      push @exclude_regions, $_;
    }
    undef @exclude;
  }
  $filts{'e'} = 1;
  return(\%filts, \@exclude_regions);
}


sub usage {
  print
"
usage: perl $Script -v variant_file.vcf [options]

svParser
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: $VERSION
description: Browse vcf output from several SV callers LUMPY, DELLY and novobreak

arguments:
  -h, --help                show this help message and exit
  -v file, --vcf            VCF input    [required]
  -p, --print               print filtered vcf and summary file to './filtered'
  -m str, --method          specify input source [default: guess from input]
                              -l = LUMPY
                              -d = DELLY
                              -n = novobreak
                              -snp = snp/snv

  -i str, --id              breakpoint id to inspect
  -d, --dump                print each variant called to screen, press any key to advance, q to exit

  -n, --normalPrint         maximum number of PON members to print to screen when using -d or -i [default: 5]
  -c str, --chromosome      limit search to chromosome and/or region (e.g. X:10000-20000)
                              can be used in conjunction with -d

  -e file, --exclude        path to .bed file containing regions to exlcude

  -f key=val, --filter
                            filter options:
                              -f su=INT    [number of tumour reads supporting var]
                              -f dp=INT    [minimum depth for both tumour normal at variant site]
                              -f rdr=FLOAT [fraction of supporting reads/tumour depth - a value of 1 would mean all reads support variant]
                              -f sq=INT    [phred-scaled variant likelihood]
                              -f chr=1     [only show chromosomes in 'chroms.txt'. [Default use Drosophila chroms: 2L 2R 3L 3R 4 X Y]
                              -f st=1      [only keep somatic tumour events]
                              -f sn=1      [only keep somatic normal events]
                              -f gp=1      [only keep germline private events]
                              -f gr=1      [only keep germline recurrent events]
                              -f, -f a     [apply default filters: -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -f chr=1 ]
"
}
