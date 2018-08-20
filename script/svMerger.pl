#!/usr/bin/perl
use strict;
use warnings;

use 5.18.2;

use File::Basename;
use File::Path qw(make_path);
use FindBin qw($Bin);
use feature qw/ say /;
use Data::Dumper;
use Getopt::Long qw/ GetOptions /;
use Cwd;

my @files;
my $help;

# foo=s{1,} indicates one or more values

GetOptions( 'files=s{1,}'   =>    \@files,
            'help'          =>    \$help
          ) or die usage();

if ($help) { exit usage() }

if (@files == 0){
  say "Exiting. Must list at least 1 file";
  exit usage();
}

my $dir = 'merged/';
# my $dir = "$Bin/../filtered/summary/merged/";
#
# eval { make_path($dir) };
# if ($@) {
#   print "Couldn't create $dir: $@";
# }

my @parts = split( /\./, basename($files[0]) );


open my $out, '>', "$dir" . $parts[0] . "_merged_SVs.txt" or die $!;

say "Writing merged files to: " . "'$dir" . $parts[0] . "_merged_SVs.txt'";


say "Merging files: ";
my %SVs;
my $header;
foreach (@files){
  say "* $_";
  open my $in, '<', $_ or die $!;

  while(<$in>){
    chomp;
    if ($. == 1){
      $header = $_;
      next;
    }
    my ($source, $chromosome1, $bp1, $chromosome2, $bp2) = (split)[ 0,2,3,4,5 ];

    push @{$SVs{$chromosome1}{$bp1}{$bp2}{$chromosome2}}, $_;

  }
}

print $out "$header\n";

for my $chr1 (sort keys %SVs){
  for my $bp1 (sort { $a <=> $b } keys %{ $SVs{$chr1} }){
    for my $chromosome2 (sort keys %{ $SVs{$chr1}{$bp1} }){
      for my $bp2 (sort keys %{ $SVs{$chr1}{$bp1}{$chromosome2} }){
        print $out "$_\n" foreach @{ $SVs{$chr1}{$bp1}{$chromosome2}{$bp2} };
      }
    }
  }
}

sub usage {
  print
"
usage: perl svMerger.pl [-h] [-f FILE]

svMerger
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v0.1
description: Merge summary output of SV calls from svParser

arguments:
  -h, --help            show this help message and exit
  -f FILE, --file
                        All summary files to be merged for a sample [required]
"
}
