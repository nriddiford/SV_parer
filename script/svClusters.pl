#!/usr/bin/env perl

use strict;
use warnings;

use autodie;
use Data::Dumper;
use File::Basename;
use feature qw/ say /;

use FindBin qw($Bin);
use Cwd;


my $in = $ARGV[0];

# my $dir = "$Bin/../filtered/summary/merged/";
my $dir = getcwd;

my @parts = split( /\_/, basename($in) );

open my $out, '>', "$dir/" . $parts[0] . "_clustered_SVs.txt" or die $!;
say "Writing clustered sv file to: " . "'$dir/" . $parts[0] . "_clustered_SVs.txt'";
#
# open my $out, '>', $parts[0] . "_clustered_SVs.txt" or die $!;
# say "Writing clustered sv file to: " . $parts[0] . "_clustered_SVs.txt'";

my $sample_count;

run($in);

sub run {
    my $file = shift;

    open my $fh, '<', $file;

    my @header = split(/\t/, scalar <$fh>);
    chomp(@header);

    my @events = ([]);

    while (<$fh>) {
      chomp;
      $sample_count++;
      my %event;
      @event{ @header } = split(/\t/);

        # $events[-1] last element of array
        unless (@{ $events[-1] }) {
            push @{ $events[-1] }, \%event;
            next;
        }

        if (same_event($events[-1][-1], \%event, 50)) {
            push @{ $events[-1] }, \%event;
            next;
        }

        push @events, [ \%event ];
    }

    print $out join("\t", event => @header), "\n";

    say $sample_count . " calls were grouped into " . scalar @events . " events for " . $parts[0];

    for my $i (1 .. @events) {
        for my $ev (@{ $events[$i - 1] }) {
            print $out join("\t", $i, @{$ev}{@header}), "\n";
        }

    }
}

sub same_event {
    my ($x, $y, $threshold) = @_;
    return if $x->{chromosome1} ne $y->{chromosome1};
    return if abs($x->{bp1} - $y->{bp1}) > $threshold;
    return if abs($x->{bp2} - $y->{bp2}) > $threshold;
    return 1;
}
