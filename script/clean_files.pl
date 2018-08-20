#!/usr/bin/perl
use strict;
use warnings;
use feature qw/ say /;
use autodie;
use File::Basename;

use Data::Printer;


sub read_files {
  my $fp_file = 'all_samples_false_calls.txt';
  my $tp_file = 'all_samples_whitelist.txt';

  my $tp = get_TP($tp_file);
  my $fp = get_FP($fp_file);

  open my $false_calls_file, '>>', 'all_samples_false_calls.txt';
  open my $true_calls_file, '>>', 'all_samples_whitelist.txt';

  open my $sv_file, '<', $ARGV[0];

  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($ARGV[0]), 2);
  my ($sample) = split(/_/, $name, 3);
  my $out_file = $sample . "_" . "cleaned_SVs.txt";

  # This allows files to be parsed irrespective of line ending. Not recommened for large files
  my @lines;
  {
  local $/ = undef;
  my $content = <$sv_file>;
  @lines = split /\r\n|\n|\r/, $content;
  }

  return($out_file, $false_calls_file, $true_calls_file, $sample, \@lines, $tp, $fp);
}

sub mark_events{
  my ($out_file, $false_calls_file, $true_calls_file, $sample, $l, $tp, $fp) = read_files();
  my @lines       = @{ $l };
  my %false_calls = %{ $fp };
  my %true_calls  = %{ $tp };
  my %events;

  open my $cleaned_file, '>', $out_file;
  my $filtered_calls = 0;
  my $whitelisted_calls = 0;
  my $new_fp = 0;
  my $new_tp = 0;

  for (@lines){
    chomp;
    my @cells = split(/\t/);

    if ($cells[0] eq 'event'){
      print $cleaned_file $_ . "\n";
      next;
    }

    my $match = join("_", $sample, @cells[3..6]);

    my $notes = "-";
    if (@cells >= 21){
      $notes = $cells[21];
      # p(@cells);
    }

    if ($notes eq 'F') {
      $filtered_calls++;
      if (not $false_calls{$match}){
        say "Adding new false call to blacklist: $match";
        $new_fp++;
        print $false_calls_file $match . "\n";
        # print fp_out
      }
    }
    elsif ($false_calls{$match}) {
      say "Unmarked false call: $match";
      $filtered_calls++;
    }

    elsif ($notes eq 'T') {
      if ( not $true_calls{$match} ){
        say "Added new call to whitelist: $match";
        print $true_calls_file $match . "\n";
        $new_tp++;
        $whitelisted_calls++;
      }
      else{
        $whitelisted_calls++;
      }
    }

    if ($notes ne 'F'){
      print $cleaned_file $_ . "\n";
    }
  }

  if ($filtered_calls>0){
    printf "Removed %s calls marked as false positives from %s [%s new]\n", $filtered_calls, $sample, $new_fp;
  }
  if ($whitelisted_calls>0){
    printf "%s calls in whitelist in %s [%s new]\n", $whitelisted_calls, $sample, $new_tp;
  }
}
sub get_TP {
  my $tp_file = shift;
  open my $true_positives, '<', $tp_file if -e $tp_file;
  my %true_calls;
  while(<$true_positives>){
    chomp;
    my $id = (split)[0];
    $true_calls{$id}++;
  }
  return(\%true_calls);
}

sub get_FP {
  my $fp_file = shift;
  open my $false_positives, '<', $fp_file if -e $fp_file;
  my %false_calls;
  while(<$false_positives>){
    chomp;
    $false_calls{$_}++;
  }
  return(\%false_calls);
}

mark_events();
