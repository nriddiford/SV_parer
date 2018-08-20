#!/usr/bin/perl

package svParser;
use strict;
use warnings;
use autodie;

use feature qw/ say /;
use Data::Printer;

use Exporter;
our @ISA = 'Exporter';
our @EXPORT = qw( $VERSION );
our $VERSION = '1.1';


sub typer {
  my ($file, $type, $exclude_regions, $chrom_keys, $filters ) = @_;

  if ( $type eq 'l' or $type =~ m/lumpy/i ){
    say "Specified $file as a Lumpy file";
    $type = 'lumpy';
    parse($file, $type, $exclude_regions, $chrom_keys, $filters);
  }

  elsif ( $type eq 'd' or $type =~ m/delly/i ){
    say "Specified $file as a Delly file";
    $type = 'delly';
    parse($file, $type, $exclude_regions, $chrom_keys, $filters);
  }

  elsif ( $type eq 'n' or $type =~ m/novobreak/i ){
    say "Specified $file as a novoBreak file";
    $type = 'novobreak';
    parse($file, $type, $exclude_regions, $chrom_keys, $filters);
  }

  elsif ($type eq 'snp'){
    say "Forcing parsing of $file";
    $type = 'snp';
    parse($file, $type, $exclude_regions, $chrom_keys, $filters);
  }

  elsif ( $type eq 'guess' ){
    if ( `grep "source=LUMPY" $file` ){
      say "Recognised $file as a Lumpy file";
      $type = 'lumpy';
      parse($file, $type, $exclude_regions, $chrom_keys, $filters);
    }
    elsif ( `grep "DELLY" $file` ){
      say "Recognised $file as a Delly file";
      $type = 'delly';
      parse($file, $type, $exclude_regions, $chrom_keys, $filters);
    }
    elsif ( `grep "bamsurgeon spike-in" $file` ){
      say "Recognised $file as a novoBreak file";
      $type = 'novobreak';
      parse($file, $type, $exclude_regions, $chrom_keys, $filters);
    }
  }
  else {
    die "This VCF can not be parsed. Try specfiying type '-t' explicitly. See -h for details. Abort";
  }
}


sub parse {
  my ($file, $type, $exclude_regions, $chrom_keys, $filter_flags) = @_;
  open my $in, '<', $file or die $!;

  my @headers;
  my %filter_flags = %{ $filter_flags };

  my (%SVs, %info, %filtered_SVs, %call_lookup);
  my ($tumour_name, $control);
  my %format_long;
  my %info_long;
  my $filter_count;
  my @samples;
  my $replacement_id = 1;

  while(<$in>){
    chomp;
    if (/^#{2}/){
      push @headers, $_;
      $filtered_SVs{$.} = $_;

      if (/##FORMAT/){
        my ($format_long) = $_ =~ /\"(.*?)\"/;
        my ($available_format_info) = $_ =~ /ID=(.*?),/;
        $format_long{$available_format_info} = $format_long;
      }

      if (/##INFO/) {
        my ($info_long) = $_ =~ /\"(.*?)\"/;
        my ($available_info) = $_ =~ /ID=(.*?),/;
        $info_long{$available_info} = $info_long;
      }
      next;
    }

    if (/^#{1}/){
      push @headers, $_;
      $filtered_SVs{$.} = $_;
      my @split = split;
      push @samples, $_ foreach @split[9..$#split];

      $tumour_name = $samples[0];
      $control = $samples[1];
      next;
    }

    my @fields = split;

    my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, @sample_info) = @fields;

    # NovoBreak doesn't assign ids
    if ($id eq 'N' or $id eq '.'){
      $id = $replacement_id++;
    }

    my %sample_parts;

    push @{$sample_parts{$samples[$_]}}, split(/:/, $sample_info[$_]) for 0..$#samples;

    my @tumour_parts   = split(/:/, $sample_info[0]);
    my @normal_parts   = split(/:/, $sample_info[1]) if @samples > 1; # In case there are no control samples...

    my @format        = split(/:/, $format_block);
    my @info_parts    = split(/;/, $info_block);

    my %sample_info;

    for my $sample (@samples){
      for my $info (0..$#format){
        $sample_info{$id}{$sample}{$format[$info]} = $sample_parts{$sample}[$info];
      }
    }

    my @normals = @samples[1..$#samples];
    my @filter_reasons;

    my %information;

    foreach(@info_parts){
      my ($info_key, $info_value);

      if (/=/){
        ($info_key, $info_value) = $_ =~ /(.*)=(.*)/;
      }
      else {
        ($info_key) = $_ =~ /(.*)/;
        $info_value = "TRUE";
      }
      $information{$id}{$info_key} = $info_value;
    }

    my ($SV_type) = $info_block =~ /SVTYPE=(.*?);/;
    my ($SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $genotype, $filter_list);

    if ($type eq 'lumpy'){
      ( $SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $genotype, $filter_list ) = lumpy( $id, $chr, $info_block, $SV_type, $alt, $start, \%sample_info, $tumour_name, $control, \@samples, \@normals, \@filter_reasons, \%filter_flags );
    }
    elsif ($type eq 'delly'){
      next if $SV_type eq 'TRA'; # temp to resolve issues with svTyper...

      ( $SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $genotype, $filter_list ) = delly( $id, $info_block, $start, $SV_type, $tumour_name, $control, \@normals, \@filter_reasons, \%filter_flags, \%sample_info );
    }
    elsif ($type eq 'novobreak'){
      @samples = qw/tumour normal/;
      my ( $sample_info_novo, $format_novo, $format_long_novo );
      ( $SV_length, $chr2, $stop, $t_PE, $t_SR, $genotype, $filter_list, $sample_info_novo, $format_novo, $format_long_novo ) = novobreak( $id, $info_block, $start, $SV_type, $tumour_name, $control, \@sample_info, \@filter_reasons, \%filter_flags );
      %sample_info = %{ $sample_info_novo };
      $ab = "-";
      @format = @{ $format_novo };
      %format_long = %{ $format_long_novo };
    }
    elsif ($type eq 'snp'){
      $chr2 = $chr;
      $filter_list = \@filter_reasons;
    }

    if ( exists $filter_flags{'chr'} ){
      $filter_list = chrom_filter( $chr, $chr2, $filter_list, $chrom_keys );
    }

    $SV_length = abs($SV_length);

    ## NEW 25.7.18
    # Now filter any var with no SR support unless PE = su*2

    if ( $filter_flags{'su'} and $t_SR == 0 and ($t_PE <= $filter_flags{'su'} * 2) ){
      push @{$filter_list}, "$SV_type has no split read support and PE support < 2*$filter_flags{'su'}=" . $t_PE;
    }


    # Don't include DELS/INVs < 1kb with split read support == 0
    # Unless there is high PE support
    # if ( $filter_flags{'su'} ){
    #   if ( ($SV_type eq "DEL" or $SV_type eq "INV") and ( $SV_length < 1000 and $t_SR == 0 ) ){
    #       push @{$filter_list}, "$SV_type < 1kb with no split read support and PE support < 2*$filter_flags{'su'}=" . $t_PE if $t_PE <= $filter_flags{'su'} * 2;
    #     }
    #   elsif ( ($SV_type eq "BND" and $chr eq $chr2) and ( $SV_length < 1000 and $t_SR == 0 ) ){
    #       push @{$filter_list}, "$SV_type < 1kb with no split read support and PE support < 2*$filter_flags{'su'}=" . $t_PE if $t_PE <= $filter_flags{'su'} * 2;
    #     }
    # }

    # Filter for vars falling in an excluded region +/- slop
    if ( exists $filter_flags{'e'} and @$filter_list == 0 ){
      $filter_list = region_exclude_filter($chr, $start, $chr2, $stop, $exclude_regions, $filter_list);
    }

    $SVs{$id}  = [ @fields[0..10], $SV_type, $SV_length, $stop, $chr2, $t_SR, $t_PE, $ab, $filter_list, $genotype, \@samples ];
    $info{$id} = [ [@format], [%format_long], [%info_long], [@tumour_parts], [@normal_parts], [%information], [%sample_info] ];

    if (@$filter_list == 0){
      $filtered_SVs{$.} = $_;
      my $lookup_id = $chr . "_" . $start . "_" . $chr2 . "_" . $stop;
      $call_lookup{$lookup_id} = $id;
    }
  }
  return (\%SVs, \%info, \%filtered_SVs, \%call_lookup);
}


sub lumpy {
  my ( $id, $chr, $info_block, $SV_type, $alt, $start, $sample_info, $tumour, $control, $samples, $normals, $filters, $filter_flags ) = @_;
  my $genotype;
  my %filter_flags   = %{ $filter_flags };
  my @normals        = @{ $normals };
  my %sample_info    = %{ $sample_info };

  my ($SV_length) = $info_block =~ /SVLEN=(.*?);/;

  # This is all very redundant with the read support assignment. Combine. 21.12.17
  my ($t_hq_alt_reads, $n_hq_alt_reads, $n_sq) = (0,0,0);
  my %PON_alt_reads;

  # Non-high quality control read support for alt var
  my $c_alt_reads = $sample_info{$id}{$control}{'SU'};

  # will fail if not svtyper
  for my $normal (@normals[1..$#normals]){
    $PON_alt_reads{$normal} = $sample_info{$id}{$normal}{'QA'};
  }

  if (exists $sample_info{$id}{$tumour}{'QA'}){
    $t_hq_alt_reads = $sample_info{$id}{$tumour}{'QA'};
    $n_hq_alt_reads = $sample_info{$id}{$control}{'QA'};
    $n_sq = $sample_info{$id}{$control}{'SQ'};

    ($genotype, $filters) = genotype( $id, $t_hq_alt_reads, $c_alt_reads, $n_sq, $n_hq_alt_reads, \%PON_alt_reads, $filters );
  }
  # if genotype is set to 'NA' (this should be because there are no quality reads,
  # but some supporting reads) then re-genotype using supporting reads
  if ($genotype eq "NA" or not exists $sample_info{$id}{$tumour}{'QA'}){
    $t_hq_alt_reads = $sample_info{$id}{$tumour}{'SR'} + $sample_info{$id}{$tumour}{'PE'};
    $n_hq_alt_reads = $sample_info{$id}{$control}{'SR'} + $sample_info{$id}{$control}{'PE'};
    ($genotype, $filters) = genotype( $id, $t_hq_alt_reads, $c_alt_reads, "NA", $n_hq_alt_reads, \%PON_alt_reads, $filters );
  }

  my %PON_info; # What is this for? 25.7.18
  my ($c_HQ, $all_c_HQ) = (0,0);

  foreach my $normal (@normals){
    if ($sample_info{$id}{$normal}{'QA'}){
      $sample_info{$id}{$normal}{'QA'} eq '.' ? $sample_info{$id}{$normal}{'QA'} = '0' : $all_c_HQ += $sample_info{$id}{$normal}{'QA'};
      $PON_info{$normal} = [ $sample_info{$id}{$normal}{'QA'}, $sample_info{$id}{$normal}{'GT'} ];
    }
    else{
      $PON_info{$normal} = [ 0, $sample_info{$id}{$normal}{'GT'} ];
    }
  }

  if ( $filter_flags{'st'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'somatic_tumour' );
  }
  elsif ( $filter_flags{'sn'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'somatic_normal' );
  }
  elsif ( $filter_flags{'gp'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'germline_private' );
  }
  elsif ( $filter_flags{'gr'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'germline_recurrent' );
  }

  # my @samples = @{ $samples };
  my ($tumour_read_support, $direct_control_read_support) = ($t_hq_alt_reads, $n_hq_alt_reads);

  if ( $genotype ne 'somatic_tumour' ){
    # switch tumour normal
    my $tum2norm = $control;
    $control = $tumour;
    $tumour = $tum2norm;
    $tumour_read_support = $n_hq_alt_reads;
    $direct_control_read_support = $t_hq_alt_reads;
  }

#  my ($t_PE, $t_SR, $c_PE, $c_SR) = (0,0,0,0);
  my ($all_c_PE, $all_c_SR) = (0,0);
  my $t_PE = $sample_info{$id}{$tumour}{'PE'};
  my $t_SR = $sample_info{$id}{$tumour}{'SR'};

#  my ($tumour_read_support, $direct_control_read_support) = (0,0);

  # Get allele balance
  my $ab = '-';
  if($sample_info{$id}{$tumour}{'AB'}){
    $ab = $sample_info{$id}{$tumour}{'AB'};
  }

  ########################
  # Read support filters #
  ########################

  # Create temp pseudo counts to avoid illegal division by 0
  my $pc_tumour_read_support         = $tumour_read_support + 0.001;
  my $pc_direct_control_read_support = $direct_control_read_support + 0.001;

  if ( $filter_flags{'su'} ){
    $filters = read_support_filter($tumour_read_support, $filter_flags{'su'}, $tumour, $filters);
    if ($genotype eq 'germline_private' or $genotype eq 'germline_recurrent' ){
      # also require same read suppport for tumour
      # Maybe this is too harsh??
      $filters = read_support_filter($pc_direct_control_read_support, $filter_flags{'su'}, $control, $filters);
    }
  }

  my @filter_reasons = @{ $filters };

  ######################
  # Read depth filters #
  ######################

  if ( $filter_flags{'dp'} and $sample_info{$id}{$tumour}{'DP'} ){
    my $t_DP =  $sample_info{$id}{$tumour}{'DP'};

    if (scalar @{$samples} > 1){ # In case there are no control samples...
      my $c_DP =  $sample_info{$id}{$control}{'DP'};
      $c_DP = 0 if $c_DP eq '.';

      # Flag if either control or tumour has depth < 10 at site
      # slightly redundant section in fun (could just call fun once for each sample?)
      $filters = read_depth_filter($tumour, $control, $t_DP, $c_DP, $filter_flags{'dp'}, \@filter_reasons);
    }
    @filter_reasons = @{ $filters };

    # Subtract control reads from tumour reads
    # If this number of SU is less than 10% of tumour read_depth then filter
    if ( exists $filter_flags{'rdr'} and ( $tumour_read_support  / ( $t_DP + 0.01 ) ) < $filter_flags{'rdr'} ){
    # This is quite harsh (particularly for germline recuurent if theres poor cov in tum/normal...)
    # Modified 21.6.18 to print once if 'rdr' filter used
      push @filter_reasons, "$tumour\_reads/$tumour\_depth<" . ($filter_flags{'rdr'}*100) . "%" . '=' . $tumour_read_support . "/" . $t_DP;
    }
  }

  ##################
  # Quality filter #
  ##################

  if ( $sample_info{$id}{$tumour}{'SQ'} and exists $filter_flags{'sq'} ){
    $sample_info{$id}{$tumour}{'SQ'} = 0 if $sample_info{$id}{$tumour}{'SQ'} eq '.';

    if ( $sample_info{$id}{$tumour}{'SQ'} <= $filter_flags{'sq'} ){
      push @filter_reasons, "SQ<" . $filter_flags{'sq'} . '=' . $sample_info{$id}{$tumour}{'SQ'};
    }
  }

  my ($chr2, $stop) = ($chr, 0);
  if ($SV_type =~ /BND|TRA/){
    $chr2 = $alt =~ s/[\[\]N]//g;
    ($chr2, $stop) = $alt =~ /(.+)\:(\d+)/;
    $SV_length = $stop - $start;
  }
  else {
      ($stop) = $info_block =~ /;END=(.*?);/;
  }
  return ($SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $genotype, \@filter_reasons);
}


sub novobreak {
  my ( $id, $info_block, $start, $SV_type, $tumour_name, $control, $info, $filters, $filter_flags) = @_;

  my %filter_flags = %{ $filter_flags };

  my @info = @{ $info };

  my %sample_info;
  my %format_long;
  my $genotype;
  # for my $i (0 .. $#info){
  #   say "$i: $info[$i]";
  # }
  # Tumour reads
  my $bp1_SR = $info[9];  # high qual split reads bp1
  my $bp2_SR = $info[19]; # high qual split reads bp2
  my $bp1_PE = $info[26]; # PE reads bp1
  my $bp2_PE = $info[28]; # PE reads bp2

  my $t_SR = $bp1_SR + $bp2_SR;
  my $t_PE = $bp1_PE + $bp2_PE;

  # Normal reads
  my $n_bp1_SR = $info[14]; # high qual split reads bp1
  my $n_bp2_SR = $info[24]; # high qual split reads bp2
  my $n_bp1_PE = $info[27]; # PE reads bp1
  my $n_bp2_PE = $info[29]; # PE reads bp2


  # say "ID: $id";
  # say "bp1SR: $bp1_SR";
  # say "bp2SR: $bp2_SR";
  # say "bp1PE: $bp1_PE";
  # say "bp2PE: $bp2_PE";
  # say "---";

  my $n_SR = $n_bp1_SR + $n_bp2_SR;
  my $n_PE = $n_bp1_PE + $n_bp2_PE;

  # $t_SR = $t_SR/2;
  # $t_PE = $t_PE/2;

  # $t_SR = int($t_SR + 0.5);
  # $t_PE = int($t_PE + 0.5);

  my ($t_DP_1, $t_DP_2) = @info[6,16];

  if ($t_PE > ($t_DP_1 + $t_DP_2) ){
    $t_PE = 0; # Don't believe PE read support if greater than depth!!
    $n_PE = 0
  }

  my $tumour_read_support = ( $t_SR + $t_PE );
  my $all_control_read_support = ( $n_SR + $n_PE );

  ########################
  # Read support filters #
  ########################

  if ( $filter_flags{'su'} ){
    $filters = read_support_filter($tumour_read_support, $filter_flags{'su'}, 'tumour', $filters);
  }

  # We don't have a PON for novobreak - so create a dummy
  my %dummy_hash;
  $dummy_hash{'dummy'} = 0;

  ($genotype, $filters) = genotype( $id, $tumour_read_support, $all_control_read_support, 10, $all_control_read_support, \%dummy_hash, $filters );

  if ( $filter_flags{'st'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'somatic_tumour' );
  }
  elsif ( $filter_flags{'sn'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'somatic_normal' );
  }

  ######################
  # Read depth filters #
  ######################

  my ($c_DP_1, $c_DP_2) = @info[11,21];

  my $c_DP = ($c_DP_1 + $c_DP_2)/2;
  my $t_DP = ($t_DP_1 + $t_DP_2)/2;

  if ( $filter_flags{'dp'}){
    $filters = read_depth_filter('tumour', 'normal', $t_DP, $c_DP, $filter_flags{'dp'}, $filters);
  }

  my @tumour_parts = @info[6..10, 16..20, 26, 28];
  my @normal_parts = @info[11..15, 21..25, 27, 29];

  my @short_format = (
  "Breakpoint 1 depth",
  "Breakpoint 1 split reads",
  "Breakpoint 1 quality score",
  "Breakpoint 1 high quality split reads",
  "Breakpoint 1 high quality quality score",
  "Breakpoint 2 depth",
  "Breakpoint 2 split reads",
  "Breakpoint 2 quality score",
  "Breakpoint 2 high quality split reads",
  "Breakpoint 2 high quality quality score",
  "Breakpoint 1 discordant reads",
  "Breakpoint 2 discordant reads"
  );

  my @format = qw / DP1 SR1 Q1 HCSR1 HCQ1 DP2 SR2 Q2 HCSR2 HCQ2 PE1 PE2 /;

  $format_long{$format[$_]} = $short_format[$_] for 0..$#format;

  $sample_info{$id}{'tumour'}{$format[$_]} = $tumour_parts[$_] for 0..$#format;
  $sample_info{$id}{'normal'}{$format[$_]} = $normal_parts[$_] for 0..$#format;

  my ($stop) = $info_block =~ /;END=(.*?);/;
  my ($SV_length) = ($stop - $start);
  my ($chr2) = $info_block =~ /CHR2=(.*?);/;

  return ($SV_length, $chr2, $stop, $t_PE, $t_SR, $genotype, $filters, \%sample_info, \@format, \%format_long );
}


sub delly {
  my ($id, $info_block, $start, $SV_type, $tumour, $control, $normals, $filters, $filter_flags, $sample_ref) = @_;

  my %filter_flags   = %{ $filter_flags };
  my %sample_info    = % { $sample_ref };
  my @normals        = @{ $normals };

  # TUM alt PE
  my $dv = $sample_info{$id}{$tumour}{'DV'};
  # TUM alt SR
  my $rv = $sample_info{$id}{$tumour}{'RV'};
  # TUM ref PE
  my $dr = $sample_info{$id}{$tumour}{'DR'};
  # TUM ref SR
  my $rr = $sample_info{$id}{$tumour}{'RR'};

  # NORM alt PE
  my $n_dv = $sample_info{$id}{$control}{'DV'};
  # NORM alt SR
  my $n_rv = $sample_info{$id}{$control}{'RV'};
  # NORM ref PE
  my $n_dr = $sample_info{$id}{$control}{'DR'};
  # NORM ref SR
  my $n_rr = $sample_info{$id}{$control}{'RR'};

  if ($sample_info{$id}{$tumour}{'FT'} ne 'PASS'){
    push @{$filters}, "LowQual flag in tumour";
  }

  my $t_hq_alt_reads = $dv + $rv;
  my $tum_ref = $dr + $rr;
  my $n_hq_alt_reads = $n_dv + $n_rv;
  my $norm_ref = $n_dr + $n_rr;

  my %PON_alt_reads;
  my $genotype;
  my $n_sq = 0;

  $n_sq = 10;  ### Need to fix...
  for my $n (@normals[1..$#normals]){
    $sample_info{$id}{$n}{'DV'} eq '.' ? $sample_info{$id}{$n}{'DV'} = 0 : $sample_info{$id}{$n}{'DV'} = $sample_info{$id}{$n}{'DV'};
    $sample_info{$id}{$n}{'RV'} eq '.' ? $sample_info{$id}{$n}{'RV'} = 0 : $sample_info{$id}{$n}{'DV'} = $sample_info{$id}{$n}{'RV'};
    $PON_alt_reads{$n} = ( $sample_info{$id}{$n}{'DV'} + $sample_info{$id}{$n}{'RV'} );
  }
  ($genotype, $filters) = genotype( $id, $t_hq_alt_reads, $n_hq_alt_reads, $n_sq, $n_hq_alt_reads, \%PON_alt_reads, $filters );

  my $pc_tumour_read_support = $t_hq_alt_reads + 0.001;
  my $pc_direct_control_read_support = $n_hq_alt_reads + 0.001;

  my $ab = $pc_tumour_read_support/($pc_tumour_read_support+$tum_ref);

  my ($stop) = $info_block =~ /;END=(.*?);/;
  my ($chr2) = $info_block =~ /CHR2=(.*?);/;
  my ($SV_length) = ($stop - $start);
  # my ($t_SR, $t_PE) = ($rv, $dv);

  my ($t_SR, $t_PE) = (0,0);

  if ($info_block =~ /;SR=(\d+);/){
    $t_SR = $1;
  }

  if ($info_block =~ /;PE=(\d+);/){
    $t_PE = $1;
  }

  my $tumour_read_support = $t_hq_alt_reads;

  # my %PON_info;
  # #
  # foreach my $normal (@normals){
  #   $PON_info{$normal} = [ ($sample_info{$id}{$normal}{'DV'} + $sample_info{$id}{$normal}{'RV'}) , $sample_info{$id}{$normal}{'GT'} ];
  # }

  if ( $filter_flags{'st'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'somatic_tumour' );
  }
  elsif ( $filter_flags{'sn'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'somatic_normal' );
  }
  elsif ( $filter_flags{'gp'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'germline_private' );
  }
  elsif ( $filter_flags{'gr'} ){
    $filters = genotype_filter( $id, $genotype, $filters, 'germline_recurrent' );
  }

  ########################
  # Read support filters #
  ########################

  if ( $filter_flags{'su'} ){
    $filters = read_support_filter($tumour_read_support, $filter_flags{'su'}, $tumour, $filters);
    if ($genotype eq 'germline_private' or $genotype eq 'germline_recurrent' ){
      # also require same read suppport for tumour
      # Maybe this is too harsh??
      $filters = read_support_filter($pc_direct_control_read_support, $filter_flags{'su'}, $control, $filters);
    }
  }

  ######################
  # Read depth filters #
  ######################

  my $tum_depth = $pc_tumour_read_support + $tum_ref;
  my $norm_depth = $pc_direct_control_read_support + $norm_ref;

  if ( $filter_flags{'dp'}){
    $filters = read_depth_filter($tumour, $control, $tum_depth, $norm_depth, $filter_flags{'dp'}, $filters);
  }

  ##################
  # Quality filter #
  ##################

  if ( $sample_info{$id}{$tumour}{'SQ'} and exists $filter_flags{'sq'} ){
    $sample_info{$id}{$tumour}{'SQ'} = 0 if $sample_info{$id}{$tumour}{'SQ'} eq '.';

    if ( $sample_info{$id}{$tumour}{'SQ'} <= $filter_flags{'sq'} ){
      push @{$filters}, "SQ<" . $filter_flags{'sq'} . '=' . $sample_info{$id}{$tumour}{'SQ'};
    }
  }
  return ($SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $genotype, $filters );
}


sub summarise_variants {
  my ( $SVs, $filter_switch, $region ) = @_;
  my ($dels, $dups, $trans, $invs, $filtered) = (0,0,0,0,0);
  my ($tds, $CNVs, $ins) = (0,0,0);

  my ( $chromosome, $query_start, $query_stop, $specified_region) = parseChrom($region) if $region;
  my %support_by_chrom;
  my $read_support;
  my %filtered_sv;
  my %sv_type_count;
  my @gens = qw/germline_recurrent germline_private somatic_tumour somatic_normal/;
  my @types = qw/DEL DUP INV BND TRA DUP:TANDEM CNV INS /;

  for my $gen (@gens){
    for my $ty (@types){
      $sv_type_count{$gen}{$ty} = 0;
    }
  }

  for (keys %{ $SVs } ){
    my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $genotype, $samples ) = @{ $SVs->{$_} };

    if ( $chromosome ){
      next if $chr ne $chromosome;
    }
    if ( $specified_region ){
      if ( ( $start < $query_start or $start > $query_stop ) and ( $stop < $query_stop or $stop > $query_stop ) ){
         next;
       }
    }
     # Need this re-assignment for novobreak - should be harmless for lumpy and delly
    $id = $_;

    my @filter_reasons = @{ $filters };

    foreach (@filter_reasons){
      my ($reason) = $_ =~ /(.+)=/;
      $filtered_sv{$reason}++ ;
    }

    if ( scalar @filter_reasons > 0 ){
      $filtered++;
      next if $filter_switch;
    }

    $read_support = ( $SR + $PE );
    $support_by_chrom{$id} = [ $read_support, $sv_type, $chr, $SV_length, $start ];
    $sv_type_count{$genotype}{$sv_type}++;
  }
  print "\n";

  if ($filter_switch){
    say "Running in filter mode: $filtered calls filtered out:";
    say " - $_: $filtered_sv{$_}" for sort {$filtered_sv{$b} <=> $filtered_sv{$a} } keys %filtered_sv;
    print "\n";
  }

  printf "%-10s %-20s %-20s %-20s %-20s\n", "type", @gens;
  my $space_pad = "%-20s";
  for my $ty (@types){
    printf $space_pad, $ty;
    for my $gen (@gens){
      printf  $space_pad, "$sv_type_count{$gen}{$ty}" ;
    }
    print "\n";
  }

  my $top_count = 0;
  my %connected_bps;

  print "\nTop SVs by read count:\n";
  for ( sort { $support_by_chrom{$b}[0] <=> $support_by_chrom{$a}[0] } keys %support_by_chrom ){
    my $bp_id = $_;

    if ($bp_id =~ /_/){
      ($bp_id) = $bp_id =~ /(.+)?_/;
    }

    # flatten connected bps into 1 id for summary
    next if $connected_bps{$bp_id}++;
    $top_count++;

    print join("\n",
    "ID: $_",
    "TYPE: $support_by_chrom{$_}[1]",
    "CHROM: $support_by_chrom{$_}[2]",
    "START: $support_by_chrom{$_}[4]",
    "READS: $support_by_chrom{$_}[0]",
    "LENGTH: $support_by_chrom{$_}[3]\n") . "\n";

    last if $top_count >= 5;
  }
}


sub get_variant {
  my ($id_lookup, $SVs, $info, $filter_flag, $PON_print) = @_;
  if (not $info->{$id_lookup}){
    say "Couldn't find any variant with ID: '$id_lookup' in file. Abort";
    exit;
  }

  my (@format)     = @{ $info->{$id_lookup}->[0]};
  my (%format_long)   = @{ $info->{$id_lookup}->[1]};
  my (%info_long)    = @{ $info->{$id_lookup}->[2]};

  my (@tumour_parts)   = @{ $info->{$id_lookup}->[3]};
  my (@normal_parts)   = @{ $info->{$id_lookup}->[4]};

  my (%information)  = @{ $info->{$id_lookup}->[5]};
  my (%sample_info)  = @{ $info->{$id_lookup}->[6]};

  my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $genotype, $samples ) = @{ $SVs->{$id_lookup} };

  my @filter_reasons = @{ $filters };

  my @samples = @{ $samples };

  if ($PON_print > (scalar @samples - 2 )){
        $PON_print = scalar @samples - 2;
    }


  # Should change so that it will only print filter reasons if user specifies them
  if (scalar @filter_reasons > 0 ){
  say "\n________________________________________________________";
  say "Variant '$id_lookup' will be filtered for the following reasons:";
  say "* $_" foreach @filter_reasons;
  say "________________________________________________________\n";
  }

  printf "%-10s %-s\n",       "ID:",     $id_lookup;
  printf "%-10s %-s\n",       "TYPE:",   $sv_type;
  printf "%-10s %-s\n",       "GENOTYPE:",   $genotype;
  $chr2 ?
  printf "%-10s %-s\n",       "CHROM1:",   $chr :
  printf "%-10s %-s\n",       "CHROM:",  $chr;
  printf "%-10s %-s\n",       "CHROM2:",   $chr2 if $chr2;
  printf "%-10s %-s\n",       "START:",   $start;
  printf "%-10s %-s\n",       "STOP:",    $stop;
  ($chr2 and ($chr2 ne $chr) ) ? printf "%-10s %-s\n",   "IGV:",     "$chr:$start" : printf "%-10s %-s\n", "IGV:", "$chr:$start-$stop";
  printf "%-10s %-s\n",       "LENGTH:",   $SV_length unless $chr2 ne $chr2;
  printf "%-10s %-s\n",       "PE:",      $PE;
  printf "%-10s %-s\n",       "SR:",    $SR;
  printf "%-10s %-s\n",       "QUAL:",     $quality_score;
  printf "%-10s %-s\n",       "FILT:",     $filt;
  printf "%-10s %-s\n",       "REF:",     $ref;
  printf "%-10s %-s\n",       "ALT:",     $alt;

  my $separating_line = "___________________" x ($PON_print+4);
  say $separating_line;
  printf "%-20s",         "INFO";
  printf "%-20s",         $_ for @samples[0..$PON_print+1];
  printf "%-s\n",         "EXPLAINER";
  say $separating_line;

  foreach my $format_block (@format){
    printf "%-20s", "$format_block";

    foreach (@samples[0..$PON_print+1]){
      printf "%-20s", "$sample_info{$id_lookup}{$_}{$format_block}";
    }
    printf "%-s", "$format_long{$format_block}";
    print "\n";
  }

  say "______________________________________________________";
  printf "%-20s %-20s %-s\n", "INFO", "VALUE", "EXPLAINER";
  say "______________________________________________________";

  for (sort keys %{$information{$id_lookup}}){
    # turn off warnings for badly formatted novobreak vcf
    no warnings;
    printf "%-20s %-20s %-s\n", $_, $information{$id_lookup}{$_}, $info_long{$_};
  }
  say "______________________________________________________";

}


sub dump_variants {
  my ( $SVs, $info, $filter_flag, $region, $type, $PON_print ) = @_;

  my ( $chromosome, $query_start, $query_stop, $specified_region) = parseChrom($region) if $region;

  say "Running in filter mode - not displaying filtered calls" if $filter_flag;
  say "\nEnter any key to start cycling through calls or enter 'q' to exit";

  for ( sort { @{ $SVs->{$a}}[0] cmp @{ $SVs->{$b}}[0] or
        @{ $SVs->{$a}}[1] <=> @{ $SVs->{$b}}[1]
      }  keys %{ $SVs } ){

    my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $genotype, $samples ) = @{ $SVs->{$_} };

    $id = $_;

    if ( $chromosome ){
      next if $chr ne $chromosome;
    }

    if ( $specified_region ){
      if ( ( $start < $query_start or $start > $query_stop ) and ( $stop < $query_stop or $stop > $query_stop ) ){
        next;
      }
    }

    my (@format)        = @{ $info->{$_}->[0] };
    my (%format_long)   = @{ $info->{$_}->[1] };
    my (%info_long)     = @{ $info->{$_}->[2] };

    my (@tumour_parts)  = @{ $info->{$_}->[3] };
    my (@normal_parts)  = @{ $info->{$_}->[4] };

    my (%information)   = @{ $info->{$_}->[5] };
    my (%sample_info)   = @{ $info->{$_}->[6] };

    my @filter_reasons  = @{ $filters };

    my @samples         = @{ $samples };

    if ($PON_print > (scalar @samples - 2 )){
        $PON_print = scalar @samples - 2;
    }

    if ( scalar @filter_reasons > 0 ){
      next if $filter_flag;
    }

    my $next_line = <>;
    say "Displaying info for variant '$id'. Enter any key to go to the next variant or type 'q' to exit\n";

    if ( $next_line ){
      chomp($next_line);
      exit if $next_line eq 'q';

      if (scalar @filter_reasons > 0 ){
        say "_______________________________________________________";
        say "Variant '$id' will be filtered for the following reasons:";
        say "* $_" foreach @filter_reasons;
        say "________________________________________________________\n";
      }
      if ($type ne 'snp'){
        printf "%-10s %-s\n",        "ID:",         $id;
        printf "%-10s %-s\n",        "TYPE:",       $sv_type;
        printf "%-10s %-s\n",       "GENOTYPE:",   $genotype;
        $chr2 ?
        printf "%-10s %-s\n",       "CHROM1:",   $chr :
        printf "%-10s %-s\n",       "CHROM:",  $chr;
        printf "%-10s %-s\n",       "CHROM2:",      $chr2 if $chr2;
        printf "%-10s %-s\n",       "START:",       $start;
        printf "%-10s %-s\n",       "STOP:",        $stop;
        ($chr2 and ($chr2 ne $chr) ) ? printf "%-10s %-s\n",   "IGV:",     "$chr:$start" : printf "%-10s %-s\n", "IGV:", "$chr:$start-$stop";
        printf "%-10s %-s\n",       "LENGTH:",   $SV_length unless $chr2 ne $chr2;
        printf "%-10s %-s\n",       "PE:",          $PE;
        printf "%-10s %-s\n",       "SR:",          $SR;
        printf "%-10s %-s\n",       "QUAL:",        $quality_score;
        printf "%-10s %-s\n",       "FILT:",        $filt;
        printf "%-10s %-s\n",       "REF:",         $ref;
        printf "%-10s %-s\n",       "ALT:",         $alt;
      }

      elsif ($type eq 'snp'){
        printf "%-10s %-s\n",    "CHROM:",  $chr;
        printf "%-10s %-s\n",    "POS:",    $start;
        printf "%-10s %-s\n",    "IGV:",    "$chr:$start";
        printf "%-10s %-s\n",    "FILT:",   $filt;
        printf "%-10s %-s\n",    "REF:",    $ref;
        printf "%-10s %-s\n",    "ALT:",    $alt;
        printf "%-10s %-s\n",    "MUT:",    "$ref>$alt";
      }

        my $separating_line = "___________________" x ($PON_print+4);

        say $separating_line;

        printf "%-20s",         "INFO";
        printf "%-20s",         $_ for @samples[0..$PON_print+1];
        printf "%-s\n",         "EXPLAINER";

        say $separating_line;

      foreach my $format_block (@format){
        printf "%-20s",       $format_block;
        foreach (@samples[0..$PON_print+1]){
          printf "%-20s",     $sample_info{$id}{$_}{$format_block};
        }
        printf "%-s",         $format_long{$format_block};
        print "\n";
      }

      say "______________________________________________________";
      printf "%-20s %-20s %-s\n", "INFO", "VALUE", "EXPLAINER";
      say "______________________________________________________";

      for (sort keys %{$information{$id}}){
        # turn off warnings for badly formatted novobreak vcf
        no warnings;
        printf "%-20s %-20s %-s\n", $_, $information{$id}{$_}, $info_long{$_};
      }
      say "______________________________________________________";

    }
  }
}


sub print_variants {

  my ( $SVs, $filtered_SVs, $name, $output_dir, $germline ) = @_;

  open my $out, '>', $output_dir . $name . ".filtered.vcf" or die $!;
  say "Writing output to " . "'$output_dir" . $name . ".filtered.vcf'";

  my %filtered_SVs = %{ $filtered_SVs };
  my $sv_count = 0;

  for (sort {$a <=> $b} keys %filtered_SVs){
    my $line = $filtered_SVs{$_};

    if ($line =~ /^#/){
      print $out $line . "\n";
    }
    else {
      $sv_count++;
      my @cols = split("\t", $line);
      print $out join("\t", @cols[0..5], "PASS", @cols[7..$#cols]) . "\n";
    }
  }
  say "$sv_count variants passed all filters";
}


sub write_summary {
  my ( $SVs, $name, $summary_out, $type, $germline ) = @_;

  my $info_file;

  if ($germline){
    open $info_file, '>', $summary_out . $name . ".germline_filtered.summary.txt" or die $!;
    say "Writing useful info to " . "'$summary_out" . $name . ".germline_filtered.summary.txt'";
  }
  else{
    open $info_file, '>', $summary_out . $name . ".filtered.summary.txt" or die $!;
    say "Writing useful info to " . "'$summary_out" . $name . ".filtered.summary.txt'";
  }

  $type = "lumpy" if $type eq 'l';
  $type = "delly" if $type eq 'd';
  $type = "novobreak" if $type eq 'n';

  my %connected_bps;

  print $info_file join("\t", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2", "split_reads", "disc_reads", 'genotype', "id", "length(Kb)", "position", "consensus|type", "microhomology", "configuration", "allele_frequency", "mechanism|log2(cnv)") . "\n";

  for ( sort { @{ $SVs->{$a}}[0] cmp @{ $SVs->{$b}}[0] or
        @{ $SVs->{$a}}[1] <=> @{ $SVs->{$b}}[1]
      }  keys %{ $SVs } ){
    my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $genotype, $samples ) = @{ $SVs->{$_} };

    if ( @$filters == 0){

      my $bp_id = $_;
      if ($bp_id =~ /_/){
        ($bp_id) = $bp_id =~ /(.+)?_/;
      }

      # flatten connected bps into 1 id for summary
      next if $connected_bps{$bp_id}++;

      my ($length_in_kb) = sprintf("%.1f", abs($SV_length)/1000);

      $ab = sprintf("%.2f", $ab) unless $type eq 'novobreak' or $ab eq '.';

      my ($consensus, $mh_length, $ct, $rdr, $rde );

      # Consensus seq
      if ($info_block =~ /CONSENSUS=(.*?);/){
         $consensus = $1;
      }
      else{
        $consensus = "-";
      }

      # Read depth ratio (delly)
      if ($info_block =~ /RDRATIO=(\d+\.?\d*)/){
        $rdr = log($1)/log(2);
        $rdr = sprintf("%.2f", $rdr)
      }
      else{
        $rdr = '-';
      }

      # Microhology length (delly)
      if ($info_block =~ /HOMLEN=(\d+);/){
        $mh_length = $1;
      }
      else{
        $mh_length = "-";
      }

      # Configuration
      if ($info_block =~ /CT=(.*?);/){
        ($ct) = $1;
      }
      elsif ($alt =~ /\[|\]/) {
        $ct = $alt;
      }
      else {
        $ct = "-";
      }

      if ( $chr2 and ($chr2 ne $chr) ){
        print $info_file join("\t", $type, $sv_type, $chr, $start, $chr2, $stop, $SR, $PE, $genotype, $_, $length_in_kb, "$chr:$start $chr2:$stop", $consensus, $mh_length, $ct, $ab, $rdr ) . "\n";
      }
      else {
        print $info_file join("\t", $type, $sv_type, $chr, $start, $chr, $stop, $SR, $PE, $genotype, $_, $length_in_kb, "$chr:$start-$stop", $consensus, $mh_length, $ct, $ab, $rdr) . "\n";
      }

    }
  }
}


sub parseChrom {
  my $region = shift;

  if ( $region =~ /:/ ){
    my ($chromosome, $query_region) = split(/:/, $region);

    if ( $query_region !~ /-/ ){
      die "Error parsing the specified region.\nPlease specify chromosome regions using the folloing format:\tchrom:start-stop\n";
    }

    my ($query_start, $query_stop) = split(/-/, $query_region);
    $query_start =~ s/,//g;
    $query_stop =~ s/,//g;
    say "Limiting search to SVs within region '$chromosome:$query_start-$query_stop'";
    return($chromosome, $query_start, $query_stop, 1);
  }
  else {
    say "Limiting search to SVs on chromosome '$region'";
    return($region, undef, undef, 0);
  }
}
# Exclude variants where either breakpoint is withing +/- 250 bps of an unmappable region
# returns @filter_reasons
sub region_exclude_filter {
  my ( $chr1, $bp1, $chr2, $bp2, $exclude_regions, $filter_reasons ) = @_;

  my $slop = 250;

  my @filter_reasons = @{ $filter_reasons };

  my @bed = @{ $exclude_regions };

  foreach(@bed){
    my ($chromosome, $start, $stop) = split;
    next if $stop - $start < 200; # don't consider unmappable regions < 200 bps

    if ( $bp1 >= ($start - $slop) and $bp1 <= ($stop + $slop) ) {
      next unless $chromosome eq $chr1;
      push @filter_reasons, 'bp1 in or very close to excluded region=' . "$chromosome:$bp1 in:" . $start . '-' . $stop;
      last;
    }
    if ( $bp2 >= ($start - $slop) and $bp2 <= ($stop + $slop) ) {
      next unless $chromosome eq $chr2;
      push @filter_reasons, 'bp2 in or very close to excluded region=' . "$chromosome:$bp2 in:" . $start . '-' . $stop;
      last;
    }
  }
  return (\@filter_reasons);
}

# Remove variants where read support < $filter
# returns @filter_reasons
sub read_support_filter {
  my ($tumour_read_support, $read_support_flag, $sample, $filter_reasons ) = @_;
  my @filter_reasons = @{ $filter_reasons };

  # Filter if tum reads below specified threshold [default=4]
  if ( $tumour_read_support < $read_support_flag ){
    push @filter_reasons, "Read support in '$sample'<" . $read_support_flag . '=' . $tumour_read_support;
  }
  return(\@filter_reasons);
}


sub read_depth_filter {
  my ($tumour_name, $control, $tum_depth, $norm_depth, $depth_threshold, $filter_reasons) = @_;

  my @filter_reasons = @{ $filter_reasons };

  if ( $norm_depth < $depth_threshold ){
    push @filter_reasons, "$control has depth < " . $depth_threshold . '=' . $norm_depth;
  }

  if ( $tum_depth < $depth_threshold ){
    push @filter_reasons, "$tumour_name has depth < " . $depth_threshold . '=' . $tum_depth;
  }

  return(\@filter_reasons);
}


sub genotype_filter {
  my ($id, $genotype, $filter_reasons, $selected_genotype ) = @_;
  my @filter_reasons = @{ $filter_reasons };
    if ( $genotype ne $selected_genotype ){
      push @filter_reasons,  "Not $selected_genotype: genotype=" . $genotype;
    }
  return(\@filter_reasons);
}


sub genotype {
  my ( $id, $t_hq_alt_reads, $c_alt_reads, $n_sq, $n_hq_alt_reads, $PON, $filter_reasons ) = @_;
  my $genotype = 'NA';
  $n_sq = 20 if $n_sq eq "NA";
  $n_sq = 0 if  $n_sq eq ".";

  my $germline_recurrent = 0;
  my $tum = 0;
  my $norm = 0;
  my @filter_reasons = @{ $filter_reasons };

  if ( $t_hq_alt_reads > 0){
    $tum = 1;
  }
  if ( $n_hq_alt_reads > 0 and $n_sq > 1 ){
    $norm = 1;
  }
  # What was this for ?
  # if ( $c_alt_reads == 0){
  #   $norm = 0;
  # }

  my %PON_alt_reads = %{ $PON };

  for my $n ( keys %PON_alt_reads ){
    if ( $PON_alt_reads{$n} > 0){
      $germline_recurrent = 1;
    }
  }

  if ( $tum and not ($norm or $germline_recurrent) ){
    $genotype = 'somatic_tumour';
  }
  elsif ( $tum and $norm and not $germline_recurrent ){
    $genotype = 'germline_private';
  }
  elsif ( $tum or $norm and $germline_recurrent ){
    $genotype = 'germline_recurrent';
  }
  elsif ( $norm and not ($tum or $germline_recurrent) ){
    $genotype = 'somatic_normal';
  }
  elsif ( $germline_recurrent and not ($tum or $norm) ){
    $genotype = 'PON_var';
    push @filter_reasons, "PON member has quality support for alternative genotype=" . 1;
  }

  return($genotype, \@filter_reasons);
}


sub chrom_filter {
  my ( $chr, $chr2, $filters, $chrom_keys ) = @_;
  my @keys = @{ $chrom_keys };
  my %chrom_filt;

  $chrom_filt{$_} = 1 for (@keys);
  my @filter_reasons = @{ $filters };

  if ($chr2 eq '0'){
    $chr2 = $chr;
  }

  if ( not $chrom_filt{$chr} ){
    push @filter_reasons, 'chrom1=' . $chr;
  }

  elsif ( not $chrom_filt{$chr2} ){
    push @filter_reasons, 'chrom2=' . $chr2;
  }
  return (\@filter_reasons);
}


1;
