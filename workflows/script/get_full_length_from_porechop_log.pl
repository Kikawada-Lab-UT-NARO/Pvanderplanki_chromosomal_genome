#!/usr/bin/env perl
use strict;
use warnings;
use G;

my $count = 0;
my $total = 0;
my ($read);
my ($start_flag, $end_flag) = (0,0);
my $adapter;
my (%start, %end);

my %adapter_which = (
   "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG" => "read_1",
  "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG" => "read_2"
);

my $input_file = shift;
open my $fh, ">", $input_file.".plus";
open my $fg, ">", $input_file.".neg";

foreach my $line ( readFile( $input_file ,1 , -format=>"plain" )) {
  print STDERR "\r".$count;
  if($line =~ /runid/ ) {
    $total++;
    $read = (split /\t/, $line)[0];
    $start_flag = 0;
    $end_flag = 0;
    %start = ();
    %end = ();
  }
  elsif($line =~ /    start alignments/) {
    $start_flag = 1;
  }
  elsif( $line =~ /      Nextera Transposase Adapters/ && $start_flag == 1 && $end_flag == 0 ) {
    $line =~ /Nextera Transposase Adapters \(read\d-leading, read\d-comp\), full score=(\d+.\d+), partial score=(\d+.\d+), read position: (\d+)-(\d+) ([ATGC]+)/;
    my ($full, $partial, $read_pos_5, $read_pos_3, $seq) = ($1, $2, $3, $4, $5);
    my $tmp = defined $adapter_which{$seq} ? $adapter_which{$seq} : "NOHIT";
    $start{$tmp} = $full;

  }
  elsif( $line =~ /    end alignments/ ) {
    $end_flag = 1;
  }
  elsif( $line =~ /      Nextera Transposase Adapters/ && $start_flag == 1 && $end_flag == 1 ) {
    $line =~ /Nextera Transposase Adapters \(read\d-leading, read\d-comp\), full score=(\d+.\d+), partial score=(\d+.\d+), read position: (\d+)-(\d+) ([ATGC]+)/;
    my ($full, $partial, $read_pos_5, $read_pos_3, $seq) = ($1, $2, $3, $4, $5);
    my $tmp = defined $adapter_which{complement($seq)} ? $adapter_which{complement($seq)} : "NOHIT";
    $end{$tmp} = $full;
  }
  elsif( $line =~ /^$/ ) {
    if( $start_flag == 1 && $end_flag == 1 ) {
#      say $read;

      my( $start_adapter, $end_adapter);
      if( ~~(keys %start) == 2 ) {
	my $read1 = $start{read_1};
	my $read2 = $start{read_2};
	$start_adapter = $start{read_1} > $start{read_2} ? "read_1" : "read_2";
      }
      else {
	$start_adapter = defined $start{read_1} ? "read_1" : "read_2";
      }
      if( ~~(keys %end) == 2 ) {
	my $read1 = $end{NOHIT};
	my $read2 = $end{read_2};
	$end_adapter = $end{NOHIT} > $end{read_2} ? "NOHIT" : "read_2";
      }
      else {
	$end_adapter = defined $end{NOHIT} ? "NOHIT" : "read_2";
      }

      my $flag = 0;
      if( $start_adapter eq "read_1" && $end_adapter eq "read_2" ) {
	$flag = 1;
	print $fh $read."\n";
      }
      elsif( $start_adapter eq "read_2" && $end_adapter eq "NOHIT" ) {
	$flag = 1;
	print $fg $read."\n";

      }
      if( $flag == 1 ) {
#	say $start_adapter."\t".$end_adapter;
	$count++;
      }
    }
  }
}
print STDERR "\r".$count."/".$total."\n";

# reads had adapters trimmed from their start

close $fh;
close $fg;
