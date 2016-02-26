#!/usr/bin/perl
#
#  SubSample_Fastq.pl
#  
#  Copyright 2013 Dylan <dylan.storey@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
####

use warnings;
use strict;
use Getopt::Long;
use List::Util 'shuffle';
use File::Basename;

my $total_reads = 0;
my $desired_reads = 0;
my $paired_end = 1;
my $pair_1 = '';
my $pair_2 = '';
my $single_end;
my $interleave_out = 0;
my $output_dir = '.'; #mkuhn, 20140305: output dir
my $name_suffix = ''; #mkuhn, 20140306: extra string
my $verbose = 0; #mkuhn, 20140305: verbose output on screen
my $help = 0;

GetOptions(	"total=i"   => 	\$total_reads,
			"desired=i" =>	\$desired_reads,
			"paired"  =>	\$paired_end,
			"1=s" 	  => 	\$pair_1,
			"2=s"	  	  =>	\$pair_2,
			"interleave_out=i" => \$interleave_out,
      "output=s"  => \$output_dir,
      "name=s"    => \$name_suffix,
      "verbose"   => \$verbose,
			"help|?" => \$help,
);

print_usage() if ( $help || ! length $pair_1 );


# mkuhn, 20140305
die "Output directory does not exist." if (! -d $output_dir);

if ($paired_end == 1) {
	die "file $pair_1 doesn't exist" if (! -e $pair_1);
	die "file $pair_2 doesn't exist" if (! -e $pair_2);
	if ($total_reads == 0){
		my $total_reads_1 = $1 if (`wc -l $pair_1` =~ /(\d+)/);
		my $total_reads_2 = $1 if (`wc -l $pair_2` =~ /(\d+)/);

		#File Sanity Checks ( Bare minimum )
		die "files not of same length" if ($total_reads_1 ne $total_reads_2);
		die "Malformed Fastq $pair_1" if ($total_reads_1 % 4 ne 0);
		die "Malformed Fastq $pair_2" if ($total_reads_2 % 4 ne 0);
		
    #mkuhn, total_reads count read pairs
		#$total_reads = ($total_reads_1 + $total_reads_2) / 4 ; #total records = total lines divide 4
    $total_reads = $total_reads_1 / 4 ; # number of paired end reads
		}
	
  # mkuhn, 20140305: do not have more reads desired than there are
  die "Too many desired reads!" if ($desired_reads > $total_reads);

  if ( $verbose ) {
    ## PRINT OPTIONS TO SCREEN #	
  	print "Current Options Are:
  	Total Reads: $total_reads
  	Desired Read: $desired_reads
  	";
    print"Output Directory:\t$output_dir\n\t";
  	print"Paired end reads:\t$pair_1 \t $pair_2 \n" if $paired_end == 1;
  	print"Single end reads:\t$single_end" if $paired_end == 0;
  }

## Get our Randomly selected record numbers 
  # mkuhn, 20140305: already considering read pairs
	###$desired_reads = $desired_reads/2 ; # taking half from each file.
	###$total_reads = $total_reads/2;
	
  open(TMP ,'>','temp.txt') || die $!;
	map {print TMP "$_\n"} sort {$a <=> $b}  (shuffle(1..$total_reads)) [0..($desired_reads-1)]; ## sorted list of records desired reads long from a list of numbers containing all reads.
	print TMP "-1"; # padding by one;
	print "Generated Random Numbers\n";
	close TMP; open (TMP ,'<', 'temp.txt') || die $!;
	
	
	open (my $FH1 , '<' , $pair_1) || die $!;
	open (my $FH2 , '<' , $pair_2) || die $!;
	
	my $base_name  = basename($pair_1 , qw".fastq .fq  .fq33  .fq64");
	my $base_name2  = basename($pair_2 , qw".fastq .fq  .fq33  .fq64");
	
	open (ONEOUT , '>' ,  $output_dir."/".$base_name.".subfq".$name_suffix) || die $!;
	open (TWOOUT , '>' ,  $output_dir."/".$base_name2.".subfq".$name_suffix) || die $!;

	my $fastq_iterator = iterator($FH1);
	my $fastq_iterator_2 = iterator($FH2);
	
	my $entry = $fastq_iterator->(); 
	my $entry2 = $fastq_iterator_2->();
	
	my $counter = 1;
	
	my $sample_num = readline(TMP);

	
	
	
	while ($sample_num != -1){
	
	if (($sample_num == $counter)){
		print  ONEOUT $entry;
		print  TWOOUT $entry2;
		$sample_num = <TMP>;
		}
	 $entry = $fastq_iterator->(); 
	 $entry2 = $fastq_iterator_2->();	
		
	$counter++;
	}
	
	
	
unlink 'temp.txt';
close $FH2;
close $FH1;
close ONEOUT;
close TWOOUT;

if ($interleave_out == 1){
	open (ONE , '<' , $base_name.".subfq") || die $!;
	open (TWO , '<' , $base_name2.".subfq") || die $!;
	open (MIX , '>' , $base_name."_inter.subfq") || die $!;

	select MIX;
	while (<ONE>){
		print $_;
		$_ = <ONE>;
		print $_;
		$_ = <ONE>;
		print $_;
		$_ = <ONE>;
		print $_;
		
		$_ = <TWO>;
		print $_;
		$_ = <TWO>;
		print $_;
		$_ = <TWO>;
		print $_;
		$_ = <TWO>;
		print $_;
		}
	}
}	
	
	
	
else{
	#implement single end;
  # mkuhn, 20140305
  die "Single-end not implemented." if $paired_end == 0;
	}


exit;


sub iterator{
	my $file_handle = shift // die $!;
	return sub{
		my $record .= readline($file_handle) // return;
		$record .= readline($file_handle);
		$record .= readline($file_handle);
		$record .= readline($file_handle);
		return $record;
		};
	}
	

sub print_usage{
	print "\nUsage is ./SubSample_Fastq <options> \n";
	print "Available Options:
	--total <int>\t total reads to select from [default is calculated from the file(s)]
	--desired  <int>\t desired number of reads [mkuhn: read pairs]
	--paired <1|0>\t paired end reads [default ON (1)]
	-1 <string>\t file name forward
	-2 <string>\t file name in reverse
  --output <string>\t output directory [mkuhn]
  --name <string>\t name suffix of output file
  --verbose <0|1>\t verbose output to stdout
	--interleave_out <1|0> interleave the output file\n";
	exit;	
	}
