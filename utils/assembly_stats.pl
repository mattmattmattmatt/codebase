#! /usr/bin/perl
use Statistics::Descriptive::Full;
use Data::Dumper;
use strict;

my $length;
my $NfreeLength;
my $totalLength;
my $totalNFreeLength;
my $k_contigs = my $oneM_contigs = my $tenM_contigs = 0;
my @arr;

while(<>){
    chomp;
    next if /^$/;
    if(/>/){
      next if $length !~ /\d/; #Account for first contig
		  push (@arr, $length);
		  $totalLength += $length;
      $totalNFreeLength += $NfreeLength;   

      if ($length > 100000) {
			  $k_contigs++;
		  }
		  if ($length > 1000000) {
			  $oneM_contigs++;
		  } 
		  if ($length > 10000000) {
			  $tenM_contigs++;
		  }
		  $length=0;
      $NfreeLength = 0;
      next;
    }
    
    $length += length($_);
    $_ =~ s/N//gi;
    $NfreeLength += length($_);
    #print "L $length NFL $NfreeLength\n";
}
#Account for last contig
push (@arr, $length);
$totalLength += $length;
$totalNFreeLength += $NfreeLength;

my ($contig_number) = @arr;

close(FH);

my @sort = sort {$b <=> $a} @arr;
my $n50_count;
my $L50;
my $N50;

foreach my $val(@sort){
    #print "$val\n";
    $n50_count+=$val;
    if($n50_count >= $totalLength/2){
    	$N50 = $val;
    	#print "Total: $totalLength N50 length is $n50 and N50 value is: $val and L50 is $L50\n";
        last;
    }
    $L50++;
}

my $stats = Statistics::Descriptive::Full->new();
$stats->add_data(@arr);
printf("Contig_num: %d\nAssembly_size: %d\nN-free_size: %d\nLongest_contig: %d\nContigs_100k+: %d\nContigs_1M+: %d\nContigs_10M+: %d\nN50: %d\nL50: %d\n",
       $stats->count,
       $totalLength,
       $totalNFreeLength,
       $stats->max,
	   $k_contigs,
	   $oneM_contigs,
	   $tenM_contigs,
	   $N50,
	   $L50
      );
