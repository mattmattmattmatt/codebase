#!/usr/bin/perl

use strict;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;
use List::Util;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "coord=s",
	   "ref=s",
	   "all=s",
	   "fail",
	   "full",
	   "full_rev",
	   "dir=s",
	   "count",
	   "binsize=i",
	   "us",
	   "verbose",
	   "out_sort",
	   "dbg",
	   "fiftyp");

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{coord} || !$OPT{ref});

=pod

=head1 SYNOPSIS

overlapcoordinates_fast.pl will read the 2 input files (ref, coord) and check which elements in ref overlap elements in coord (-fail will report the regions that have no match)

Required flags: -ref -coord

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

    -fail  reports only non-matching elements

    -all   reports everything, whether matched or not

    -full  reports only matching that are 100% overlapping

    -full_rev reports only elements in -coord, that are fully overlapped by elements in -ref

    -dir   Will look throught this directory and overlap with every file matching "ref" pattern

    -count Will only report the number of elements in coord that overlap with each element in ref

    -binsize Will use bins of this size to split up our lookup table (default == 1000)

    -us Will replace all spaces in the strings in the coord descriptions with '_'

    -verbose Instead of joining all overlaps on one line, they get split into multiple lines

    -dbg Will display some notes to stderr while running
		
    -fiftyp Makes matches only if 50% of the ref region is overlapped by a coord region

=head1 NAME

overlapcoordinates_fast.pl

=head1 DESCRIPTION

August, 2007

This script will read each element from both the supplied ref and coord and report which, if any elements are overlapping between them

The input files should have the following format

chr start end (anything after these fields will be kept as the descriptive string to be reported upon match)

Usually, ref will be a list of genes, etc, while coord is a list of clone, for example

=head1 AUTHOR

Richard Corbett

=cut

#Load coord
open (FILE1, $OPT{coord}) || die "Can't open file $OPT{coord} \n";
my %gene_lists;
my %gene_locs;
my $cnt = 1;
my $binsize = exists($OPT{binsize}) ? $OPT{binsize} : 1000;
my $lineNum = 0;

while(<FILE1>) {
  chomp;
  my ($chr,$start,$end,$flagthis) = $_ =~ /(\w+)\s(\d+)\s(\d+)\s(.*)/;
  if($flagthis eq '') {
    ($chr,$start,$end) = $_ =~ /(\w+)\s(\d+)\s(\d+)/;
    $flagthis = $cnt++;
  }
  if(exists($OPT{us})) {
    $flagthis =~ s/ /_/g;
  }
  while(exists($gene_locs{$flagthis})) {
    $flagthis = $flagthis."-DUPKEY".$cnt;
    $cnt++;
  }
#  print "$chr $start $end $flagthis\n";

#  print "$flagthis\n";

  #save the info for this element in a hash
  $gene_locs{$flagthis}{set}{"$start-$end"} = 1;
  $gene_locs{$flagthis}{chr} = $chr;
  $gene_locs{$flagthis}{start} = $start;  #in case we want to sort the output

  #Add this element to a hash table to speed up the lookups later on
  my $start_ind = int($start/$binsize);
  my $stop_ind = int($end/$binsize);
  foreach( $start_ind-1 ..$stop_ind+1 ) {
    $gene_lists{$chr}{$_}{$flagthis} = 1;
  }

  if(exists($OPT{dbg}) && ($lineNum % 1000000) == 0) {
    print STDERR "Loading line number $lineNum\n"
  }
  $lineNum++;

}
close(FILE1);

#print Dumper %gene_locs;
#print Dumper %gene_lists;
#print "Finished Loading....\n";

if(exists($OPT{dir})) {
  my $d = $OPT{dir};
  my $ref = $OPT{ref};

  #Make sure the output directory is all set up, and cleaned out
  my $outpath = "$d/";
  #system("rm -Rf $outpath");
  system("mkdir -p $outpath");

  #Get a list of the hole files in this directory
  opendir(DIR,$d) || die "Can't open $d $!\n";
  my @goodfiles = grep {/$ref/} readdir DIR;

  foreach my $file (@goodfiles) {
    ##Load the stuff in ref 
    open (FILE2, "$d/$file") || die "Can't open file $d/$file \n";
    my $outfile = "$outpath/$file"."_ovlp";
    open (OUTFILE,">$outfile") || die "Can't open $outfile";
    while (<FILE2>) {
      chomp;
      my ($chr,$start,$end,$flagthis) = $_ =~ /(\w+)\s(\d+)\s(\d+)\s(.*)/;
      my %overlaps = getGenesInRange($chr, $start, $end);
      my $outstr = '';
      foreach my $gene (keys %overlaps) {
		$outstr = $outstr." $gene";
      }
      if(!exists($OPT{fail})){
		if(length($outstr) >0 ) {
	  		print OUTFILE "$chr $start $end $flagthis $outstr\n";
		} else {
	  		if(exists($OPT{all})) {
	    		my $non_mark = $OPT{all};
	    		print OUTFILE "$chr $start $end $flagthis $outstr $non_mark\n";
	  		}
		}
      }
      #only print the unmatched stuff if asked for
      if(exists($OPT{fail}) && length($outstr)==0){
		print OUTFILE "$chr $start $end $flagthis unmatched\n";
      }
    }
    close(OUTFILE);
  }

} else {

  ##Load the stuff in ref - usually a list of clones
  open (FILE2, $OPT{ref}) || die "Can't open file $OPT{ref} \n";
  while (<FILE2>) {
    chomp;
    my ($chr,$start,$end,$flagthis) = $_ =~ /(\w+)\s(\d+)\s(\d+)\s(.*)/;
    if($flagthis eq '') {
       ($chr,$start,$end) = $_ =~ /(\w+)\s(\d+)\s(\d+)/;
       $flagthis = '';
    }

    my %overlaps = getGenesInRange($chr, $start, $end);

    #only report the number of overlaps
    if(exists($OPT{count})) {
      print "$chr $start $end $flagthis ".scalar(keys %overlaps)."\n";
      next;
    }

    my @ovlps=();
    if(exists($OPT{out_sort})) {
      @ovlps =  sort { $gene_locs{$a}{start}  cmp $gene_locs{$b}{start} } keys %overlaps; 
    } else {
      @ovlps = keys %overlaps;
    }

    my $outstr = '';
    foreach my $gene (@ovlps) {
      	if(exists($OPT{verbose})) {
			print "$chr $start $end $flagthis $gene\n";
      	} else {
			$outstr = $outstr."^^^$gene";
      	}
    }

    next if exists($OPT{verbose});

    if(!exists($OPT{fail})){
      	if(length($outstr) >0 ) {
			print "$chr $start $end $flagthis $outstr\n";
      	} else {
			if(exists($OPT{all})) {
		  		my $non_mark = $OPT{all};
	  			print "$chr $start $end $flagthis $outstr $non_mark\n";
			}
      	}
    }

    #only print the unmatched stuff if asked for
    if(exists($OPT{fail}) && length($outstr)==0){
      print "$chr $start $end $flagthis unmatched\n";
    }
  }
}





##returns clones whose  alignments intersect with span
sub getGenesInRange {

	my $chr = $_[0];
	my $start = $_[1];
	my $end = $_[2];
	
	my %return_hash;
	
	#find all the clones that overlap with the given region
	my $strt = int($start/$binsize);
	my $stp = int($end/$binsize);
	foreach( $strt ..$stp ) {
	  	foreach my $gene (keys %{ $gene_lists{$chr}{$_} } ) {
			foreach my $coords (keys %{ $gene_locs{$gene}{set} }) {
				my ($ref_start, $ref_end) = split('-', $coords);
					
				if(exists $OPT{full}) {
					#check the "ref" (list being looped) is fully overlapped by "coord" (the reference hash)
					if ($ref_start <= $start && $ref_end >= $end) {
						$return_hash{$gene} = 1;
					}
				} elsif (exists $OPT{full_rev}) {
					#check the "coord" (the reference hash) is fully overlapped by "ref" (list being looped)
					if ($start <= $ref_start && $end >= $ref_end) {
						$return_hash{$gene} = 1;
					}
				} else {
					if (exists $OPT{fiftyp}) { #need at least half of the region to overlap
						if( $ref_start <= $end && $ref_end >= $start) {
							my $my_size = $end-$start;
							my $ref_size = $ref_end-$ref_start;
							my $left = $end-$ref_start;
							my $right = $ref_end-$start;
							my $mx = ($left > $right) ? $left : $right;
							$mx = ($mx < $ref_size) ? $mx : $ref_size;
							#print "Some Overlap: $mx $my_size ". $mx/$my_size." \n";
							if( $mx > (0.5 * $my_size)) {
								#print "Overlap\n";
								$return_hash{$gene} = 1;
							}
						}
					} elsif ($ref_start <= $end && $ref_end >= $start) {
						#Just looking for any overlap		
						$return_hash{$gene} = 1;
					}
				}
			}
	  	}
	}

  return %return_hash;
}
