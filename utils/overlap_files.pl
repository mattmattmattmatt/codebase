#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Overlap;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "ref=s",
	   "coord=s",
	   "ignore=s",
	   "exact",
	   "fail",
	   "silent",
	   "fiftyp",
	   "max_distance=i",
	   "full",
	   "full_rev",
	   "all",
	   "just_overlap",
	   "debug"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{coord} || !$OPT{ref});



=pod

=head1 SYNOPSIS

overlap.pl -ref ref_file -coord coord_file -ignore ignore_coord_entries [options]

Required flags: -ref -coord

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

overlap.pl -> Test the Overlap module for overlapping coordinate sets

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./overlap.pl -ref list1 -coord list2

=cut

my $overlap = modules::Overlap->new();
my $ref = $OPT{ref};
my $coord = $OPT{coord};

if (!-e $ref || !-e $coord) {
    print "Problem with input files $!\n";
    exit;
}

my $overlap_struct;
my $count;

my %args = (
			-ref=>$ref,
			-coord=>$coord
			);

if ($OPT{ignore}) {
	$args{-ignore} = $OPT{ignore};
} 

if ($OPT{max_distance}) {
	$args{-max_distance} = $OPT{max_distance};
}

if ($OPT{full}) {
	$args{-full} = 1;
} 

if ($OPT{full_rev}) {
	$args{-full_rev} = 1;
} 

if ($OPT{all}) {
	$args{-all} = 1;
} 

if ($OPT{fail}) {
	$args{-fail} = 1;
} 

if ($OPT{exact}) {
	$args{-exact} = 1;
} 

if($OPT{silent}) {
	$args{-silent} = 1
} 

if ($OPT{fiftyp}) {
	$args{-fiftyp} = 1;
}

($overlap_struct, $count) = $overlap->overlap(%args);

print Dumper $overlap_struct if $OPT{debug};
print "OVERLAP COUNT: $count\n" unless $OPT{just_overlap};

