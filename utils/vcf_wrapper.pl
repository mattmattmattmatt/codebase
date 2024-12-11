#! /usr/bin/perl -w

use strict;
use FindBin qw( $RealBin );
use lib $RealBin;
use modules::Exception;
use modules::Vcf;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;
use File::Basename;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"vcf_in=s",
		"vcf_out=s"
    	);

	   
pod2usage(1) if ($OPT{help} || !$OPT{vcf_in});


=pod

=head1 SYNOPSIS

./vcf_wrapper.pl -vcf_in -vcf_out vcf_out(default=replace_vcf_suffix_with_txt)

Required flags: -vcf_in

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

vep_wrapper.pl -> Get vep info from variant outside pipeline

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./sample_list.pl -sample_type normal_cell_line

=cut

my $vcf_file = $OPT{vcf_in};

if ( !-e $vcf_file ) {
	modules::Exception->throw("File $vcf_file doesn't exist");
}

my $vcf = Vcf->new(-vcf=>$vcf_file);

#For output files
my $vcf_out;
if (defined $OPT{vcf_out}) {
	$vcf_out = basename($OPT{vcf_out});
} else {
	my ($vcf_short) = basename($vcf_file);
	($vcf_out = $vcf_short) =~ s/.vcf/.txt/;
}

my $dir = dirname($vcf_file);
$vcf->check_vcf(-vcf_file=>$vcf_file);
$vcf->parse_vcf(-vcf_file=>$vcf_file);
$vcf->write_normalised(-vcf_file=>$vcf_file,-vcf_out=>"$dir/$vcf_out");
