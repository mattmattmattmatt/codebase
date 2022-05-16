#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);
use Bio::SearchIO;


GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vartrix=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vartrix});



=pod

=head1 SYNOPSIS

vartrix_count.pl -vatrix snv_matrix.tsv(from_R) [options]

Required flags: -vartrix

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

template.pl -> Script to do stuff

=head1 DESCRIPTION

Feb 5th, 2015

This script ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

>template.pl -arg1 arg1 -arg2 arg2

=cut

my $snv_matrix = defined $OPT{vartrix}?$OPT{vartrix}:'./snv_matrix.tsv';

open(MATRIX,$snv_matrix) || modules::Exception->throw("Can't open file $\n");

my %map = (
			"0" => "No Call",
			"1" => "Ref",
			"2" => "Hom",
			"3" => "Het"
				);

print join("\t",
			"Coord",
			"No data",
			"Ref",
			"Hom",
			"Het",
			"% Calls",
			"% Variant"
			) ."\n\n";

while (<MATRIX>) {
	next unless /chr/;
	chomp;
	my @fields = split("\t");
	my $coord = shift @fields;
	$coord =~ s/"//g;
	my %line_count = ();
	for my $call (@fields) {
		$line_count{$map{$call}}++;
	}
	
	my $nd = defined $line_count{'No Call'}?$line_count{'No Call'}:0;
	my $ref = defined $line_count{'Ref'}?$line_count{'Ref'}:0;
	my $hom = defined $line_count{'Hom'}?$line_count{'Hom'}:0;
	my $het = defined $line_count{'Het'}?$line_count{'Het'}:0;
	my $var_sum = $het+$hom;
	my $called_sum = $var_sum + $ref;
	
	
	my $pc_total = sprintf("%.2f",$var_sum/4958 *100);
	my $pc_called = $called_sum>0?sprintf("%.2f",$var_sum/$called_sum *100):0.00;
	
	print join("\t", "$coord", $nd, $ref, $hom, $het, $pc_called, $pc_total) . "\n";
	#print Dumper \%line_count;
}