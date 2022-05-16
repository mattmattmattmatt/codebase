#! /usr/bin/perl -w

use strict;
use modules::Adaptors::Sample;
use modules::Exception;
use modules::Pipeline;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"sample_type=s"
    );

	   
pod2usage(1) if ($OPT{help});


=pod

=head1 SYNOPSIS

sampole_lists.pl -sample_type get_specific_sample_list(one_of_ENUM_values;default=ALL)

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

sample_lists.pl -> Get sample group files from the database

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./sample_list.pl -sample_type normal_cell_line

=cut

my $report_conf = modules::Pipeline->get_report_conf();

my %sample_types = map {$_=>1} split(",",$report_conf->read('common','sample_types'));


my $sample_type_flag = 0;

if (defined $OPT{sample_type}) {
	if (!exists $sample_types{$OPT{sample_type}}) {
		modules::Exception->throw("ERROR: sample type must be one ENUM values ". join(' ',keys %sample_types). "\n");
	} 
	$sample_type_flag = 1;
} 

my @samples = modules::Adaptors::Sample->search_all();
my %sample_type_types = ();
for my $sample ( @samples ) {
    my $sample_name = $sample->sample_name;
    my $sample_type = $sample->sample_type;
    push @{$sample_type_types{$sample_type}}, $sample_name;
}

for my $sample_type_type ( keys %sample_type_types ) {
	if ($sample_type_flag) {
		next unless $sample_type_type eq $OPT{sample_type};
	} 
	my $file = $sample_type_type . '_list';
    open(FILE,">$file") || modules::Exception->throw("Can't open file $file\n");
    print FILE join("\n",sort @{$sample_type_types{$sample_type_type}}) ."\n";
    close FILE;
}
	


