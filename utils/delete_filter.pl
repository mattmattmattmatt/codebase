#! /usr/bin/perl -w

use strict;
use modules::Adaptors::VariantDB;
use DBI;
use modules::Exception;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "run_id=i",
	   "filter_name=s",
	   "indel"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{run_id} || !$OPT{filter_name});

=pod

=head1 SYNOPSIS

delete_filter.pl -run_id run_id -filter_name filter_name_to_delete -indel delete_indel_filters [options]

Required flags: -run_id -filter_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_delete.pl -> delete snv_filter or variant_filter for a specfic filter

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./delete_filter.pl -run_id 10 -filter_name filter_dbsnp_snv 
./delete_filter.pl -run_id 10 -filter_name filter_dbsnp_indel -indel


=cut

my $run_id = $OPT{run_id};
my $delete_filter_name = $OPT{filter_name};
my $indel = defined $OPT{indel}?1:0;

# Connect to database

my $dbh = DBI->connect(modules::DBAdaptor::VariantDB->getConfig) 
    or die "Unable to connect: $DBI::errstr\n";

# Retrieve filter id

my $filterid_sth = $dbh->prepare("SELECT id FROM filters WHERE name = ?");

$filterid_sth->execute($delete_filter_name);

my $filter_id;

if ($filterid_sth->rows == 1){
    ($filter_id) = $filterid_sth->fetchrow_array;
} else {
    modules::Exception->throw("Fetching id for filter [$delete_filter_name] returned " 
	. $filterid_sth->rows .  " rows");
}

if (!$indel) {
	# Search for snv ids directly with DBI
	my $snvid_sth = $dbh->prepare("SELECT id FROM snvs WHERE run_id = ?");
	$snvid_sth->execute($run_id);
	my @snv_ids;
	
	while (my ($snv_id) = $snvid_sth->fetchrow_array){
	    push @snv_ids, $snv_id;
	}
	
	# Delete snv_filter rows
	my $sf_delete_sth = $dbh->prepare('DELETE FROM snv_filters WHERE snv_id IN (' . join(',', @snv_ids) . ') AND filter_id = ' . $filter_id . ";\n");
	$sf_delete_sth->execute;
} else {
	# Search for snv ids directly with DBI
	my $variantid_sth = $dbh->prepare("SELECT id FROM variants WHERE run_id = ?");
	$variantid_sth->execute($run_id);
	my @variant_ids;
	
	while (my ($variant_id) = $variantid_sth->fetchrow_array){
	    push @variant_ids, $variant_id;
	}
	
	# Delete variant_filter rows
	my $sf_delete_sth = $dbh->prepare('DELETE FROM variant_filters WHERE variant_id IN (' . join(',', @variant_ids) . ') AND filter_id = ' . $filter_id . ";\n");
	$sf_delete_sth->execute;
}
