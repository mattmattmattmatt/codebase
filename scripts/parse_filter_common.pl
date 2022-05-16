#! /usr/bin/perl -wall

use strict;
use modules::Adaptors::SNV;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;
use modules::Adaptors::Filter;
use modules::Adaptors::Run;
use modules::Exception;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Benchmark;
use vars qw(%OPT);
use Cwd; 

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "outdir=s",
		   "ref=s",
		   "source_type=s"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

parse_filter_common.pl -outdir output_directory -ref reference_genome(used_in_name) -source_type source_type(default=mouse_single) [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_filter_common.pl -> Script to generate list of commonly occuring variants

=head1 DESCRIPTION

March 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE


=cut

#only used currently for mouse
my $ref = defined $OPT{ref}?$OPT{ref}:"mm10";
my $outdir = defined $OPT{ourdir}?$OPT{outdir}:'.';
my $source_type = defined $OPT{source_type}?$OPT{source_type}:'mouse_single';
# Retrieve previous experiments from database

my @runids = modules::Adaptors::Run->search(production=>1);
my ($filter_snv) = modules::Adaptors::Filter->search(name => 'filter_snv');
my $filterid = $filter_snv->id;
my %all_snvs = ();
my %all_indels = ();
my %sample_names = ();

for my $run_obj ( @runids ) {
	my $runid = $run_obj->id;
	#check the sample is the correct organism
	my ($sample_obj) = modules::Adaptors::Sample->search(sample_id => $run_obj->sample_id);
	my $sample_name = $sample_obj->sample_name;

	push @{$sample_names{$sample_obj->sample_name}}, $runid;

	#Check we don't have multiple production runs for a single igl id
	if (@{$sample_names{$sample_obj->sample_name}} > 1) {
		my $run_str = join(",",@{$sample_names{$sample_obj->sample_name}});
		modules::Exception->throw("Sample lib $sample_name has multiple production runs $run_str");			
	}
	
	#Retrieve snvs
	
	print STDERR "Retrieving snps for run: $runid igl: $sample_name";
	my @source_snps = modules::Adaptors::SNV->search_by_experiment_and_passed_filter($run_obj->run_id, $filter_snv->id);
	

	my $snp_number = @source_snps;
	print STDERR "Found $snp_number snps";

	if ($snp_number == 0) {
		print STDERR "WARNING: run: $runid igl: $sample_name has no SNVs";
	}
    
	for my $snv_obj ( @source_snps ) {
	   	$all_snvs{$snv_obj->chr}{$snv_obj->coord}{$sample_name}++;
	}
	
	#Retrieve indels
	my @variant_types = qw(filter_insertion filter_deletion);
	for my $variant_type (@variant_types) {
		my @source_indels = modules::Adaptors::Variant->search_by_experiment_and_passed_filter_and_type($run_obj->run_id, $filterid,$variant_type);
		my $indel_number = @source_indels;
		print STDERR "Found $indel_number snps";

		if ($indel_number == 0) {
			print STDERR "WARNING: run: $runid igl: $sample_name has no INDELs";
		}
		for my $indel_obj ( @source_indels ) {
			$all_indels{$indel_obj->chr}{$indel_obj->start_coord}{$indel_obj->end_coord}{$sample_name}++;
		}
		
	}
	
	  
}

for my $chr ( sort keys %all_snvs ) {
	my $chr_output = $outdir.'/'.$ref."common.overlap.snv.".$chr;
	open(FILE,">$chr_output") ||  modules::Exception->throw("Can't open file $chr_output\n");
    for my $coord (sort {$a<=>$b} keys %{$all_snvs{$chr}}) {
    	my $run_str = join(",",keys %{$all_snvs{$chr}{$coord}});
    	print FILE join("\t",
    					$chr,
    					$coord,
    					$coord,
    					$run_str,
    					);
    }
    close FILE;
}


for my $chr ( sort keys %all_indels ) {
	my $chr_output = $outdir.'/'.$ref."common.overlap.indel.".$chr;
	open(FILE,">$chr_output") ||  modules::Exception->throw("Can't open file $chr_output\n");
    for my $start_coord (sort {$a<=>$b} keys %{$all_indels{$chr}}) {
    	for my $end_coord (sort {$a<=>$b} keys %{$all_indels{$chr}{$start_coord}}) {
    	
    		my $run_str = join(",",keys %{$all_indels{$chr}{$start_coord}{$end_coord}});
    		print FILE join("\t",
    					$chr,
    					$start_coord,
    					$end_coord,
    					$run_str,
    					);
    	}
    }
    close FILE;
}





