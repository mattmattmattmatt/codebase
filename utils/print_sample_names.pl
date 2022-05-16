#! /usr/bin/perl -w

use strict;
use modules::Exception;
use modules::Pipeline;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "single_id=s",
		   "cohort_id=s",
		   "project_name=s"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{cohort_id} & !$OPT{project_name} &!$OPT{single_id});

	   
=pod

=head1 SYNOPSIS

print_sample_names.pl -cohort_id cohort_id -project_name project_name [options]

Required flags: -cohort_id or -project_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

print_sample_names.pl -> get sample names from db

=head1 DESCRIPTION

Sep 10, 2014

a script that ...

=head1 AUTHOR

Vicky Cho

=head1 EXAMPLE
./print_sample_names.pl -single_id APOSLE_single88
./print_sample_names.pl -cohort_id RARE_cohort2
./print_sample_names.pl -project_name APOSLE

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my @sample_names;
my $outfile;

if($OPT{single_id}){
	my $single_id = $OPT{single_id};
	my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$single_id);
	if(!defined $source_obj){
		modules::Exception->throw("ERROR:Can't retrieve source from db for $single_id");
	}	
	my ($sg_obj) = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
	my ($sample_obj) = modules::Adaptors::Sample ->search(source_group_id=>$sg_obj->id);
	push @sample_names, $sample_obj->sample_name;

	$outfile = "/g/data/u86/snv_pipeline_runs/v2.1_human_genome/sample_names/$single_id\_sample_names.txt";
}

if($OPT{cohort_id}){
	my $cohort_id = $OPT{cohort_id};
	my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$cohort_id);
	if(!defined $source_obj){
		modules::Exception->throw("ERROR:Can't retrieve source from db for $cohort_id");
	}
	my ($sg_obj) = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
	my (@sample_objs) = modules::Adaptors::Sample->search(source_group_id=>$sg_obj->id);
	
	for my $sample (@sample_objs){
		push @sample_names, $sample->sample_name;
	}

	$outfile = "/g/data/u86/snv_pipeline_runs/v2.1_human_genome/sample_names/$cohort_id\_sample_names.txt";

}

if($OPT{project_name}){
	my $project_name = $OPT{project_name};
	my @sources;
	my @source_groups;

	my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
	if (!defined $project_obj) {
		modules::Exception->throw("ERROR: Can't retrieve project from db for $OPT{project}");
	}
	my @source_objs = modules::Adaptors::Source->search(project_id=>$project_obj->id);
	if (!@source_objs) {
		modules::Exception->throw("ERROR: Couldn't retrieve any sources for project $OPT{project}");
	}		
	
	for my $source_obj (@source_objs) {
		my ($sg_obj) = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
		my (@sample_objs) = modules::Adaptors::Sample->search(source_group_id=>$sg_obj->id);
		
		for my $sample (@sample_objs){
			push @sample_names, $sample->sample_name;
		}
		
	}

	$outfile = "/g/data/u86/snv_pipeline_runs/v2.1_human_genome/sample_names/$project_name\_sample_names.txt";
}

open(OUTFILE,">$outfile") || modules::Exception->throw("Can't open file to write $outfile\n");

for my $sample (@sample_names){
	print OUTFILE "$sample\n";
}

close OUTFILE;
