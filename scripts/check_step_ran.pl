#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Cluster;
use modules::Pipeline;
use modules::PipelineSample;
use modules::Adaptors::Sample;
use modules::Adaptors::Pipeline_Step;
use modules::Adaptors::Pipeline_Step_Run;
use modules::SystemCall;
use modules::Adaptors::Run;
use Pod::Usage;
use Cwd;

use vars qw(%OPT);


GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"sample=s",
	   		"step=s",
	   		"xml=s"
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{step} || !$OPT{sample} || !$OPT{xml});

	   
=pod

=head1 SYNOPSIS

./check_step_ran.pl -step step_name -sample sample_name -xml xml_config

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

process_sample.pl -> Manage samples from config file

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./process_samples.pl -single_sample RA_single1 -add -submit

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

my $sample_name = $OPT{sample};
my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);	

my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);

my $xml_config_file = $OPT{xml};
my ($xml_full_path) = `readlink -f $xml_config_file`;
chomp $xml_full_path;
if (!-e $xml_config_file) {
	modules::Exception->throw("ERROR:  XML file $xml_config_file doesn't exist $!\n");
}

my $run_config = modules::ConfigXML->new($xml_config_file);

my $steps = $run_config->read('steps_order', 'step');



if (!defined $sample_obj) {
	modules::Exception->throw("Cannot retrieve any sample with name $sample_name");
}

my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_obj->id);
if (!defined $run_obj) {
	modules::Exception->throw("Cannot retrieve any runs for $sample_name");
}
my $run_id = $run_obj->id;

my $step_name = $OPT{step};

my $required_steps = 1;

if ($run_config->read('steps', 'step', $step_name, 'by_chr')) {
	my $pipe_config = modules::Pipeline::get_pipe_conf();
	my @chr = split(" ",$pipe_config->read($source_type,'annotation_version','chr'));
	$required_steps = @chr;
} 



my ($step_obj) = modules::Adaptors::Pipeline_Step->search(name=>$step_name);

if (!defined $step_obj) {
	modules::Exception->throw("Can't find step $step_name in db");
}
my $step_id = $step_obj->id;


my @pipeline_step_objs = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_id);
my $steps_run_count = @pipeline_step_objs;

if ($steps_run_count != $required_steps) {
	modules::Exception->throw("Step $step_name didn't run for sample $sample_name; Expecting $required_steps and got $steps_run_count steps");
}


