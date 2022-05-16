#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"bam=s",
	   	"source_name=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{bam} || !$OPT{source_name});

	   
=pod

=head1 SYNOPSIS

copy_bam.pl -bam bam_file -source_name source_name(needed_for_destdir) [options]

Required flags: -bam -source_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

copy_bam.pl -> Add lanes to existing sample

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

copy_bam.pl 

=cut

my $bam_file = $OPT{bam};
my $bam_index = $bam_file . '.bai';

if (!-e $bam_file) {
	modules::Exception->throw("ERROR: File $bam_file doesn't exist");
}

if (!-e $bam_index) {
	modules::Exception->throw("ERROR: File $bam_index doesn't exist");
}

my $source_name = $OPT{source_name};

my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source_name);

if (!defined $source_obj) {
	modules::Exception->throw("ERROR: Can't retrieve source obj from the db for $source_name");
}

my $source_type = $source_obj->source_type;

my $cluster_config = modules::Pipeline::get_cluster_conf();
my $destdir = $cluster_config->read($source_type,'base_directories','base_results_directory') . '/' . $source_name;

#Add the source_name directory if required
mkdir($destdir) if !-d $destdir;

my $command = "cp $bam_file $bam_index $destdir";
print STDERR "$command\n";
my $sys_call = modules::SystemCall->new();
$sys_call->run($command);

my $mdss = $cluster_config->read('common','mdss','mdss_flag');
if ($mdss) {
	my $mdss_results_dir = $cluster_config->read('common','mdss','mdss_results');
	my $mdss_mkdir_cmd = "mdss mkdir $mdss_results_dir/$source_type/$source_name";
	$sys_call->run($mdss_mkdir_cmd);
	my $mdss_command = "mdss put $bam_file $bam_index $mdss_results_dir/$source_type/$source_name";
	$sys_call->run($mdss_command);  
}








