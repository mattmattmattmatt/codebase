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
     	"runid=i",
     	"run"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid});

	   
=pod

=head1 SYNOPSIS

cleanup.pl -runid runid -run run_commands(print_otherwise) [options]

Required flags: -runid

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

cleanup.pl -> remove large files not needed

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

copy_bam.pl 

=cut

my $runid = $OPT{runid};
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);

my $gatk_directory  = $run_dir .'/gatk';
my $bam_directory = $run_dir . '/bam';
my @commands;

my $sys_call = modules::SystemCall->new();

opendir(BAM,$bam_directory) || modules::Exception->throw("Can't open directory $bam_directory");
my @bam_files = grep {/ba[im]/} readdir BAM;

#Check the merge gatk file is ready
my @merge_gatk = grep {/gatk/} @bam_files; 

if (!@merge_gatk) {
	modules::Exception->throw("ERROR: Can't find a merged GATK file in $gatk_directory");
}

my @markdups = grep {/dupmark/} @bam_files;
for my $bam_file ( @markdups ) {
    push @commands, "rm $bam_directory/$bam_file";
}

my @mergebams = grep {/merge_bam/} @bam_files;
for my $bam_file ( @mergebams ) {
    push @commands, "rm $bam_directory/$bam_file";
}

opendir(GATK,$gatk_directory) || modules::Exception->throw("Can't open directory $gatk_directory");
my @bam_gatk_files = grep {/\.*.ba[im]/} readdir GATK;
for my $gatk_file ( @bam_gatk_files ) {
    push @commands, "rm $gatk_directory/$gatk_file";
}

for my $command (@commands) {
	print STDERR "$command\n";
	$sys_call->run($command) if $OPT{run};	
}










