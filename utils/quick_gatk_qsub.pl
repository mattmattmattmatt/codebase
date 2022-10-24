#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use modules::Cluster;
use modules::Adaptors::Release_File;
use Pod::Usage;
use Cwd 'abs_path';
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"readdir=s",
		"outdir=s",
		"sample_list=s",
		"templates=s",
		"qsubdir=s",
		"softdir=s",
		"run"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{readdir} || !$OPT{outdir} || !$OPT{sample_list} || !$OPT{templates});

	   
=pod

=head1 SYNOPSIS

quick_gatk_qsub.pl -template template_qsub -qsubdir qsub_dir -readdir readdir(fastq_must_be_named_'sample_R[12]_fastq.gz') -outdir outdir -sample_list sample_list_file -run submit_single_jobs(except_joint.qsub) -softdir software_directory

Required flags: -readdir -outdir -sample_list

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

quick_gatk_qsub.pl -> Create gvcf for fastq

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

scp_results.pl 

=cut

my @templates = split(",",$OPT{templates});

for my $template ( @templates ) {
	if (!-e $template) {
	        modules::Exception->throw("Template $template doesn't exist");
	}
}

my $sample_file = $OPT{sample_list};

my $readdir = $OPT{readdir};
$readdir = abs_path($readdir);

if (!-d $readdir) {
        modules::Exception->throw("Dir $readdir doesn't exist");
}

my $outdir = $OPT{outdir};
$outdir = abs_path($outdir);

if (!-d $outdir) {
        modules::Exception->throw("Dir $outdir doesn't exist");
}

my $qsubdir = defined $OPT{qsubdir}?$OPT{qsubdir}:$outdir;
$qsubdir = abs_path($qsubdir);


my @samples = ();
open(SAMPLE, $sample_file) || modules::Exception->throw("File $sample_file doesn't exist");

while (<SAMPLE>) {
	chomp;
	my ($sample) = $_ =~ /(\S+)/;
	push @samples, $sample;
}


my $softdir = defined $OPT{softdir}?$OPT{softdir}:'/g/data/u86/software/';
$softdir = abs_path($softdir);

if (!-d $softdir) {
        modules::Exception->throw("Dir $softdir doesn't exist");
}


my @g_vcfs = ();

#print Dumper \@samples;


for my $sample (@samples) {
	my $template_count = 1;
	for my $template (@templates) {
		my $qsub = "$sample.${template_count}.qsub";
		open(QSUB, ">$qsubdir/$qsub") || modules::Exception->throw("Can't open qsub for $sample");
		open(TEMPLATE,$template) || modules::Exception->throw("Can't open file $template");
		
		my $next_qsub;
		if ($qsub =~ /(\d+)\.qsub/) {
			my $qsub_count = $1 + 1;
			
			($next_qsub = $qsub) =~ s/\d+.qsub/$qsub_count.qsub/;	
		}
		
		
		
		while (<TEMPLATE>) {
			$_ =~ s/QSUBDIR/$qsubdir/g;
			$_ =~ s/SAMPLE/$sample/g;
			$_ =~ s/OUTDIR/$outdir/g;
			$_ =~ s/READDIR/$readdir/g;
			$_ =~ s/QSUBNEXT/$next_qsub/g;
			$_ =~ s/SOFT/$softdir/g;
			
			print QSUB $_;
			if ($_ =~ /(\S+g.vcf.gz$)/) {
				push @g_vcfs, $1;
			}
		}	
		close QSUB;
		$template_count++;		
	}
	
	
	
}

open(JOINT,">$qsubdir/joint_call.qsub") || modules::Exception->throw("Can't open file $qsubdir/joint_call.qsub\n");
print JOINT "#!/bin/bash\n#PBS -P pq84\n#PBS -q normal\n#PBS -l walltime=48:00:00,mem=64GB,jobfs=100GB,ncpus=4\n#PBS -l storage=gdata/u86+scratch/u86+gdata/xx92+gdata/pw1\n";
print JOINT "\nmodule load java/jdk-8.40\n";
my $gvf_str = join(' -V ',@g_vcfs);
print JOINT "java -Xms16g -Xmx16g -jar /g/data/u86/variantdb/v2.3/ext/bin/GATK-3.6.0/GenomeAnalysisTK.jar -T GenotypeGVCFs -nt 1 -R /g/data/u86/variantdb/v2.38/conf/human/GRCh38/fasta/single_file/GRCh38d1_noalt.fa -V $gvf_str -o $outdir/joint_calls.vcf\n";






