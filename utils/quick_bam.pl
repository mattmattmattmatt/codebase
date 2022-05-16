#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT $PBS_JOBFS @RG);
use File::Basename;
use modules::SystemCall;
use Cwd 'abs_path';
use File::Basename;
use modules::Exception;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"bam_name=s",
		"fq1=s",
		"fq2=s",
		"outdir=s",
		"run",
		"ref=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{fq1});

	   
=pod

=head1 SYNOPSIS

./quick_bam.pl -bam_name bam_basename(default=fq1_Without_R1) -fq1 fastq1_name -fq2 fastq2_name(default=swap_R1_with_R2) -run -outdir output_directory(default=current) -ref genome_ref_name(default=hs37d5)[options]

Required flags: -fq1 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

quick_bam.pl -> Generate quick BAM file on raijin outside pipe

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

=cut

#Get binaries
my $base = abs_path($0);
(my $software_bin_full = $base) =~ s/utils\/quick_bam.pl/ext\/bin/;
my $ref = defined $OPT{ref}?$OPT{ref}:'hs37d5';

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $bwa_index = "${svndir}/conf/human/${ref}/bwa_index/${ref}";

my $samtools_bin = $software_bin_full . '/samtools';
my $bwa_bin = $software_bin_full . '/bwa_mem';
my $novosort_bin = $software_bin_full . '/novosort';

#Get outdir
my $outdir = defined $OPT{outdir}?`readlink -f $OPT{outdir}`:`pwd`;
chomp $outdir;

my $fq1 = $OPT{fq1};
my $fq2;
if ($OPT{fq2}) {
	$fq2 = $OPT{fq2};
} else {
	($fq2 = $fq1) =~ s/R1/R2/;
}


if (!-e $fq1 || !-e $fq2) {
	modules::Exception->throw("ERROR: fastq file problems R1 $fq1 R2 $fq2");
}


my $bam_name;

if (defined $OPT{bam_name}) {
	$bam_name = $OPT{bam_name};
} else {
	($bam_name = $fq1) =~ s/R1.*//;
	$bam_name =~ s/_$//;
} 


my $qsub = $outdir .'/' . basename($bam_name) . '.qsub';
open(QSUB,">$qsub") || modules::Exception->throw("Can't open $qsub");

print "QSUB $qsub BAM $bam_name\n";

print QSUB <<EOF;
#!/bin/bash
#PBS -P u86
#PBS -q normal
#PBS -l walltime=48:00:00,mem=32GB,jobfs=200Gb,ncpus=16
#PBS -l other=gdata2
#PBS -W umask=0007
#PBS -W group_list=u86
module load samtools
EOF

my $align_command = "$bwa_bin mem -M -R \"\@RG\\tID:${bam_name}_l1\\tSM:${bam_name}\\tPL:ILLUMINA\" -t 16 $bwa_index $fq1 $fq2 |  $samtools_bin view -u -S - | $novosort_bin --md --kt -c 8 -x 6 -m 10G -t \$PBS_JOBFS -i -o ${bam_name}.bam -";
print QSUB $align_command ."\n";
close QSUB;

system("qsub $qsub") if $OPT{run};









