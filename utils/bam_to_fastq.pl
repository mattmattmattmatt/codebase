#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use File::Basename;
use modules::SystemCall;
use Cwd 'abs_path';
use File::Basename;


GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"bam=s",
		"fq1=s",
		"fq2=s",
		"outdir=s",
		"ref=s",
		"run"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{bam});

	   
=pod

=head1 SYNOPSIS

./bam_to_fastq.pl -bam bam_file -fq1 fastq1_name(default=bam.R1.fq) -fq2 fastq2_name(default=bam.R2.fq) -outdir output_directory(default=current) -ref reference_genome(default=GRCh38_NCI)[options]

Required flags: -bam

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

bam_to_fastq.pl -> Convert BAM file to fastq to allow consistent analysisSychronise v2.1 and v2.2

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

=cut

#Get binaries
my $base = abs_path($0);


#Get outdir
my $outdir = defined $OPT{outdir}?`readlink -f $OPT{outdir}`:`pwd`;
chomp $outdir;

my $ref = defined $OPT{ref}?$OPT{ref}:'/g/data/u86/variantdb/v2.38/conf/human/GRCh38/fasta/single_file/GRCh38d1_noalt.fa';

#test here
my $bam = $OPT{bam};
my $short_bam = basename($bam);

if (!-e $bam) {
	modules::Exception->throw("ERROR: Bam file $bam doesn't exist");
}

#Make sure BAM is full path
my $bam_full = abs_path($bam);

my $fq1 = defined $OPT{fq1}?$OPT{fq1}:$short_bam.".R1.fq";
my $fq2 = defined $OPT{fq2}?$OPT{fq2}:$short_bam.".R2.fq";

my @commands = ();
#assumes 48 threads
push @commands, "module load samtools";
push @commands, "samtools  sort -n -@ 48 JT_product.RG.bam| samtools bam2fq -@ 48 -1 $outdir/$fq1  -2 $outdir/$fq2 --reference $ref  -s /dev/null -0 /dev/null -";


my $sys_call = modules::SystemCall->new();

for my $cmd (@commands) {
   print "$cmd\n";
   $sys_call->run($cmd) if $OPT{run}; 
}










