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
		"tumor_list=s",
		"outdir=s",
		"prefix=s",
		"normal_list=s",
		"template=s",
		"qsubdir=s",
		"software=s",
		"bychr"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{normal_list} || !$OPT{outdir} || !$OPT{tumor_list} || !$OPT{template});

	   
=pod

=head1 SYNOPSIS

quick_mutect4.pl -template template_qsub -qsubdir qsub_dir -outdir outdir -tumor_list -normal_list normal_list -bychr 

Required flags: -outdir -tumor_list -normal_list -template

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

my $template = $OPT{template};

if (!-e $template) {
        modules::Exception->throw("Template $template doesn't exist");
}

my $tumor_file = $OPT{tumor_list};
my $normal_file = $OPT{normal_list};

my $outdir = $OPT{outdir};
$outdir = abs_path($outdir);

if (!-d $outdir) {
        modules::Exception->throw("Dir $outdir doesn't exist");
}

my $qsubdir = defined $OPT{qsubdir}?$OPT{qsubdir}:$outdir;
$qsubdir = abs_path($qsubdir);

my $pref = defined $OPT{pref}?$OPT{pref}:'mutect';
my $bychr = defined $OPT{bychr}?1:0;

my @samples = ();
my @normals = ();
open(TUMOR, $tumor_file) || modules::Exception->throw("File $tumor_file doesn't exist");

while (<TUMOR>) {
	chomp;
	my ($sample) = $_ =~ /(\S+)/;
	if ($sample =~ /l1$/) {
		$sample =~ s/_l1//;
	}
	push @samples, $sample;
}

open(NORMAL, $normal_file) || modules::Exception->throw("File $normal_file doesn't exist");


while (<NORMAL>) {
	chomp;
	my ($sample) = $_ =~ /(\S+)/;
	
	if ($sample =~ /l1$/) {
		push @normals, $sample;
		$sample =~ s/_l1//;
	} else {
		push @normals, $sample."_l1";
	}
	
	push @samples, $sample;
}






my $softdir = defined $OPT{softdir}?$OPT{softdir}:'/g/data/u86/software/';
$softdir = abs_path($softdir);

if (!-d $softdir) {
        #modules::Exception->throw("Dir $softdir doesn't exist");
}



#print Dumper \@samples;

		
				
if ($bychr) {
	my @chr = (1..22,'X');
	for my $chr ( @chr ) {
		my $qsub_chr = "${pref}.${chr}.qsub";
		open(QSUB,">$qsubdir/$qsub_chr") || modules::Exception->throw("Can't open file $qsub_chr\n");
		
		open(TEMPLATE,$template) || modules::Exception->throw("Can't open file $template");
		
		
		while (<TEMPLATE>) {
			$_ =~ s/QSUBDIR/$qsubdir/g;
			$_ =~ s/OUTDIR/$outdir/g;
			$_ =~ s/SOFT/$softdir/g;
			$_ =~ s/CHR/$chr/g;
			$_ =~ s/PREF/$pref/g;
			my $sample_bams;
			my $normal_samples;
			
			for my $sample (@samples) {
				$sample_bams .= " -I ${outdir}/${sample}_md_realign_recal.bam";
			}
			
			for my $normal (@normals) {
				$normal_samples .= " -normal $normal";
			}
			
			
			$_ =~ s/SAMPLE_BAMS/$sample_bams/;
			$_ =~ s/NORMAL_SAMPLES/$normal_samples/;
			$_ =~ s/BYCH/-L \/g\/data\/u86\/variantdb\/v2.38\/conf\/human\/GRCh38\/fasta\/intervals\/${chr}.intervals/;
			print QSUB $_;
		}
		
		
		
		
		close QSUB;
	}
} else {
	my $qsub = "${pref}.qsub";
	open(QSUB, ">$qsubdir/$qsub") || modules::Exception->throw("Can't open qsub");
	open(TEMPLATE,$template) || modules::Exception->throw("Can't open file $template");
	while (<TEMPLATE>) {
		$_ =~ s/QSUBDIR/$qsubdir/g;
		$_ =~ s/OUTDIR/$outdir/g;
		$_ =~ s/SOFT/$softdir/g;
		$_ =~ s/\.CHR//g;
		$_ =~ s/BYCH //;
		$_ =~ s/PREF/$pref/g;
		my $sample_bams;
		my $normal_samples;
		
		for my $sample (@samples) {
			$sample_bams .= " -I ${outdir}/${sample}_md_realign_recal.bam";
		}
		
		for my $normal (@normals) {
			$normal_samples .= " -normal $normal";
		}
		
		
		$_ =~ s/SAMPLE_BAMS/$sample_bams/;
		$_ =~ s/NORMAL_SAMPLES/$normal_samples/;
		print QSUB $_;
		
	}
	close QSUB;
}
		








