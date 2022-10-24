#!/usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
#use modules::Pipeline;
use modules::Cluster;
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
			"template=s",
			"qsubdir=s",
			"read_length=i",
			"notrim",
			"run",
			"ref=s",
			"rlib=s"
		   );
		   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{readdir} || !$OPT{outdir} || !$OPT{sample_list});

	   
=pod

=head1 SYNOPSIS

quick_RNASeq_qsub.pl -template template_qsub(default=RNASEQ_template_trim.qsub in cwd) -ref ref_genome(default=GRCh38) -notrim skip_trim_step -read_length read_length(default=150bp) -qsubdir qsub_dir(default=outdir) -readdir readdir(fastq_must_be_named_'sample_R[12]_fastq.gz') -outdir outdir -sample_list sample_list_file -run submit_single_jobs(except_joint.qsub) -rlib path_to_use

Required flags: -readdir -outdir -sample_list

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

quick_gatk_qsub.pl -> Create RNASeq jobs for trimming, STAR aligner, and consensusDE

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

scp_results.pl 

=cut



my $template;

my $refdir = defined $OPT{refdir}?$OPT{refdir}:'/g/data/u86/mxf221/ref_genomes/';

if ($OPT{template}) {
	$template = $OPT{template};
} else {
	if ($OPT{notrim}) {
		$template = '../templates/RNA_template_notrim.qsub';
	} else {
		$template = '../templates/RNA_template_trim.qsub';
	}
}

if (!-e $template) {
        modules::Exception->throw("Template $template doesn't exist");
}

my $ref = defined $OPT{ref}?$OPT{ref}:'GRCh38';

my $rlib = defined $OPT{rlib}?$OPT{rlib}:'/g/data/u86/mxf221/R';

my $sample_file = $OPT{sample_list};

my $readdir = $OPT{readdir};
$readdir = abs_path($readdir);
$readdir .= '/' if $readdir !~ /\/$/;

if (!-d $readdir) {
        modules::Exception->throw("Dir $readdir doesn't exist");
}

my $outdir = $OPT{outdir};
$outdir = abs_path($outdir);
$outdir .= '/' if $outdir !~ /\/$/;

if (!-d $outdir) {
        modules::Exception->throw("Dir $outdir doesn't exist");
}

my $qsubdir = defined $OPT{qsubdir}?$OPT{qsubdir}:$outdir;
$qsubdir = abs_path($qsubdir);

my $read_length = defined $OPT{read_length}?$OPT{read_length}:150;
my $sj = $read_length - 1;


my $gtf = $refdir .'/'.$ref .'/'.$ref .'.gtf';

if ( !-e $gtf ) {
	#modules::Exception->throw("File $gtf doesn't exist");	
}


my @samples = ();
open(SAMPLE, $sample_file) || modules::Exception->throw("File $sample_file doesn't exist");

while (<SAMPLE>) {
	chomp;
	my ($sample) = $_ =~ /(\S+)/;
	push @samples, $sample;
}


#print Dumper \@samples;
my @sample_groups = ();

for my $sample (@samples) {
	open(QSUB, ">$qsubdir/$sample.qsub") || modules::Exception->throw("Can't open qsub for $sample");
	open(TEMPLATE,$template) || modules::Exception->throw("Can't open file $template");
	while (<TEMPLATE>) {
		$_ =~ s/SAMPLE/$sample/g;
		$_ =~ s/OUTDIR/$outdir/g;
		$_ =~ s/READDIR/$readdir/g;
		$_ =~ s/READLENGTH/$read_length/g;
		$_ =~ s/SJDB_value/$sj/g;
		$_ =~ s/REFDIR/$refdir/g;
		$_ =~ s/REF/$ref/g;
		print QSUB $_;
	}	
	push @sample_groups, '"Group1"';
	close QSUB;
	system("qsub $qsubdir/$sample.qsub") if $OPT{run};
	
}



open(ENV,">$qsubdir/env.txt");
print ENV "export R_LIBS_USER=\"$rlib\"\n";


my $sample_string = join(",",@sample_groups);

open(DE_R,">$qsubdir/DE.r");
print DE_R <<EOF;
library(consensusDE)
library(GenomicFeatures)
library(org.Hs.eg.db)
EOF

print DE_R "setwd(\"$outdir\")\n";
print DE_R 'file_list <- list.files(path=getwd(),pattern = ".bam$",full=TRUE)'."\n";
print DE_R "sample_table<-data.frame(\"file\"=basename(file_list),\"group\"=c(".$sample_string ."))\n";
print DE_R "bam_dir <- as.character(gsub(basename(file_list)[1], \"\", file_list[1]))\n";
print DE_R "summarized_dm3 <- buildSummarized(sample_table = sample_table, bam_dir = bam_dir,read_format = \"paired\",output_log=\"".$outdir."/\",gtf=\"".$gtf."\")\n";
close DE_R;




open(DE,">$qsubdir/ConsensusDE.qsub") || modules::Exception->throw("Can't open file $qsubdir/consensusDE.qsub\n");
print DE <<EOF;
#!/bin/bash
#PBS -P pq84
#PBS -q normalbw
#PBS -l walltime=24:00:00,mem=128GB,jobfs=100GB,ncpus=16
#PBS -l storage=gdata/pq84+gdata/u86+scratch/u86+gdata/xx92

module load R/3.6.1
EOF

print DE "source $qsubdir/env.txt\n";
print DE "R --vanilla < $qsubdir/DE.r > $qsubdir/DE.out\n";
close DE;








