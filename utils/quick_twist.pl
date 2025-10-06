#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
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
		"pref=s",
		"template=s",
		"qsubdir=s",
		"softdir=s",
		"run"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{readdir} || !$OPT{outdir} || !$OPT{sample_list} || !$OPT{template});

	   
=pod

=head1 SYNOPSIS

quick_twist.pl -template template_qsub -qsubdir qsub_dir -readdir readdir(fastq_must_be_named_'sample_R[12]_fastq.gz') -outdir outdir -sample_list sample_list_file -prefix file_prefix -softdir software_directory

Required flags: -readdir -outdir -sample_list -template

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

my $sample_file = $OPT{sample_list};

my $readdir = $OPT{readdir};
$readdir = abs_path($readdir);

if (!-d $readdir) {
        modules::Exception->throw("Dir $readdir doesn't exist");
}

my $pref = defined $OPT{pref}?$OPT{pref}:'twist';


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


my $softdir = defined $OPT{softdir}?$OPT{softdir}:'/g/data/pq84/software/';
$softdir = abs_path($softdir);

if (!-d $softdir) {
        #modules::Exception->throw("Dir $softdir doesn't exist");
}



#print Dumper \@samples;


for my $sample (@samples) {
		my $qsub = "${sample}_${pref}.qsub";
		open(QSUB, ">$qsubdir/$qsub") || modules::Exception->throw("Can't open qsub for $sample");
		open(TEMPLATE,$template) || modules::Exception->throw("Can't open file $template");
		
		while (<TEMPLATE>) {
			$_ =~ s/SAMPLE/$sample/g;
			$_ =~ s/OUTDIR/$outdir/g;
			$_ =~ s/READDIR/$readdir/g;
			$_ =~ s/SOFT/$softdir/g;
			
			print QSUB $_;
			
		}	
		close QSUB;
	
}
	







