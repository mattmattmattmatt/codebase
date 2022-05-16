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
	   	"bamfiles=s",
	   	"bamdir=s",
	   	"output=s",
	   	"source_name=s",
     	"runid=i",
     	"byread",
     	"mdss"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || (!$OPT{bamfiles} && !$OPT{bamdir}) || !$OPT{output} || !$OPT{source_name});

	   
=pod

=head1 SYNOPSIS

merge_bam.pl -bamfiles individual_lane_bam_file -bamdir directory_with_bam_files -source_name source_name(needed_for_destdir) -byread sort_by_read_name -output merge_out_bam  -mdss mdss_copy(default=off) [options]

Required flags: (-bamfiles || -bamdir) -source_name -output

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

merge_bam.pl -> Merge individual bam files together; if single bam file then copies

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

copy_bam.pl 

=cut

my $source_name = $OPT{source_name};
my $outbam = $OPT{output};
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $source_type = modules::Pipeline::get_source_type(-source_name=>$source_name);
my $samtools_bin = $pipe_config->read($source_type,'binaries','samtools','binary');
my $java_bin = $pipe_config->read($source_type,'binaries','javabin','binary');
my $picard_jar = $pipe_config->read($source_type,'binaries','picard','jar');
my $runid = $OPT{runid};
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $byread = defined $OPT{byread}?1:0;

my @bam_files;
if ($OPT{bamfiles}) {
 	@bam_files = split(",",$OPT{bamfiles});
} else {
	my $bamdir = $OPT{bamdir};
	if ( !-d $bamdir ) {
		modules::Exception->throw("File $bamdir doesn't exist");	
	}
	opendir(DIR,$bamdir) || modules::Exception->throw("Can't open bam directory");
	my @bams = grep {/.bam$/} readdir DIR;
	for my $bam ( sort {my ($a_chr) = $a =~ /\.([0-9XY]+)\./; my ($b_chr) = $b =~ /\.([0-9XY]+)\./; $a_chr <=> $b_chr} @bams ) {
	    next unless $bam =~ /recal/;
	    push @bam_files, $bamdir . '/' . $bam;
	}
	#Put X and Y last
	push @bam_files, $bam_files[1], $bam_files[0];
	shift @bam_files for 1..2;
}

if (!@bam_files) {
	modules::Exception->throw("ERROR: No bam files found");
}


for my $bam_file (@bam_files) {
	if (!-e $bam_file) {
		modules::Exception->throw("ERROR: File $bam_file doesn't exist");
	}
}

my @commands = ();

#Here it's a single lane so we don't need to merge; we just copy
if (@bam_files == 1) {
	my $bam_str = $bam_files[0];
	push @commands, "cp $bam_str $outbam";
	push @commands, "cp $bam_str.bai $outbam.bai";
} else {
	#my $bam_str = join(" ",@bam_files);
	#push @commands, "$samtools_bin merge $outbam $bam_str";
	my $bam_str = 'I='.join(" I=",@bam_files);
	if ($byread) {
		push @commands, "$java_bin -jar $picard_jar MergeSamFiles SO=queryname VALIDATION_STRINGENCY=LENIENT TMP_DIR=$run_dir O=$outbam $bam_str";
	} else {
		push @commands, "$java_bin -Djava.io.tmpdir=$run_dir -jar $picard_jar MergeSamFiles VALIDATION_STRINGENCY=LENIENT TMP_DIR=$run_dir O=$outbam $bam_str";
	}
	push @commands, "$samtools_bin index $outbam";
}

my $sys_call = modules::SystemCall->new();

for my $command (@commands) {
	print STDERR "$command\n";
	$sys_call->run($command);	
}



my $mdss_pipe = $cluster_config->read('common','mdss','mdss_flag');
if ($mdss_pipe && $OPT{mdss}) {
	my $mdss_results_dir = $cluster_config->read('common','mdss','mdss_results');
	my $mdss_mkdir_cmd = "mdss mkdir $mdss_results_dir/$source_type/$source_name";
	$sys_call->run($mdss_mkdir_cmd);
	my $mdss_command = "mdss put $outbam ${outbam}.bai $mdss_results_dir/$source_type/$source_name";
	$sys_call->run($mdss_command);  
}








