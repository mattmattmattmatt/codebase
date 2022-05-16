#! /usr/bin/perl -w

use strict;
use modules::Overlap;
use modules::Exception;
use modules::ConfigXML;
use modules::SystemCall;
use modules::Pipeline;
use modules::Vcf;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "sample_name=s"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sample_name});

	   
=pod

=head1 SYNOPSIS

get_geno_snp_validation.pl -sample_name sample_name [options]

Required flags: -sample_name 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

get_geno_snp_validation -> Script to call genotypes for follow-up snp validation

=head1 DESCRIPTION

Sep 10, 2014

a script that ...

=head1 AUTHOR

Vicky Cho

=head1 EXAMPLE

./get_geno_snp_validation.pl -sample_name APOSLE_single87_sg1_humansingle1

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $sample_name = $OPT{sample_name};
my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);
my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);
if (!$sample_obj) {
	modules::Exception->throw("ERROR: Cannot retrieve sample object for $sample_name");
}
my $sample_type = $sample_obj->sample_type;
my $sample_id = $sample_obj->id;
my $external_sample_name = $sample_obj->external_sample_name;

my $source_name = modules::Pipeline::get_source_name($sample_name);

my $source_group_name = modules::Pipeline::get_source_group_name($sample_name);

my $clus_conf = modules::Pipeline::get_cluster_conf();
my $resultdir = $clus_conf->read($source_type,'base_directories','base_results_directory');
my $outdir = "/g/data/u86/snv_pipeline_runs/v2.1_human_genome/human_snp_validation";
my $outfile = "$outdir/$sample_name.genotype.tsv";
my $sys_call = modules::SystemCall->new();


#Get the runid from the database
my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_id);
if (!defined $run_obj) {
	modules::Exception->throw("Cannot retrieve any runs for $sample_name");
}
my $runid = $run_obj->id;

my $bam = "$resultdir/$source_name/$sample_name\_$runid.merge_bam.out.bam";
if ( !-e $bam ) {
	modules::Exception->throw("File $bam doesn't exist");	
}

my $ref = "$svndir/conf/human/GRCh37/fasta/single_file/GRCh37.fa";
if ( !-e $ref ) {
	modules::Exception->throw("File $ref doesn't exist");	
}

my $vcf_file = "$outdir/$sample_name.vcf";

#Get the xml files and create the pipeline object
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $samtools_bin = $pipe_config->read($source_type,'binaries','samtools','binary');
my $bcftools_bin = $pipe_config->read($source_type,'binaries','bcftools','binary');
my $rundir = $cluster_config->read($source_type,'base_directories','base_run_directory');
my $fasta_ref = $cluster_config->read($source_type,'svn','fasta_file');
my $mpileup_args = $pipe_config->read($source_type,'binaries','samtools','mpileup','args');	#"-uAEf"
my $bcftools_args = "view -g";

my $snp_list_file = "$svndir/conf/human_snp_validation_list.tsv";
 
my $pileup_command = "$samtools_bin $mpileup_args $fasta_ref -l $snp_list_file $bam | $bcftools_bin $bcftools_args - > $vcf_file";;

 
#Mpileup command to create vcf
$sys_call->run($pileup_command);
	
open(VCF,$vcf_file) || modules::Exception->throw("Can't open file $vcf_file\n");
my %variant_data = ();
    
while (<VCF>) {
    	next if /^#/;
    	
    	if ($_ !~ /^chr[0-9XYM]/) {
    		modules::Exception->throw("ERROR: Error with VCF format at line $_");
    	}
    	
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,$normal_gt) = split("\t");
    	$chr =~ s/^chr//;
	my @gt = split(":",$normal_gt);
    	
    	my $genotype;
    	
    	if($gt[0] eq '0/1'){
    		$genotype = "$ref/$var_str";
    	}elsif($gt[0] eq '1/1'){
    		$genotype = "$var_str/$var_str";
    	}elsif($gt[0] eq '0'){
    		$genotype = "$ref/$ref";
    	}
    	
    	if (!$qual) {
    		modules::Exception->throw("ERROR: Error with VCF format at line $_");
    	}
		$variant_data{$chr}{$first_coord} = $genotype;		
}

#print Dumper \%variant_data;


open(SNPLIST,$snp_list_file) || modules::Exception->throw("Can't open file $snp_list_file\n");
my %snp_list = ();

while (<SNPLIST>) {    	
    	my ($chr,$first_coord,$snpID) = split("\t");
	$chr =~ s/^chr//;
    	$snp_list{$chr}{$first_coord} = $snpID;
}
close SNPLIST;

open(OUTFILE,">$outfile") || modules::Exception->throw("Can't open file to write $outfile\n");

print OUTFILE "chr\tcoord\t$sample_name:$external_sample_name\n";		#header

for my $chr (sort {$a<=>$b} keys %snp_list){
	for my $coord (sort {$a<=>$b} keys %{$snp_list{$chr}}){
		if($variant_data{$chr}{$coord}){
			print OUTFILE "$chr\t$coord\t$variant_data{$chr}{$coord}\n";
		}else{
			print OUTFILE "$chr\t$coord\tNA\n";
		}
	}
	
}

close OUTFILE;

#Don't need this file
$sys_call->run("rm $vcf_file");


