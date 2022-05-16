#! /usr/bin/perl -w

use strict;
use modules::Adaptors::Variant;
use modules::Adaptors::Filter;
use modules::Overlap;
use modules::Exception;
use modules::Adaptors::BulkInsert;
use modules::ConfigXML;
use modules::SystemCall;
use modules::Vcf;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "overlap_outfile_snv=s",
		   "overlap_outfile_indel=s",
		   "vcf_file=s",
		   "snv_infile=s",
		   "indel_infile=s",
		   "sample_name=s",
		   "phred64",
		   "gauss"	
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf_file}  || !$OPT{overlap_outfile_snv} || !$OPT{overlap_outfile_indel} || !$OPT{snv_infile} || !$OPT{indel_infile} || !$OPT{sample_name});

	   
=pod

=head1 SYNOPSIS

filter_vcf.pl -vcf_file merge_vcf_file -phred64 set_phred64_gatk_arguments -gauss run_with_max_gaussian_4 -overlap_outfile_snv output_snv_overlap_file -overlap_outfile_indel output_indel_overlap_file -indel_infile -snv_infile  [options]

Required flags: -vcf_file -overlap_outfile_snv -overlap_outfile_indel -snv_infile -indel_infile merge_vcf_snv_file -indel_infile merge_vcf_indel_file -sample_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_vcf.pl -> Script to parse and filter a vcf file

=head1 DESCRIPTION

Feb 7, 2012

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_vcf.pl

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $sys_call = modules::SystemCall->new();
my $phred64 = defined $OPT{phred64}?1:0;
#Get the flags for the particular library
my $outfile_indel = $OPT{overlap_outfile_indel};
my $outfile_snv = $OPT{overlap_outfile_snv};
my $sample_name = $OPT{sample_name};
my $source_name = modules::Pipeline::get_source_name($sample_name);
my $sg_name = modules::Pipeline::get_source_group_name($sample_name);
my $source_type = modules::Pipeline::get_source_type(-source_name=>$source_name);
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $base_rundir = $cluster_config->read($source_type,'base_directories','base_run_directory');
my $vcf_file = $OPT{vcf_file};
my $full_vcf = $vcf_file .'.all.vcf';
my $rundir = $cluster_config->read($source_type,'base_directories','base_run_directory');
my $tumour_flag = modules::Pipeline::get_tumour_flag(-sample_name=>$sample_name);

 # Some GATK-specific config
my $java_bin = $pipe_config->read($source_type,'binaries','javabin','binary');
my $gatkjar = $pipe_config->read($source_type,'binaries','gatk','jar');
my $gatksnv = $pipe_config->read($source_type,'cutoffs','gatk_snv_vqsr');
my $gatkindel = $pipe_config->read($source_type,'cutoffs','gatk_indel_vqsr');
my $genome_fasta = $cluster_config->read($source_type,'svn','fasta_file');
my $filenamestub = $rundir.'/'.$source_name.'/'.$sg_name;
my $gatkdatafilesdir = $cluster_config->read($source_type,'svn','gatk_resources');
my $intervals_dir =  $cluster_config->read($source_type,'svn','chr_intervals_dir');

my $vqsr_snv_out = $vcf_file.'.vqsr.vcf.snv';
my $vqsr_indel_out = $vcf_file.'.vqsr.vcf.indel'; 

my $snv_recal_args = $pipe_config->read($source_type,'binaries','gatk','args_snv_var_recal');
my $snv_apply_recal = $pipe_config->read($source_type,'binaries','gatk','args_snv_apply_recal');
my $snv_recal_fallback_args = $pipe_config->read($source_type,'binaries','gatk','args_snv_var_recal_fallback'); #marcin
my $indel_recal_args = $pipe_config->read($source_type,'binaries','gatk','args_indel_var_recal');
my $indel_recal_fallback_args = $pipe_config->read($source_type,'binaries','gatk','args_indel_var_recal_fallback'); #marcin


my $indel_apply_recal = $pipe_config->read($source_type,'binaries','gatk','args_indel_apply_recal');
my $hapmap = "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $gatkdatafilesdir/". $pipe_config->read($source_type,'binaries','gatk','hapmap_file');
my $omni = "-resource:omni,known=false,training=true,truth=true,prior=12.0 $gatkdatafilesdir/". $pipe_config->read($source_type,'binaries','gatk','omni_file');  
my $genome_snv = "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $gatkdatafilesdir/".$pipe_config->read($source_type,'binaries','gatk','genome_snv_file');
my $dbsnp = "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $gatkdatafilesdir/".$pipe_config->read($source_type,'binaries','gatk','dbsnp_file');
my $mills = "-resource:mills,known=false,training=true,truth=true,prior=12.0 $gatkdatafilesdir/".$pipe_config->read($source_type,'binaries','gatk','mills_indel_file');
my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_name);


my $targeted = 0;

#Skip this step with targeted as variant sets aren't big enough
if ($sequence_type eq 'targeted') {
	$sys_call->run("cat $full_vcf | grep -v ^# | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\"PASS\",\$8,\$9,\$10}\' OFS='\t' > $vqsr_snv_out");
	$sys_call->run("cp $vqsr_snv_out $vqsr_indel_out");
	$targeted = 1;
}



#marcin:
#my $snv_recal = "$java_bin -jar $gatkjar $snv_recal_args -R $genome_fasta -input $full_vcf -recalFile $vcf_file.snv.recal -tranchesFile $vcf_file.snv.tranches $hapmap $omni $genome_snv $dbsnp";
my $snv_recal          = "-R $genome_fasta -input $full_vcf -recalFile $vcf_file.snv.recal -tranchesFile $vcf_file.snv.tranches $hapmap $omni $genome_snv $dbsnp";
my $snv_recal_fallback = "$java_bin -jar $gatkjar $snv_recal_fallback_args $snv_recal";
$snv_recal             = "$java_bin -jar $gatkjar $snv_recal_args $snv_recal";

my $snv_apply = "$java_bin -jar $gatkjar $snv_apply_recal -R $genome_fasta -input $full_vcf -recalFile $vcf_file.snv.recal -tranchesFile $vcf_file.snv.tranches -o $vqsr_snv_out -ts_filter_level $gatksnv";

#marcin:
#my $indel_recal = "$java_bin -jar $gatkjar $indel_recal_args -R $genome_fasta -input $full_vcf -recalFile $vcf_file.indel.recal -tranchesFile $vcf_file.indel.tranches $mills $dbsnp";
my $indel_recal = "-R $genome_fasta -input $full_vcf -recalFile $vcf_file.indel.recal -tranchesFile $vcf_file.indel.tranches $mills $dbsnp";
my $indel_recal_fallback = "$java_bin -jar $gatkjar $indel_recal_fallback_args $indel_recal";
$indel_recal = "$java_bin -jar $gatkjar $indel_recal_args $indel_recal";

my $indel_apply = "$java_bin -jar $gatkjar $indel_apply_recal -R $genome_fasta -input $full_vcf -recalFile $vcf_file.indel.recal -tranchesFile $vcf_file.indel.tranches -o $vqsr_indel_out -ts_filter_level $gatkindel";

if ($phred64) {
	$snv_recal .= " --fix_misencoded_quality_scores";
	$snv_recal_fallback .= " --fix_misencoded_quality_scores";  #marcin
	$snv_apply .= " --fix_misencoded_quality_scores";
	$indel_recal .= " --fix_misencoded_quality_scores";
	$indel_recal_fallback .= " --fix_misencoded_quality_scores"; #marcin
  $indel_apply .= " --fix_misencoded_quality_scores";
}

#marcin commented out, it seems to be dead and non-compatible code:
#if ($OPT{gauss}) {
#	$snv_recal .= " --maxGaussians 4";
#}

#Now run the vqsr commands

#marcin:
#$sys_call->run($snv_recal);

if (!$targeted) {
	my $ret = $sys_call->run($snv_recal, 1); #allow for fallback
	if($ret != 1) #1 means all good, anything but 1 or 0 mean command failed; 0 means no command at all
	{
		print STDERR "Running fall back command:\n";
		$sys_call->run($snv_recal_fallback); #no fallback this time, it was our only second chance
	}
	
	$sys_call->run($snv_apply);
	
	#marcin
	#$sys_call->run($indel_recal);
	$ret = $sys_call->run($indel_recal, 1);
	if($ret != 1)
	{
		print STDERR "Running fall back command:\n";
		$sys_call->run($indel_recal_fallback);
	}
	
	$sys_call->run($indel_apply);
}

my $vcf_obj = modules::Vcf->new(-sample_name=>$sample_name);


#Parse the filtered vcf
$vcf_obj->parse_vcf(-vcf_file => $vqsr_indel_out,-tumour_flag=>$tumour_flag);
my $filter_indel_data = $vcf_obj->filter_vcf(-vcf_file=>$vqsr_indel_out, -snv_depth_file=>$OPT{snv_infile},-indel_depth_file=>$OPT{indel_infile},-tumour_flag=>$tumour_flag,-gatk=>1);

if (keys %{$filter_indel_data} == 0) {
	modules::Exception->throw("ERROR: No variants returned filter_vcf call");
}

open(VAR,">$outfile_indel") || modules::Exception->throw("Can't open file $outfile_indel\n");

#Parse the filtered deletions; these will be added to the database
if (exists $filter_indel_data->{DEL}) {
	for my $chr (sort keys %{$filter_indel_data->{DEL}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$filter_indel_data->{DEL}{$chr}}) {
			for my $end_coord (keys %{$filter_indel_data->{DEL}{$chr}{$start_coord}}) {
				for my $bases (keys %{$filter_indel_data->{DEL}{$chr}{$start_coord}{$end_coord}}) {
					if ($bases =~ /N/i) {
						next;
					}
					my $rest = $filter_indel_data->{DEL}{$chr}{$start_coord}{$end_coord}{$bases};
									
		  			print VAR join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							"DEL^^^$bases^^^$rest"
		  							) . "\n";
				}
			}
		}
	}
}

#Parse the insertions
if (exists $filter_indel_data->{INS}) {
	for my $chr (sort sort keys %{$filter_indel_data->{INS}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$filter_indel_data->{INS}{$chr}}) {
			for my $end_coord (keys %{$filter_indel_data->{INS}{$chr}{$start_coord}}) {
				for my $bases (keys %{$filter_indel_data->{INS}{$chr}{$start_coord}{$end_coord}}) {
					if ($bases =~ /N/i) {
						next;
					}
					my $rest = $filter_indel_data->{INS}{$chr}{$start_coord}{$end_coord}{$bases};
					
		  			 print VAR join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							"INS^^^$bases^^^$rest"
		  							) . "\n";
		  			 
				}
			}
		}
	}
}

$vcf_obj->parse_vcf(-vcf_file => $vqsr_snv_out,-tumour_flag=>$tumour_flag);
my $filter_snv_data = $vcf_obj->filter_vcf(-vcf_file=>$vqsr_snv_out, -snv_depth_file=>$OPT{snv_infile},-indel_depth_file=>$OPT{indel_infile},-tumour_flag=>$tumour_flag,-gatk=>1);

if (keys %{$filter_snv_data} == 0) {
        modules::Exception->throw("ERROR: No variants returned filter_vcf call");
}


#Parse the vcf and generate the pileup coords file
open(SNV,">$outfile_snv") || modules::Exception->throw("Can't open file $outfile_snv\n");

#Parse the filtered snvs; these will be added to the database
if (exists $filter_snv_data->{SNV}) {
	for my $chr (sort keys %{$filter_snv_data->{SNV}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$filter_snv_data->{SNV}{$chr}}) {
			for my $end_coord (keys %{$filter_snv_data->{SNV}{$chr}{$start_coord}}) {
				for my $bases (keys %{$filter_snv_data->{SNV}{$chr}{$start_coord}{$end_coord}}) {
					if ($bases =~ /N/i) {
						next;
					}
					my $rest = $filter_snv_data->{SNV}{$chr}{$start_coord}{$end_coord}{$bases};
					
		  			print SNV join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							"SNV^^^$bases^^^$rest"
		  							) . "\n";
				}
			}
		}
	}
}




