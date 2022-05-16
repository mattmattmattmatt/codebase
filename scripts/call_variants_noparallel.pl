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
		   "sample_name=s",
		   "output=s",
		   "outvcf=s",
		   "ref=s",
		   "bam=s",
		   "no_cleanup_files",
		   "rmbam=s",
		   "paired"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sample_name} || !$OPT{output}  || !$OPT{ref} || !$OPT{bam} || !$OPT{outvcf} || !$OPT{rmbam});

	   
=pod

=head1 SYNOPSIS

call_variants_noparallel.pl -sample_name sample_name -output output_file -ref fasta_ref -bam input_bam -rmbam merge_bam_to_replace_with_symlink -outvcf outvcf_file -no_cleanup_file leave_sam_sai_files(default=delete) [options]

Required flags: -sample_name -output -ref -bam -outvcf -rmbam

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

call_variants_noparallel.pl -> Script to calls variants on a single bam file (not per chromosome)

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./call_variants_noparallel.pl  

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

my $source_name = modules::Pipeline::get_source_name($sample_name);

my $source_group_name = modules::Pipeline::get_source_group_name($sample_name);



my $outfile = $OPT{output};
my $sys_call = modules::SystemCall->new();

my $bam = $OPT{bam};
if ( !-e $bam ) {
	modules::Exception->throw("File $bam doesn't exist");	
}

my $ref = $OPT{ref};
if ( !-e $ref ) {
	modules::Exception->throw("File $ref doesn't exist");	
}
#run the commands to generate the vcf
my $vcf = modules::Vcf->new(-sample_name=>$sample_name);
my $all_vcf_file = $OPT{outvcf};

#Get the xml files and create the pipeline object
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $samtools_bin = $pipe_config->read($source_type,'binaries','samtools','binary');
my $bcftools_bin = $pipe_config->read($source_type,'binaries','bcftools','binary');
my $rundir = $cluster_config->read($source_type,'base_directories','base_run_directory');
my $fasta_ref = $cluster_config->read($source_type,'svn','fasta_file');
my $mpileup_args = $pipe_config->read($source_type,'binaries','samtools','mpileup','args');
my $bcftools_args;


my $paired = defined $OPT{paired}?1:0;
 
 
my $pileup_command;
my $bam_unaffected_full;
my $unaffected_sample_name;

if ($paired) {
	
 	if ($source_type ne 'mouse_cancer') {
 		modules::Exception->throw("ERROR: Can only call paired for human related sources or mouse cancer");
 	}
	
	#Get rid of human_related as calling variants in pairs causes snv scores to change as it uses all bases at that position  
	
	#by default mouse cancer has paired flag but if normal set to unpaired
	if (!modules::Pipeline::get_tumour_flag(-sample_name=>$sample_obj->sample_name)) {
		$paired = 0;
		$bcftools_args = $pipe_config->read($source_type,'binaries','bcftools','view','unpaired_args');
	}
} else {
	$bcftools_args = $pipe_config->read($source_type,'binaries','bcftools','view','unpaired_args');
}

if ($paired) {
 	$bcftools_args = $pipe_config->read($source_type,'binaries','bcftools','view','paired_args');
 	
 	my $bam_unaffected;
 	#Get the control bam
 	if ($source_type eq 'mouse_cancer') {
	 	$unaffected_sample_name = $source_group_name . '_normal1';
 		my $bam_norm = $source_group_name . '_normal1.bam';
		$bam_unaffected = $bam_norm;
 	}
 	
	$bam_unaffected_full = $rundir . '/' . $source_name . '/bam_links/'. $bam_unaffected;
	if (!-e $bam_unaffected_full ) {
		modules::Exception->throw("File $bam_unaffected_full doesn't exist");	
	}
	$pileup_command = "$samtools_bin $mpileup_args $fasta_ref $bam_unaffected_full $bam | $bcftools_bin $bcftools_args - | grep -v NNNNN  > $all_vcf_file";
 } else {
 	$pileup_command = "$samtools_bin $mpileup_args $fasta_ref $bam | $bcftools_bin $bcftools_args - | grep -v NNNNN  > $all_vcf_file";
 }
 
#Mpileup command to create vcf
$sys_call->run($pileup_command);
	
my $vcf_obj = modules::Vcf->new(-sample_name=>$sample_name);

if (!$vcf_obj->check_vcf(-vcf=>$all_vcf_file,-chromosome=>'all')) {
	modules::Exception->throw("ERROR: $all_vcf_file fails\n");
}	



(my $snv_vcf_file = $all_vcf_file) =~ s/all.vcf/snv.vcf/;
(my $indel_vcf_file = $all_vcf_file) =~ s/all.vcf/indel.vcf/;

my $snv_out = $outfile.'.snv';
my $indel_out = $outfile .'.indel';

#Create vcfs for the different variant types
my $snv_command = "grep -v INDEL $all_vcf_file > $snv_vcf_file";
print STDERR "$snv_command\n";
$sys_call->run($snv_command);


my $out = `grep INDEL $all_vcf_file 2>/dev/null | head`;

if ($out =~ /INDEL/) {
	my $indel_command = "grep INDEL $all_vcf_file > $indel_vcf_file";
	print STDERR "$indel_command\n";
	$sys_call->run($indel_command);
}

#Parse the unfiltered vcf
$vcf->parse_vcf(-vcf_file => $all_vcf_file);
my $vcf_data = $vcf->get_vcf(-vcf_file => $all_vcf_file);

my $pileup_indel_file = $indel_out . '.pileupcoord';
my $pileup_snv_file = $snv_out . '.pileupcoord';

#Open the coord file for writing; needed for hom/het calls later
open(SNVCOORD,">$pileup_snv_file") || modules::Exception->throw("Can't open file to write $pileup_snv_file\n");
if (exists $vcf_data->{SNV}) {

	for my $chr (sort keys %{$vcf_data->{SNV}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{SNV}{$chr}}) {
			print SNVCOORD "$chr\t$start_coord\n";
		}
	}
}

close SNVCOORD;


#Open the coord file for writing; needed for hom/het calls later
open(INDELCOORD,">$pileup_indel_file") || modules::Exception->throw("Can't open file to write $pileup_indel_file\n");
if (exists $vcf_data->{DEL}) {

	for my $chr (sort keys %{$vcf_data->{DEL}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{DEL}{$chr}}) {
			$start_coord--; #Report from reference base before first deleted base
			print INDELCOORD "$chr\t$start_coord\n";
		}
	}
}

if (exists $vcf_data->{INS}) {

	for my $chr (sort keys %{$vcf_data->{INS}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{INS}{$chr}}) {
			#$start_coord--;
			print INDELCOORD "$chr\t$start_coord\n";
		}
	}
}

close INDELCOORD;


#Now generate the mpileup files
my $mpileup_snv_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_snv_file $bam  > $snv_out.pileup";
print STDERR "$mpileup_snv_command\n";
$sys_call->run($mpileup_snv_command);
	
my $mpileup_indel_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_indel_file $bam  > $indel_out.pileup";
print STDERR "$mpileup_indel_command\n";
$sys_call->run($mpileup_indel_command);

#Need normal mpileups for each tumour or affected
if ($paired) {
	my $snv_out_control = my $indel_out_control;
	if ($source_type eq 'mouse_cancer') {
 		$snv_out_control = $snv_out.'.normal.pileup';
		$indel_out_control = $indel_out.'.normal.pileup';
 	}
	my $mpileup_normal_snv_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_snv_file $bam_unaffected_full  > $snv_out_control";
	print STDERR "$mpileup_normal_snv_command\n";
	$sys_call->run($mpileup_normal_snv_command);
	
	my $mpileup_normal_indel_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_indel_file $bam_unaffected_full  > $indel_out_control";
	print STDERR "$mpileup_normal_indel_command\n";
	$sys_call->run($mpileup_normal_indel_command);
} 

#Get the depths from the pileup; this is needed for filtering in the next step
my %snv_depths = ();
my %indel_depths = ();

open(SNVPILEUP,"$snv_out.pileup") || modules::Exception->throw("ERROR: Can't open file $snv_out.pileup");
while (<SNVPILEUP>) {
	my @fields = split("\t");
	$snv_depths{$fields[0]}{$fields[1]} = $fields[3];
}

open(INDELPILEUP,"$indel_out.pileup") || modules::Exception->throw("ERROR: Can't open file $indel_out.pileup");
while (<INDELPILEUP>) {
	my @fields = split("\t");
	my $indel_start = $fields[1] + 1;
	$indel_depths{$fields[0]}{$indel_start} = $fields[3];
}

open(INDEL,">$indel_out") || modules::Exception->throw("Can't open file to write $indel_out\n");
open(SNV,">$snv_out") || modules::Exception->throw("Can't open file to write $snv_out\n");


#Write out unfiltered variants to text files
if (exists $vcf_data->{DEL}) {
	for my $chr (sort keys %{$vcf_data->{DEL}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{DEL}{$chr}}) {
			for my $end_coord (keys %{$vcf_data->{DEL}{$chr}{$start_coord}}) {
				for my $bases (keys %{$vcf_data->{DEL}{$chr}{$start_coord}{$end_coord}}) {
					if ($bases =~ /N/i) {
						next;
					}
					
					my $depth = 0;
					#If not defined then tumour has no reads
					if (defined $indel_depths{$chr}{$start_coord}) {
						$depth = $indel_depths{$chr}{$start_coord};
					} 
					
					my $rest = 'DEL^^^' . $bases . '^^^' .$vcf_data->{DEL}{$chr}{$start_coord}{$end_coord}{$bases}.'^^^D'.$depth;
					print INDEL join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							$rest
		  							) ."\n";
				}
			}
		}
	}
}

if (exists $vcf_data->{INS}) {
	for my $chr (sort keys %{$vcf_data->{INS}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{INS}{$chr}}) {
			for my $end_coord (keys %{$vcf_data->{INS}{$chr}{$start_coord}}) {
				for my $bases (keys %{$vcf_data->{INS}{$chr}{$start_coord}{$end_coord}}) {
					if ($bases =~ /N/i) {
						next;
					}
					my $depth = 0;
					my $new_start = $start_coord + 1; #Account for base we add to deletion calls
					if (defined $indel_depths{$chr}{$new_start}) {
						$depth = $indel_depths{$chr}{$new_start};
					} 
					my $rest = 'INS^^^' . $bases . '^^^' .$vcf_data->{INS}{$chr}{$start_coord}{$end_coord}{$bases}.'^^^D'.$depth;
					print INDEL join("\t",
		  							$chr,
		  							$start_coord,
		  							$end_coord,
		  							$rest
		  							) ."\n";
				}
			}
		}
	}
}

close INDEL;

if (exists $vcf_data->{SNV}) {
	for my $chr (sort keys %{$vcf_data->{SNV}}) {
		for my $start_coord (sort {$a<=>$b} keys %{$vcf_data->{SNV}{$chr}}) {
			for my $end_coord (keys %{$vcf_data->{SNV}{$chr}{$start_coord}}) {
				for my $base_change (keys %{$vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}}) {
					if ($base_change =~ /N/i) {
						next;
					}
					my $depth = 0;
					if (defined $snv_depths{$chr}{$start_coord}) {
						$depth = $snv_depths{$chr}{$start_coord};
					} 
					my $rest = 'SNV^^^'.$base_change .'^^^' . $vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}{$base_change} . '^^^D'.$depth;
					print SNV join("\t",
									$chr,
									$start_coord,
									$end_coord,
									$rest
									)."\n";
				}
			}
		}
	}
}

close SNV;

#Clean up lanes .sam and .sai files and change merge_bam to symlinks
unless ($OPT{no_cleanup_files}) {
	my $rmbam = $OPT{rmbam};
	my $rmbam_index = $rmbam .'.bai';
	if ( !-e $rmbam ) {
		modules::Exception->throw("File $rmbam doesn't exist");	
	}
	
	my $rmbam_just_file = basename($rmbam);
	
	my $destdir = $cluster_config->read($source_type,'base_directories','base_results_directory') . '/'. $source_name;
	
	my $bam_mass = $destdir . '/' . $rmbam_just_file;
	my $bam_mass_index = $bam_mass . '.bai';
	
	if ( !-e $bam_mass ) {
		modules::Exception->throw("File $bam_mass doesn't exist");	
	}
	
	$sys_call->run("rm $rmbam; rm $rmbam_index");
	$sys_call->run("ln -s $bam_mass $rmbam; ln -s $bam_mass_index $rmbam_index");
	
	
	my $lane_dir = $rundir . '/' . $source_name . '/' . $sample_name . '_lanes';
	opendir(LANES,$lane_dir) || modules::Exception->throw("ERROR: Can't open lanes dir $lane_dir");
	my @delete_files = grep {/[bs]a[im]$/} readdir LANES;
	closedir LANES;
	
	for my $file ( @delete_files ) {
	    my $full_delete_file = $lane_dir . '/' . $file;
	    if ( !-e $full_delete_file ) {
	    	modules::Exception->throw("File $full_delete_file doesn't exist");	
	    }
	    $sys_call->run("rm $full_delete_file");
	}
	
	#Remove the merge_bam files and replace with symlinks
	if ($paired) {
		#Clean up the normal bam files as well
		#Get the normal bam merge_bam step file by reading the symlink and changing the rmdup file to merge_bam file
		my $rmbam_normal = readlink($bam_unaffected_full);
		$rmbam_normal =~ s/remove_duplicates/merge_bam/;
		my $rmbam_normal_index = $rmbam_normal .'.bai';
		if ( !-e $rmbam_normal ) {
			modules::Exception->throw("File $rmbam_normal doesn't exist");	
		}
		
	
		my $rmbam_normal_just_file = basename($rmbam_normal);
		
		my $bam_normal_mass = $destdir . '/' . $rmbam_normal_just_file;
		my $bam_normal_mass_index = $bam_normal_mass . '.bai';
		if ( !-e $bam_normal_mass ) {
			modules::Exception->throw("File $bam_normal_mass doesn't exist");	
		}
		
		$sys_call->run("rm $rmbam_normal; rm $rmbam_normal_index");
		$sys_call->run("ln -s $bam_normal_mass $rmbam_normal; ln -s $bam_normal_mass_index $rmbam_normal_index");
		
		#Clean up the normal lanes as well
		my $normal_lane_dir = $rundir . '/' . $source_name . '/' . $unaffected_sample_name . '_lanes';
		opendir(NLANES,$normal_lane_dir) || modules::Exception->throw("ERROR: Can't open lanes dir $normal_lane_dir");
		my @normal_delete_files = grep {/sa[im]$/} readdir NLANES;
		closedir NLANES;
	
		for my $file ( @normal_delete_files ) {
	    	my $full_delete_file = $normal_lane_dir . '/' . $file;
	    	if ( !-e $full_delete_file ) {
	    		modules::Exception->throw("File $full_delete_file doesn't exist");	
	    	}
	    	$sys_call->run("rm $full_delete_file");
		}	
		
	} 
}


#Don't need this file
$sys_call->run("rm $pileup_indel_file");
$sys_call->run("rm $pileup_snv_file");


