#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;
use modules::CLI_Options;
use modules::Adaptors::Source;

GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m"
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

create_cluster_xml.pl [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

create_cluster_xml.pl -> generate a config file for all paths, binaries, etc on a cluster; only need to run once

=head1 DESCRIPTION

Nov 30, 2011

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

create_cluster_xml.pl 

=cut


my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
	print "SVN directory is $svndir...\n";
}

#Get the source types for iteration
my @source_objs = modules::Adapators::Source->search_all();

my %source_types = ();

for my $source_obj (@source_objs) {
	$source_types{$source_obj->source_type}++;
}



my $options = modules::CLI_Options->new();

my $xml_cluster = "$svndir/conf/cluster.xml";
my $xml_pipe_template = "$svndir/conf/pipe.xml";

open(CLUSTER,">$xml_cluster") || modules::Exception->throw("ERROR: Cannot open xml $xml_cluster file");
print CLUSTER "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n<!DOCTYPE page [\n";

#Output from get_option command
my $error_flag = my $error_msg = my $input_variable;

#Get the pipeline version
my @pipeline_options = qw(v1.0 v1.1 v2.0 local);
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'pipeline_version', -default => 'v2.0', -not_empty => 1, -possibilities => \@pipeline_options);
my $pipeline_version = $input_variable;

if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'pipeline_version', -default => 'v2.0', -not_empty => 1, -possibilities => \@pipeline_options, -retry => 1);
}
print CLUSTER "<!ENTITY pipeline_version \"$pipeline_version\">\n]>\n";

print CLUSTER "<cluster_resources>\n";

print CLUSTER "\t<common>\n";
print CLUSTER "\t\t<pipeline_version>$pipeline_version</pipeline_version>\n\n";


my @schedulers = qw(PBS SGE);
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'scheduler', -default => 'PBS', -not_empty => 1, -possibilities => \@schedulers);
if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'scheduler', -default => 'PBS', -not_empty => 1, -possibilities => \@schedulers, -retry => 1);
}
print CLUSTER "\t\t<scheduler>$input_variable</scheduler>\n";

#Check if we scp or cp the results to an archive directory
my @archive_options = qw(cp scp);
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'archive with scp or cp', -default => 'cp', -yes=>1, -possibilities => \@archive_options);
if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'archive with scp or cp', -default => 'cp', -yes=>1, -possibilities => \@archive_options, -retry => 1);
}
my $archive_command = $input_variable;

print CLUSTER "\t\t<archive>\n";
print CLUSTER "\t\t\t<archive_command>$archive_command</archive_command>\n";

if ($archive_command eq 'scp') {
	print STDERR "NOTE: keyless ssh must be set up for this option to work...\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'remote_server', -default=> '150.203.159.204',-not_empty => 1);
	print CLUSTER "\t\t\t<file_server>$input_variable</file_server>\n";
} else {
	print CLUSTER "\t\t\t<file_server>N/A</file_server>\n"
}
print CLUSTER "\t\t</archive>\n\n";

print CLUSTER "\t\t<svndir>$svndir</svndir>\n";
#Set all the svn paths
print CLUSTER "\t\t<script_dir>$svndir/scripts</script_dir>\n";


if (!-d "$svndir/scripts") {
	modules::Exception->throw("ERROR: Conf directory $svndir/scripts doesn't exist\n");
}

print CLUSTER "\t\t<qsub_vars>\n";
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'max threads available / node', -default => 8, -regex => '^\d+$', -not_empty => 1);
if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'max threads available / node', -default => 8, -regex => '^\d+$', -not_empty => 1, -retry => 1);
}
print CLUSTER "\t\t\t<thread_num>$input_variable</thread_num>\n";

($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'max memory available/ node in GBs', -default => 16, -regex => '^\d+$', -not_empty => 1);
if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'max memory available / node in GBs', -default => 16, -regex => '^\d+$', -not_empty => 1, -retry => 1);
}
my $gb_mem = $input_variable.'GB';
print CLUSTER "\t\t\t<max_memory>$gb_mem</max_memory>\n";
print CLUSTER "\t\t</qsub_vars>\n\n";

print CLUSTER "\t</common>\n";


($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'external binary directory', -default => "$svndir/ext/bin", -not_empty => 1, -directory=>1);
my $ext_bin = $input_variable;
if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'external binary directory', -default => 'v1.0', -not_empty => 1, -possibilities => \@pipeline_options, -retry => 1);
}

open(PIPE,">$xml_pipe_template") || modules::Exception->throw("ERROR: Cannot open xml $xml_pipe_template file");
print PIPE "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n<!DOCTYPE page []>\n";
print PIPE "<pipeline>\n";
print PIPE "\t<common>\n";
print PIPE "\t\t<pipeline_version>$pipeline_version</pipeline_version>\n\n";

#Get db info
my $database_name_default = 'melanoma_production_'.$pipeline_version;
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'database name', -default => $database_name_default);
print PIPE "\t\t<database>\n\t\t<db_name>$input_variable</db_name>\n";
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'database host', -default => 'localhost');
print PIPE "\t\t\t<host>$input_variable</host>\n";
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'database user', -default => 'melanoma_user');
print PIPE "\t\t\t<user>$input_variable</user>\n\t\t</database>\n\n";
print PIPE "\t\t<run_subdirs>overlap,vcf,conf,log,polyphen,summary,bam</run_subdirs>\n\n";
($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'sequencing_centres', -default => "BRF,AGRF,RAM,BGI,Macrogen", -not_empty => 1);
if ($error_flag) {
	print "Error with value entered: $error_msg\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'sequencing_centres', -default => 'BRF,AGRF,RAM,BGI,Macrogen', -not_empty => 1);
}
print PIPE "\t\t<sequencing_centres>$input_variable</sequencing_centres>\n";


print PIPE "\t</common>\n";


#Now the specific source type stuff
for my $source_type (keys %source_types) {
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'reference genome', -default => 'GRCh37', -not_empty => 1);
	my $ref_genome = $input_variable;
	
	print "\nNext questions for source_type $source_type....\n\n";
	print CLUSTER "\t<$source_type>\n";
	print PIPE "\t<$source_type>\n";
	
	print CLUSTER "\t\t<svn>\n";
	
	my $conf_dir = "$svndir/conf/human/$input_variable";
	
	if (!-d $conf_dir) {
		modules::Exception->throw("ERROR: Conf directory $conf_dir doesn't exist\n");
	}
	
	my $fasta_file = "$conf_dir/fasta/single_file/$input_variable.fa";
	if (!-e $fasta_file) {
		modules::Exception->throw("ERROR: No fasta file for reference genome $fasta_file\n");
	}
	
	print CLUSTER "\t\t\t<conf_dir>$conf_dir</conf_dir>\n";
	print CLUSTER "\t\t\t<fasta_file>$fasta_file</fasta_file>\n";
	print CLUSTER "\t\t\t<bwa_index>$conf_dir/bwa_index/$input_variable</bwa_index>\n";
	print CLUSTER "\t\t\t<vep_index>$conf_dir/vep_index/homo_sapiens</vep_index>\n";
	
	
	print CLUSTER "\t\t</svn>\n\n";
	
	
	
	
	
	#Get the base directories
	
	print CLUSTER "\t\t<base_directories>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base read directory', -default => '/massdata/u86/melanoma/reads', -not_empty => 1, -directory=>1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base read directory', -default => '/massdata/u86/melanoma/reads', -not_empty => 1, -directory=>1, -retry => 1);
	}
	print CLUSTER "\t\t\t<base_read_directory>$input_variable</base_read_directory>\n";
	
	if ($archive_command eq 'cp') {
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base results directory', -default => "/massdata/u86/melanoma/results/$pipeline_version", -not_empty => 1, -directory=>1);
		if ($error_flag) {
			print "Error with value entered: $error_msg\n";
			($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base results directory', -default => "/massdata/u86/melanoma/results/$pipeline_version", -not_empty => 1, -directory=>1, -retry => 1);
		}
		#Make sure this directory has pipeline version at the end
		if ($input_variable !~ /$pipeline_version$/) {
			$input_variable .= "/$pipeline_version";
		}
		print CLUSTER "\t\t\t<base_results_directory>$input_variable</base_results_directory>\n";
	} else {
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'remote base results directory', -default => "/raid2/melanoma_results/$pipeline_version", -not_empty => 1);
		if ($error_flag) {
			print "Error with value entered: $error_msg\n";
			($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'remote base results directory', -default => "/raid2/melanoma_results/$pipeline_version", -not_empty => 1, -retry => 1);
		}
		#Make sure this directory has pipeline version at the end
		if ($input_variable !~ /$pipeline_version$/) {
			$input_variable .= "/$pipeline_version";
		}
		print CLUSTER "\t\t\t<base_results_directory>$input_variable</base_results_directory>\n";
	}
	
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base qsub directory', -default => "/g/data/u86/melanoma/qsub/$pipeline_version", -not_empty => 1, -directory=>1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base qsub directory', -default => "/g/data/u86/melanoma/qsub/$pipeline_version", -not_empty => 1, -directory=>1, -retry => 1);
	}
	#Make sure this directory has pipeline version at the end
	if ($input_variable !~ /$pipeline_version$/) {
		$input_variable .= "/$pipeline_version";
	}
	print CLUSTER "\t\t\t<base_qsub_directory>$input_variable</base_qsub_directory>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base runs directory', -default => "/g/data/u86/melanoma/runs/$pipeline_version", -not_empty => 1, -directory=>1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'base runs directory', -default => "/g/data/u86/melanoma/runs/$pipeline_version", -not_empty => 1, -directory=>1, -retry => 1);
	}
	#Make sure this directory has pipeline version at the end
	if ($input_variable !~ /$pipeline_version$/) {
		$input_variable .= "/$pipeline_version";
	}
	print CLUSTER "\t\t\t<base_run_directory>$input_variable</base_run_directory>\n";
	
	print CLUSTER "\t\t</base_directories>\n\n";

	
	#Get the latest gene, ccds, and dbsnp entries
	print PIPE "\t\t<annotation_version>\n";
	print PIPE "\t\t\t<ref_genome>$ref_genome</ref_genome>\n";
	
	#Get the highest number dbsnp
	opendir(DBSNP,"$conf_dir/dbsnp/") || modules::Exception->throw("ERROR: Cannot open dbsnp directory $conf_dir/dbsnp/");
	my @dbsnp_files = grep {/^\d/} readdir DBSNP;
	closedir DBSNP;
	my ($dbsnp_default) = reverse(sort {$a<=>$b} @dbsnp_files);
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'dbsnp version', -default => $dbsnp_default, -not_empty => 1);
	my $dbsnp_dir = $conf_dir.'/dbsnp/'.$input_variable;
	if (!-d $dbsnp_dir) {
		modules::Exception->throw("ERROR: Conf directory $dbsnp_dir doesn't exist\n");		
	}
	print PIPE "\t\t\t<dbsnp_version>$input_variable</dbsnp_version>\n";
	
	#Get the highest number ensembl
	my $vep_index_dir = "$conf_dir/vep_index/human"; #Use this for ensembl version
	opendir(VEP,$vep_index_dir) || modules::Exception->throw("ERROR: Cannot open dbsnp directory $vep_index_dir");
	my @ensembl_files = grep {/^\d/} readdir VEP;
	closedir VEP;
	my ($ensembl_default) = reverse(sort {$a<=>$b} @ensembl_files);
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'ensembl version', -default => $ensembl_default, -not_empty => 1);
	my $vep_dir = $vep_index_dir.'/'.$input_variable;
	if (!-d $vep_dir) {
		modules::Exception->throw("ERROR: Conf directory $vep_dir doesn't exist\n");		
	}
	print PIPE "\t\t\t<ensembl_version>$input_variable</ensembl_version>\n";
	
	#Get the latest exon directory
	opendir(EXON,"$conf_dir/exon/") || modules::Exception->throw("ERROR: Cannot open exon directory $conf_dir/exon/");
	my @exon_files = grep {/^\d/} readdir EXON;
	closedir EXON;
	my ($exon_default) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @exon_files);
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'exon directory', -default => $exon_default, -not_empty => 1);
	my $exon_dir = $conf_dir.'/exon/'.$input_variable;
	if (!-d $exon_dir) {
		modules::Exception->throw("ERROR: Conf directory $exon_dir doesn't exist\n");		
	}
	print PIPE "\t\t\t<exon_version>$input_variable</exon_version>\n";
	
	#Get the latest ccds directory
	opendir(COSMIC,"$conf_dir/cosmic/") || modules::Exception->throw("ERROR: Cannot open ccds directory $conf_dir/cosmic/");
	my @cosmic_files = grep {/^\d/} readdir COSMIC;
	closedir COSMIC;
	my ($cosmic_default) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @cosmic_files);
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'cosmic directory', -default => $cosmic_default, -not_empty => 1);
	my $cosmic_dir = $conf_dir.'/cosmic/'.$input_variable;
	if (!-d $cosmic_dir) {
		modules::Exception->throw("ERROR: Conf directory $cosmic_dir doesn't exist\n");		
	}
	print PIPE "\t\t\t<cosmic_version>$input_variable</cosmic_version>\n";
	
	#Get the latest gene directory
	opendir(GENE,"$conf_dir/gene/") || modules::Exception->throw("ERROR: Cannot open gene directory $conf_dir/gene/");
	my @gene_files = grep {/^\d/} readdir GENE;
	closedir GENE;
	my ($gene_default) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @gene_files);
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'gene directory', -default => $gene_default, -not_empty => 1);
	my $gene_dir = $conf_dir.'/gene/'.$input_variable;
	if (!-d $gene_dir) {
		modules::Exception->throw("ERROR: Conf directory $gene_dir doesn't exist\n");		
	}
	
	print PIPE "\t\t\t<gene_version>$input_variable</gene_version>\n";
	
	print PIPE "\t\t\t<chr>1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y</chr>\n";
	
	print PIPE "\t\t</annotation_version>\n\n";
	
	
	
	#Get the binary info; these are written to pipeline template files
	print PIPE "\t\t<binaries>\n";
	
	if ($source_type =~ /mouse/) {
		print PIPE "\t\t\t<polyphen>\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'polyphen base', -default => "", -not_empty => 1);
		if (!-d $input_variable) {
			modules::Exception->throw("ERROR: Problem with polyphen base $input_variable\n");
		}
		my $polyphen_base = $input_variable;
		
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'polyphen binary', -default => "$polyphen_base/bin/pph", -not_empty => 1);
		if (!-x $input_variable) {
			modules::Exception->throw("ERROR: Problem with polyphen binary $input_variable\n");
		}
		print PIPE "\t\t\t\t<binary>$input_variable</binary>\n";
		
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'polyphen version', -default => "2.1.0", -not_empty => 1);
		print PIPE "\t\t\t\t<version>$input_variable</version>\n";		
		my $polyphen_version  = $input_variable;
		
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'polyphen conf directory', -default => "$polyphen_base/mouse_config");
		print PIPE "\t\t\t\t<polyphen_conf>$input_variable</polyphen_conf>\n";
		print PIPE "\t\t\t</polyphen>\n";
	}
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'bwa binary', -default => "$ext_bin/bwa", -not_empty => 1);
	if (!-x $input_variable) {
		modules::Exception->throw("ERROR: Problem with bwa executable $input_variable\n");
	}
	print PIPE "\t\t\t\t<binary>$input_variable</binary>\n";
	my $bwa_version = &Get_Binary_Version($input_variable);
	print PIPE "\t\t\t\t<version>$bwa_version</version>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'bwa aln arguments (other than -I for phred 64)', -default => '');
	print PIPE "\t\t\t\t<aln>\n\t\t\t\t\t<args>aln $input_variable</args>\n\t\t\t\t</aln>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'bwa sampe arguments', -default => '');
	print PIPE "\t\t\t\t<sampe>\n\t\t\t\t\t<args>sampe $input_variable</args>\n\t\t\t\t</sampe>\n";
	print PIPE "\t\t\t</bwa>\n\n";
	
	print PIPE "\t\t\t<samtools>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'samtools binary', -default => "$ext_bin/samtools", -not_empty => 1);
	if (!-x $input_variable) {
		modules::Exception->throw("ERROR: Problem with samtools executable $input_variable\n");
	}
	print PIPE "\t\t\t\t<binary>$input_variable</binary>\n";
	my $samtools_version = &Get_Binary_Version($input_variable);
	print PIPE "\t\t\t\t<version>$samtools_version</version>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'samtools mpileup arguments', -default => '-C50 -uDEf');
	print PIPE "\t\t\t\t<mpileup>\n\t\t\t\t\t<args>mpileup $input_variable</args>\n\t\t\t\t</mpileup>\n";
	print PIPE "\t\t\t</samtools>\n\n";
	
	print PIPE "\t\t\t<bcftools>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'bcftools binary', -default => "$ext_bin/bcftools", -not_empty => 1);
	if (!-x $input_variable) {
		modules::Exception->throw("ERROR: Problem with bcftools executable $input_variable\n");
	}
	print PIPE "\t\t\t\t<binary>$input_variable</binary>\n";
	my $bcf_version = &Get_Binary_Version($input_variable);
	print PIPE "\t\t\t\t<version>$bcf_version</version>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'bcftools view arguments for unmatched samples', -default => '-vcg');
	print PIPE "\t\t\t\t<view>\n\t\t\t\t\t<unpaired_args>view $input_variable</unpaired_args>\n\t\t\t\t</view>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'bcftools view arguments for matched sample', -default => '-vcg -T pair');
	print PIPE "\t\t\t\t<view>\n\t\t\t\t\t<paired_args>view $input_variable</paired_args>\n\t\t\t\t</view>\n";
	print PIPE "\t\t\t</bcftools>\n\n";
	
	print PIPE "\t\t\t<variant_predictor>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'variant predictor binary', -default => "$ext_bin/variant_effect_predictor.pl", -not_empty => 1);
	if (!-x $input_variable) {
		modules::Exception->throw("ERROR: Problem with variant predictor executable $input_variable\n");
	}
	print PIPE "\t\t\t\t<binary>$input_variable</binary>\n";
	my $variant_predictor_version = &Get_Binary_Version($input_variable);
	print PIPE "\t\t\t\t<version>$variant_predictor_version</version>\n";
	print PIPE "\t\t\t</variant_predictor>\n\n";
	
	print PIPE "\t\t</binaries>\n\n";
	
	#Now get the cutoff
	print PIPE "\t\t<cutoffs>\n";
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'mpileup snv quality cutoff', -default => 40, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'mpileup snv quality cutoff', -default => 40, -regex => '^\d+$', -not_empty => 1, -retry => 1);
	}
	print PIPE "\t\t\t<snv_quality_cutoff>$input_variable</snv_quality_cutoff>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'mpileup clr cutoff', -default => 60, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'mpileup clr cutoff', -default => 60, -regex => '^\d+$', -not_empty =>1,  -retry => 1);
	}
	print PIPE "\t\t\t<clr_cutoff>$input_variable</clr_cutoff>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'mpileup indel quality cutoff', -default => 40, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'mpileup indel quality cutoff', -default => 40, -regex => '^\d+$', -not_empty => 1, -retry => 1);
	}
	print PIPE "\t\t\t<indel_quality_cutoff>$input_variable</indel_quality_cutoff>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'max variant cover', -default => 1000, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'max variant cover', -default => 1000, -regex => '^0.\d+$', -not_empty => 1, -retry => 1);
	}
	print PIPE "\t\t\t<max_variant_cutoff>$input_variable</max_variant_cutoff>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'min variant cover', -default => 2, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'min variant cover', -default => 2, -regex => '^\d+$', -not_empty => 1, -retry => 1);
	}
	print PIPE "\t\t\t<min_variant_cutoff>$input_variable</min_variant_cutoff>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'min read number needed to be called covered', -default => 4, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'min read number needed to be called covered', -default => 4, -regex => '^\d+$', -not_empty => 1, -retry => 1);
	}
	print PIPE "\t\t\t<min_cover_cutoff>$input_variable</min_cover_cutoff>\n";
	
	($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'splice site bp cutoff', -default => 10, -regex => '^\d+$', -not_empty => 1);
	if ($error_flag) {
		print "Error with value entered: $error_msg\n";
		($input_variable,$error_flag,$error_msg) = $options->get_option(-name => 'splice site bp cutoff', -default => 10, -regex => '^\d+$', -not_empty => 1, -retry => 1);
	}
	print PIPE "\t\t\t<splice_length_cutoff>$input_variable</splice_length_cutoff>\n";
	
	print PIPE "\t\t</cutoffs>\n";
	

	print PIPE "\t</$source_type>\n";
	
	print CLUSTER "\t</$source_type>\n";
	
}

print CLUSTER "</cluster_resources>";
close CLUSTER;

print PIPE "</pipeline>\n";
close PIPE;









#Sub routine that tries to extract the binary given the executable
sub Get_Binary_Version {
	my ($full_bin_path) = @_;
	my @lines = split("\n", `$full_bin_path 2>&1`);
	my ($out) = grep {/version/i} @lines;
	my $version;
	#This regex happens to match all out binaries
	if (defined $out && $out =~ /Version:\s+(.+)$/) {
		$version = $1;
	} elsif (defined $out && $out =~ /^version\s+(.+)/) {
		$version = $1;
	} else {
		my ($usage) = grep {/Usage/} @lines;
		#If it doesn't spit out a version then just check the Usage statement exists
		if (defined $usage) {
			$version ="NO_VERSION_GIVEN";
		} else {
    		modules::Exception->throw("ERROR: Cannot get version for binary $full_bin_path")				
		}
	}
	return $version;
}


