#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;
use modules::Pipeline;
use modules::Adaptors::Sample;
use modules::Adaptors::Source;
use modules::Adaptors::Source_Group;
use modules::Adaptors::Release_File;
use modules::Adaptors::Lane;
use modules::Adaptors::Group_Summary;
use modules::Adaptors::Pipeline_Step;
use modules::Adaptors::Human_Related_Sample;
use modules::SystemCall;
use modules::SummaryReport;
use modules::Utils;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "source_list=s",
	   "source=s",
	   "sample_list=s",
	   "project=s",
	   "gene_list=s",
	   "coord=s",
	   "no_fails",
	   "no_db",
	   "production",
	   "chrom=s",
	   "no_pileups",
	   "all_related",
	   "all_db",
	   "denovo",
	   "min_phase_block=i",
	   "phase_var_num=i",
	   "inheritance=s",
	   "comhet",
	   "max_allele_freq=s",
	   "sift=s",
	   "polyphen=s",
	   "min_num_aff=i",
	   "debug=s",
	   	"score",
		"report_xml=s",
		"skip_family",
		"sv"
	   ) || modules::Exception->throw("Invalid command-line option for script\n");;

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

summarize_report.pl -source input_single_source_name_on_command_line -min_num_aff filter_on_number_of_affected_samples_containing_variants -polyphen filter_on_polyphen_category(default=probably_damaging) -sift filter_on_sift_category(default=deleterious) -comhet only_report_comhet -max_allele_freq filter_on_allele_freq -phase_var_num min_number_variants_to_make_block(default=2) -min_phase_block min_phase_block_size(default=10kb) -denovo only_report_denovo -inheritance filter_on_inheritance_type(options=ad,ar,xd,xr) -all_related run_on_all_related -coord genomic_coordinates(format=chr:start-end) -all_db run_on_all_in_db -gene_list only_report_on_gene_in_file -chrom only_report_chr -no_pileups skip_reporting_pileup_data -source_list list_of_sources(default=all) -production production_run -sample_list file_with_list_of_samples_for_report -project project_to_combine_data_from -single_source_group run_on_first_source_group_only -no_fails ignore_fails -skip_family skip_family_for_combined_cohort_cases -sv run_for_SVs [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

summarize_report.pl -> generate report summary for cohorts, projects, etc

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./summarize_report.pl 

=cut

my $clus_conf = modules::Pipeline::get_cluster_conf();

my $report_conf = modules::Pipeline::get_report_conf();

my $pipe_config = modules::Pipeline::get_pipe_conf();

#Default is to include failed snvs
my $fails = defined $OPT{no_fails} ? 0 : 1;

my $sv = defined $OPT{sv}?1:0; #for SVs

if ($OPT{all_related}) {
	modules::Exception->throw("ERROR: Can't use this option SV");
}


#production flag used for setting directories
my $production = defined $OPT{production}?1:0;

#Determines whether related samples are being analysed as extra things are reported; default of yes until proven otherwise, mixed data gets set to no
my $family_flag = defined $OPT{skip_family}?0:1; 

#Records the type of report (project,source, all_related, all_db or sample)
my $report_type = 'all_db';

#Used when report_type is sample
my %included_samples = ();

#By default use the results dirs to get the all candidate sources;  @sources is reset if -source_list or -project or -sample_list is used
my @sources = ();

#Get the source_names and set report_type variable as this affects behavior later on
if ($OPT{source_list}) {
	$report_type = 'source';
	open(SOURCES,$OPT{source_list}) || modules::Exception->throw("Can't open file $OPT{source_list}\n");
	while (<SOURCES>) {
		my ($source_name) = $_ =~ /(\S+)/;
		push @sources, $source_name;
		my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source_name);
		if (!defined $source_obj) {
			modules::Exception->throw("ERROR: Can't retrieve source from db for $source_name");
		}
		
	}
} elsif ($OPT{source}) {
	$report_type = 'source';
	push @sources, $OPT{source};
	my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$OPT{source});
	if (!defined $source_obj) {
		modules::Exception->throw("ERROR: Can't retrieve source from db for $OPT{source}");
	}
} elsif ($OPT{project}) {
	$report_type = 'project';
	
	#here we use the db project entry 
	my ($project_obj) = modules::Adaptors::Project->search(project_name=>$OPT{project});
	if (!defined $project_obj) {
		modules::Exception->throw("ERROR: Can't retrieve project from db for $OPT{project}");
	}
	my @source_objs = modules::Adaptors::Source->search(project_id=>$project_obj->id);
	if (!@source_objs) {
		modules::Exception->throw("ERROR: Couldn't retrieve any sources for project $OPT{project}");
	}		
	
	for my $source_obj (@source_objs) {
		push @sources, $source_obj->external_source_name;
	}
	
	
} elsif ($OPT{sample_list}) {
	$report_type = 'sample';
	
	open(SAMPLES,$OPT{sample_list}) || modules::Exception->throw("Can't open file $OPT{sample_list}\n");
        while (<SAMPLES>) {
                my ($sample_name) = $_ =~ /(\S+)/;
                $included_samples{$sample_name} = 1;

                my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);
                if (!defined $sample_obj) {
                        modules::Exception->throw("ERROR: Can't retrieve sample from db for $sample_name");
                }

                my ($source_group_obj) = modules::Adaptors::Source_Group->search(id=>$sample_obj->source_group_id);
                if (!defined $source_group_obj) {
                        modules::Exception->throw("ERROR: Can't retrieve source_group from db for $sample_name");
                }

                my @source_objs = modules::Adaptors::Source->search(id=>$source_group_obj->source_id);
                if (!@source_objs) {
                        modules::Exception->throw("ERROR: Couldn't retrieve source form db for $sample_name");
                }

                for my $source_obj (@source_objs) {
                        push(@sources, $source_obj->external_source_name) unless grep{$_ eq $source_obj->external_source_name} @sources;
                }

        }
} elsif ($OPT{all_related}) {
	
	$report_type = 'all_related';
	my @possible_source_types = qw(human_related_gatk);
	
	for my $pos_source_type ( @possible_source_types ) {
	    my $result_dir = $clus_conf->read($pos_source_type,'base_directories','base_results_directory'); 
	    if (!-d $result_dir) {
			modules::Exception->throw("Directory $result_dir problem");
		}
	    opendir(DIR,$result_dir) || modules::Exception->throw("Can't open directory $result_dir");
		my @tmp_sources = grep {/\w/} readdir DIR;
		push @sources, @tmp_sources; 
	}
} else {
	my @possible_source_types = qw(human_related_gatk human_single_gatk);
	
	for my $pos_source_type ( @possible_source_types ) {
	    my $result_dir = $clus_conf->read($pos_source_type,'base_directories','base_results_directory'); 
	    if (!-d $result_dir) {
			modules::Exception->throw("Directory $result_dir problem");
		}
	    opendir(DIR,$result_dir) || modules::Exception->throw("Can't open directory $result_dir");
		my @tmp_sources = grep {/\w/} readdir DIR;
		push @sources, @tmp_sources; 
	}
}

if (!@sources) {
	modules::Exception->throw("ERROR: No sources to analyse");
}

#Iterate through sources to determine if there are non related samples in source list -> this sets family flag to 0
for my $source_name (@sources) {
	my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source_name);
	my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
		
	if (! @source_group_objs ) {
		modules::Exception->throw("ERROR: Can't retrieve source_group with source id $source_obj->id");
	}
	
	for my $source_group_obj (@source_group_objs) {
		my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);	
		
		my $complete_sample_count = 0;
		
		for my $sample_obj ( @samples ) {
			my $sample_type = $sample_obj->sample_type;
			if ($sample_type !~ /related/) {
				$family_flag = 0;  #Don't do family reporting if the samples are mixed
			}
		}
	}
	
	
}

if (!$family_flag) {
	#Check inapprorpriate filters aren't turned on
	if (defined $OPT{min_num_aff} || defined $OPT{inheritance} || defined $OPT{comhet} || defined $OPT{denovo} || defined $OPT{min_phase_block}) {
		modules::Exception->throw("Cannot utilise filters for pedigrees without pedigrees being analysed");
	}
}


#db flag; write output reports to database; default is to write to db
my $write_db = defined $OPT{no_db}?0:1;
$write_db = 0 if $report_type eq 'all_db';

my $pileups = defined $OPT{no_pileups}?0:1;
$pileups = 0 if $report_type eq 'all_db'; #never want to do this for the whole database
$pileups = 0 if $sv;
#$family_flag = 0 if $sv;

my (undef,$snv_gene_col_name) = split(",",$report_conf->read('snv','common','human_file_annotations'));
my (undef,$indel_gene_col_name) = split(",",$report_conf->read('indel','common','human_file_annotations'));

#default args for object
my %args = (
			-sources=>\@sources,
			-family=>$family_flag,
			-fail=>$fails,
			-report_type=>$report_type,
			-production=>$production,
			-indel_gene_col_name=>$indel_gene_col_name,
			-snv_gene_col_name=>$snv_gene_col_name,
			-pileups=>$pileups #For large projects we don't want to calculate every allele
			);

if ($sv) {
	$args{-sv} = $sv;
	if ($OPT{min_phase_block} || $OPT{phase_var_num}  || $OPT{polyphen} || $OPT{sift}) {
		modules::Exception->throw("Not a suitable feature for SVs");
	}
} 

$args{-sv_gene_col_name} = $snv_gene_col_name; #Same as snvs for gene name

#Finally deal with all the optional filtering 

#Filtering by gene
if (defined $OPT{gene_list}) {
	my %genes = ();
	open(GENES,$OPT{gene_list}) || modules::Exception->throw("Can't open file $OPT{gene_list}\n");
	while (<GENES>) {
		next unless /\S/;
		my ($gene_name) = $_ =~ /(\S+)/;
		$genes{$gene_name} = 1;
	}
	$args{-gene_list} = \%genes;
}

if (defined $OPT{chrom} && defined $OPT{coord}) {
	modules::Exception->throw("ERROR: Only filter by -chrom or -coord, not both");
}

#Filter by chrom
if (defined $OPT{chrom}) {
	$args{-chrom} = $OPT{chrom};
}

#Filter by chrom
if (defined $OPT{coord}) {
	if ($OPT{coord} =~ /[0-9XY]+:\d+\-\d+/) {
		$args{-coord} = $OPT{coord};
	} elsif ($OPT{coord} =~ /[0-9XY]+:(\d+)/) {
		my $coord = $OPT{coord} .'-'.$1;
		$args{-coord} = $coord;		
	} else {
		modules::Exception->throw("ERROR: -coord argument must be of the form chr:start-end OR chr:coord");
	}
}

if (defined $OPT{min_phase_block}) {
	$args{-phase_block_size} = $OPT{min_phase_block}
} 

if (defined $OPT{phase_var_num}) {
	$args{-phase_var_num} = $OPT{phase_var_num};
} 


if (defined $OPT{denovo} && $OPT{inheritance} ) {
	modules::Exception->throw("ERROR: Only filter by -denovo or -inheritance. not both");
}

if (defined $OPT{denovo}) {
	$args{-denovo} = 1;
}

my @dis_inh_options = qw(ad ar xd xr);
if (defined $OPT{inheritance}) {
	my $input = lc($OPT{inheritance});
	
	my $match = 0;
	for my $dis_type ( @dis_inh_options ) {
	    if ($dis_type eq $input) {
	    	$match = 1;
	    }
	}
	
	
		
	if (!$match) {
		modules::Exception->throw("ERROR: problem with inheritance; must be ad,ar,xd,or xr");
	} else {
		my %inheritance_map = (
								ad => 'auto-dominant',
								ar => 'auto-recessive',
								xd => 'x-dominant',
								xr => 'x-recessive'
							  );
		$args{-inheritance} = $inheritance_map{$input};
	}
}

if ($OPT{comhet}) {
	$args{-comhet} = 1;
}

if ($OPT{max_allele_freq}) {
	if ($OPT{max_allele_freq} >= 0 && $OPT{max_allele_freq} <= 1) {
		$args{-max_allele_freq} = $OPT{max_allele_freq};
	} else {
		modules::Exception->throw("ERROR: max_allele_freq must be between 0 and 1");
	}
}

if ($OPT{min_num_aff}) {
	$args{-min_num_aff} = $OPT{min_num_aff};
}


if ($OPT{polyphen}) {
	my $match = 0;
	my @polyphen_options = qw(probably_damaging possibly_damaging benign);
	for my $polyopt ( @polyphen_options ) {
	    if ($polyopt eq $OPT{polyphen}) {
	    	$match = 1;
	    	$args{-polyphen} = $OPT{polyphen};
	    }
	}
	if (!$match) {
		modules::Exception->throw("ERROR: polyphen input must be probably_damaging,possibly_damaging, or benign");
	}
}

if ($OPT{sift}) {
	my @sift_options = qw(deleterious tolerant);
	my $match = 0;
	for my $siftopt ( @sift_options ) {
	    if ($siftopt eq $OPT{sift}) {
	    	$match = 1;
	    	$args{-sift} = $OPT{sift};
	    }
	}
	if (!$match) {
		modules::Exception->throw("ERROR: sift input must be tolerated or deleterious");
	}
}


if (defined $OPT{score}){
	if(defined $OPT{report_xml}) {
		$args{-score} = 1;
		$args{-report_xml} = $OPT{report_xml}; 
	} else{
		modules::Exception->throw("ERROR: provide -report_xml file to use option -score");
	}
}


my $summary_reports = modules::SummaryReport->new(%args);
												
my $sample_to_process;												
												
#Load up all the sample info
if ($report_type eq 'sample') {
	$sample_to_process = $summary_reports->set_sample_info(\%included_samples);
} else {
	$sample_to_process = $summary_reports->set_sample_info('all');
}

print Dumper $summary_reports;


if (!$sample_to_process) {
	modules::Exception->warning("No samples available for analysis");
	exit;
}

print "Set master headers\n";
$summary_reports->set_master_headers();

print "Parse reports\n";
if (!$summary_reports->parse_reports()) {
	modules::Exception->throw("Error: found no lines that match your criteria");
}



if ($pileups) {
	print "Generate pileups\n";
	$summary_reports->generate_pileups();
}

if ($sv) {
	$summary_reports->generate_sv_report();
	$summary_reports->write_to_db() if $write_db;
	exit;
}

#Only for pedigrees
if ($family_flag && !$sv) {
	print "Get parent alleles\n";
	$summary_reports->get_parent_alleles();
	
	print "Get Com het\n";
	$summary_reports->get_compound_het();
	
	print "Get Phasing\n";
	$summary_reports->get_phasing_blocks();
}



print "Get line data\n";
$summary_reports->generate_line_data();

print "Write files\n";
$summary_reports->write_to_files();

$summary_reports->write_to_db() if $write_db;

#mdss results if necessary
my $cluster_config = modules::Pipeline::get_cluster_conf();

my $mdss = $cluster_config->read('common','mdss','mdss_flag');
if ($mdss && $write_db && $family_flag) {
	$summary_reports->mdss();
	
}

