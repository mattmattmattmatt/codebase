#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Cluster;
use modules::Pipeline;
use modules::PipelineSample;
use modules::Adaptors::Project;
use modules::Adaptors::Source;
use modules::Adaptors::Source_Group;
use modules::Adaptors::Sample;
use modules::Adaptors::Human_Related_Sample;
use modules::Adaptors::Lane;
use modules::Adaptors::Read_File;
use modules::Adaptors::Run;
use Pod::Usage;
use Cwd;

use vars qw(%OPT);


GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"external_name=s",
	   		"affected_readdirs=s",
	   		"unaffected_readdirs=s",
	   		"affected_relations=s",
	   		"unaffected_relations=s",
	   		"proband_sex=s",
	   		"unaffected_sex=s",
	   		"sequencing_centre=s",
	   		"read1_pattern=s",
	   		"cluster_xml=s",
	   		"steps_xml=s",
	   		"skip_quality",
	   		"new_source_group",
	   		"test",
	   		"only_qsubs",
	   		"source_group_number=i",
	   		"project_name=s",
	   		"no_unaffected",
	   		"sequence_type=s",
	   		"affected_sample_external_names=s",
	   		"unaffected_sample_external_names=s",
	   		"no_mdss",
			"submit",
			"gatk",
			"apf_req_id=s",
			"update_source",
			"new_member_external_names=s",
			"score",
			"report_xml=s"     
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{external_name} || !$OPT{affected_readdirs} || !$OPT{affected_relations});

	   
=pod

=head1 SYNOPSIS

add_cohort_sample.pl -new_source_group new_source_group_to_existing_patient(default=new_patient) -sequence_type sequence_experiment_type (exome,genome, or targeted; default=exome) -external_name external_name(e.g.tb412) -read1_pattern read1_regex_pattern(default=R1) -steps_xml xml_steps_file(default=../conf/steps.xml) -cluster_xml xml_cluster_file(default=../conf/cluster.xml) -sequencing_center sequencing_centre(option=BRF,AGRF,RAM;default=AGRF) -proband_sex (male,female,default=unknown) -affected_relations comma_delim_sample_type -unaffected_relations comma_delim_sample_type  -affected_readdirs comma_delim_fq_directories -unaffected_readdirs comma_delim_fq_directories -skip_quality Skip_read_check_and_assume_phred33_qualities -only_qsubs only_create_qsubs -source_group_number source_group_number_for_only_qsubs_flag -sample_external_names comma_delim_external_names(eg APOP64,APOP63) [options]

# eg of adding an extra family member to an existing source : use -update_source and -new_member_external_names

add_cohort_sample.pl -project_name 'MONA' -external_name MONA_cohort1 -read1_pattern R1 -sequencing_centre Macrogen -sequence_type genome -affected_sample_external_names 'MONA4,MONA1,MONA5' -affected_readdirs 'MONA_cohort1_affected_sister,MONA_cohort1_affected_proband,MONA_cohort1_affected_other' -affected_relations 'sister,proband,other' -unaffected_sample_external_names 'MONA3,MONA2' -unaffected_readdirs 'MONA_cohort1_unaffected_mother,MONA_cohort1_unaffected_father' -unaffected_relations 'mother,father' -proband_sex female -apf_req_id 227^227^257^227^227 -no_mdss -update_source -new_member_external_names MONA5

Required flags: -external_name -affected_readdirs  -affected_relations ((-unaffected_readdirs && -unaffected_relations) || -no_unaffected) 

    -help  brief help message

    -man   full documentation

=head1 NAME

add_cohort_sample.pl -> Add new family cohort

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE


=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

if ($svndir =~ /trunk/) {
	modules::Exception->throw("ERROR: Do not use trunk for svn directory") unless $OPT{test};
}

my $unaffected;

if ($OPT{no_unaffected}) {
	if (defined $OPT{unaffected_readdirs}  || defined $OPT{unaffected_relations}) {
		modules::Exception->throw("ERROR: Can't pass -no_unaffected and other unaffected_flags");
	}
	$unaffected = 0;
} else {
	if (!defined $OPT{unaffected_readdirs} && !defined $OPT{unaffected_relations}) {
		modules::Exception->throw("ERROR: Need -unaffected_readdirs and -unaffected_relations");
	}
	$unaffected = 1;
}

if($OPT{update_source}){	#update existing source to add additional members of the cohort, need to provide sample ids of new members!
	if (!defined $OPT{new_member_external_names}){
		modules::Exception->throw("ERROR: need to provide -new_member_external_names");
	}
}

my $seq_type = defined $OPT{sequence_type}?$OPT{sequence_type}:'exome';


my $report_conf = modules::Pipeline::get_report_conf();
my $writedb = defined $OPT{only_qsubs}?0:1;

#Get the arguments
my $external_name = $OPT{external_name};

my $parallel = 1; #Always use this mode...
my $project_name = defined $OPT{project_name}?$OPT{project_name}:'Melanoma';
my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
if (!defined $project_obj) {
	modules::Exception->throw("ERROR: Can't find project name $project_name");
}

my $source_type = 'human_related';
my $use_gatk = defined $OPT{gatk}?1:0;
$source_type .= '_gatk' if $use_gatk;

my $cluster_config = modules::Pipeline::get_cluster_conf();
my $pipe_config = modules::Pipeline::get_pipe_conf();


my $read_directory = $cluster_config->read($source_type,'base_directories','base_read_directory').'/'.$external_name;
my $read1_regex = defined $OPT{read1_pattern}?$OPT{read1_pattern}:'R1';

#Read the cluster specific variables 
my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory').'/'.$external_name;
system("mkdir $qsub_dir") if (!-d $qsub_dir);

if ( !-d $read_directory ) {
	modules::Exception->throw("File $read_directory doesn't exist");	
}
my $outdir = $cluster_config->read($source_type,'base_directories','base_run_directory').'/'.$external_name;
system("mkdir -p $outdir") if (!-d $outdir);

my $bam_link_dir = $outdir . '/bam_links';
system("mkdir -p $bam_link_dir") if (!-d $bam_link_dir);


my $threadnum = $cluster_config->read('common','qsub_vars','thread_num');
my $scheduler = $cluster_config->read('common','scheduler');


my $mdss = $cluster_config->read('common','mdss','mdss_flag');
my $sys_call = modules::SystemCall->new();

if ($mdss && $writedb) {
	#Create a backup read files in this case
	my $mdss_read_dir = $cluster_config->read('common','mdss','mdss_reads');
	
	#marcin
	#my $mdss_command = "mdss put -r $read_directory $mdss_read_dir";
	my $mdss_command = "$svndir/scripts/mdss_put.sh $read_directory $mdss_read_dir";
	$sys_call->run($mdss_command) unless $OPT{no_mdss};
}

#Check reads match if alread
if ( !-d $read_directory ) {
	modules::Exception->throw("File $read_directory doesn't exist");	
}

my @unaffected_relations = ();
if ($unaffected) {
	@unaffected_relations = split(",",$OPT{unaffected_relations});
}
my @affected_relations = split(",",$OPT{affected_relations});

my %possible_relations = map {$_ => 1} split(",", $pipe_config->read('common','related_types'));

#Check all the sample_types
for my $relation_type (@unaffected_relations, @affected_relations) {
	if (!exists $possible_relations{$relation_type}) {
		modules::Exception->throw("ERROR: $relation_type is not a correct relation_type option");
	}    
}



my @unaffected_readdirs = ();
if ($unaffected) {
	@unaffected_readdirs = split(",",$OPT{unaffected_readdirs});
	if (@unaffected_readdirs != @unaffected_relations) {
		modules::Exception->throw("ERROR: Number of arguments must match for readdirs and relations");
	}
}
my @affected_readdirs = split(",",$OPT{affected_readdirs});

for my $readdir (@affected_readdirs, @unaffected_readdirs) {
	my $output = `ls $read_directory/$readdir/*$read1_regex* 2>/dev/null`;
	if (!$output && $writedb) {
		modules::Exception->throw("ERROR: No read files match  $read_directory/$readdir/*$read1_regex*");
	}
	
}

my @new_member_external_names;

if($OPT{new_member_external_names}){
	@new_member_external_names = split(",",$OPT{new_member_external_names});	
}

#Create the pipeline object
#Get the xml files to get variables we need first
my $steps_xml = defined $OPT{steps_xml}?$OPT{steps_xml}:"$svndir/conf/".$source_type.".steps.xml";

if ( !-e $steps_xml ) {
	modules::Exception->throw("File $steps_xml doesn't exist");	
}



my %sequencing_centres = map{ $_ => 1 } split(",",$pipe_config->read('common','sequencing_centres'));

my $sequencing_centre = defined $OPT{sequencing_centre}?$OPT{sequencing_centre}:'AGRF';
if (!exists $sequencing_centres{$sequencing_centre}) {
	my $seq_centre_str = join (',',keys %sequencing_centres);
#	modules::Exception->throw("ERROR: Sequencing centre must be $seq_centre_str from pipe.xml");
}

my @affected_external_sample_names;

if (defined $OPT{affected_sample_external_names}) {
	(my $remove_spaces = $OPT{affected_sample_external_names}) =~ s/ /_/g;
	@affected_external_sample_names = split(",",$remove_spaces);
	if (@affected_relations != @affected_external_sample_names) {
		modules::Exception->throw("ERROR: Number of arguments must match for readdirs and external_name");
	}
}

my @unaffected_external_sample_names;

if ($unaffected) {
	if (defined $OPT{unaffected_sample_external_names}) {
		(my $remove_spaces = $OPT{unaffected_sample_external_names}) =~ s/ /_/g;
		@unaffected_external_sample_names = split(",",$remove_spaces);
		if (@unaffected_readdirs != @unaffected_external_sample_names) {
			modules::Exception->throw("ERROR: Number of arguments must match for readdirs and external_name");
		}
	}
}

my $proband_sex = defined $OPT{proband_sex}?$OPT{proband_sex}:'unknown';

my $normal_num = @unaffected_readdirs;
my $tumour_num = @affected_readdirs;


#Check all the argument numbers match
if (@affected_relations != @affected_readdirs) {
	modules::Exception->throw("ERROR: Number of arguments must match for readdirs and relations");
}

#First create the patient entry
my %source_info = (
					external_source_name=>$external_name,
					project_id=>$project_obj->id,
					source_type=>$source_type
				   );

my $source_db_id;
my ($source_obj) = modules::Adaptors::Source->search('external_source_name'=>$external_name);


if ($OPT{new_source_group}) {
	$source_db_id = $source_obj->id;
	print STDERR "Add sample group to source $external_name with id $source_db_id\n";
} elsif ($writedb) {
	#normal case; writing out new objects
	if (defined $source_obj && !defined $OPT{update_source}) {
		modules::Exception->throw("ERROR: A source already exists with this external_name ($external_name). Must delete source first");
	}
	if (defined $OPT{update_source}){
		$source_db_id = $source_obj->id;
		print STDERR "Source $source_obj->external_name is found! Adding new members for this source\n";
	}else {
		#normal case here; insert objects into database
		$source_db_id = modules::Adaptors::Source->insert(\%source_info);
		($source_obj) = modules::Adaptors::Source->search('id'=>$source_db_id);
		print STDERR "Create source $external_name with id $source_db_id\n";
	}
} else {
	$source_db_id = $source_obj->id;
}



if (!defined $source_obj) {
	modules::Exception->throw("ERROR: Can't retrieve source object for $external_name");
}

my $source_group_count;
#If we're adding a source_group to existing patient
if ($OPT{new_source_group}) {
	if ($writedb) {
		my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
		if (!@source_group_objs) {
			modules::Exception->throw("ERROR: No existing source_groups for $external_name");
		}
		#Get the highest sample group
		my $max_sg = 1;
		for my $source_group_obj ( @source_group_objs ) {
		    if ($source_group_obj->source_group_number > $max_sg) {
		    	$max_sg = $source_group_obj->source_group_number;
		    }
		}
		$source_group_count = $max_sg  + 1;
	} else {
		if (! defined $OPT{source_group_number}) {
			modules::Exception->throw("ERROR: Need source_group_number input for this sample");
		} else {
			$source_group_count = $OPT{source_group_number};
		}
	}
} elsif ($writedb) {
	$source_group_count = 1;
} else {
	#Writing out qsubs in patient; confirm there aren't multiple source_groups
	my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
	if (@source_group_objs > 1) {
		modules::Exception->throw("ERROR: Multiple source groups for source $external_name");
	}
	$source_group_count = 1;
}

my $source_group_name = $external_name.'_sg'.$source_group_count;
my $snv_call_dir = $outdir .'/' .$source_group_name.'_snvcalls';
system("mkdir -p $snv_call_dir") if (!-d $snv_call_dir);

#my $first_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first');
#my $first_lane_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first_lane');
#my $last_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'last');
my $first_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first'); 
my $last_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'last');
my $last_align_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'finish_align');
my $parvar_first = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_first');
my $parvar_merge = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_merge');
my $nonpar_start = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'start_non_parvar');
my $first_lane_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first_lane');

#Next create the sample group
my $sample_total = @affected_readdirs + @unaffected_readdirs;
my @all_samples = @affected_external_sample_names;
push(@all_samples, @unaffected_external_sample_names);

my %new_samples_mapper;
if($OPT{new_member_external_names}){
	for my $sample (@new_member_external_names){
		$new_samples_mapper{$sample}=1;
	} 
}


#Optional apf request id field, applies to samples from APF only
#Check apf_req_id numbers match then map ids
my %apf_req_id_mapper;

if (defined $OPT{apf_req_id}){
	my $apf_req_id_full = $OPT{apf_req_id}; 
	my @apf_req_ids = split('\^',$apf_req_id_full);

	if(@apf_req_ids != @all_samples) {
		modules::Exception->throw("ERROR: Number of apf_req_ids must match for number of samples");
	}else{
		@apf_req_id_mapper{@all_samples} = @apf_req_ids;
	}
}

my $score = 0;
my $report_xml;

if(defined $OPT{score}){
	if(!defined $OPT{report_xml}){
		modules::Exception->throw("ERROR: report_xml containing scoring scheme musth be provided with -score option");
	}
	$report_xml = $OPT{report_xml};	
	$score = 1;	
}

my $sourcegroup_db_id;
if ($writedb && !defined $OPT{update_source}) {
	my %source_group_info = (
							total_samples=>$sample_total,
							source_group_number=>$source_group_count,
							source_group_name=>$source_group_name,
							source_id=>$source_db_id
							);
	$sourcegroup_db_id = modules::Adaptors::Source_Group->insert(\%source_group_info);
	print STDERR "\tCreate source_group $source_group_name with id $sourcegroup_db_id\n";
} else {
	my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$source_group_name,source_id=>$source_db_id);
	if (!defined $source_group_obj) {
		modules::Exception->throw("No source_group object for $source_group_name and source id $source_db_id");
	}
	$sourcegroup_db_id = $source_group_obj->id;
}



my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
my $pipeline_sample = modules::PipelineSample->new(-source_name=>$external_name,-cluster_obj=>$cluster_obj,-source_type=>$source_type);

#Next create the samples

#Acceptable suffices
my %suffices = (
				gz=>1,
				bz2=>1
				);

my %qsub_commands = ();
#Now create the samples..
#First unaffected
&Create_Samples('unaffected',$normal_num) if $unaffected;

#Next affected
&Create_Samples('affected',$tumour_num);

if (keys %qsub_commands) {
	print STDERR "Run the following commands...\n";
	for my $command ( keys %qsub_commands ) {
	    print "\t$command\n";
	    system("$command") if $OPT{submit};
	}
	print "\n";
}

my @commands = ();
push @commands, "$svndir/scripts/summarize_reports.pl -source $external_name -production -sv";


#create the cohort summary file
if($score){
	push @commands, "$svndir/scripts/summarize_reports.pl -source $external_name -production -score -report_xml $report_xml";
}else{
	push @commands, "$svndir/scripts/summarize_reports.pl -source $external_name -production";
}

my $qsub_file = $external_name .'_cohort.qsub';
$cluster_obj->create_single_qsub(-qsub_dir=>"$qsub_dir/runs/",-qsub_file=>$qsub_file,-commands=>\@commands,-queue=>'copyq',-walltime=>'10:00:00');


#Does all the work creating samples, lanes, and read_files
sub Create_Samples () {
	my ($affected_unaffected,$sample_total) = @_;	
	my $affected_flag = $affected_unaffected eq 'affected'?1:0;

	for ( my $sample_count = 0 ; $sample_count < $sample_total ; $sample_count++ ) {
	   
		my $sample_number = $sample_count+1;
	
		my $sample_name = $source_group_name.'_'.$affected_unaffected.$sample_number;
		#if ($parallel) {
        	#	$pipeline_sample->create_parallel_snvcall_qsubs(-sample=>$sample_name,-source_group_name=>$source_group_name);
		#}	
		#Figure out how many lanes there are
		my $lane_data;
		
		#Get the local sample_directory for getting reads later
		my $sample_readdir;
			
		
		
		if ($affected_flag) {
			$sample_readdir = $affected_readdirs[$sample_count];
		} else {
			$sample_readdir = $unaffected_readdirs[$sample_count];
			#($lane_data) = &Get_Lane_Info($unaffected_readdirs[$sample_count]);
		}
		
#		($lane_data) = modules::Pipeline::get_lane_info($read_directory .'/'.$sample_readdir,$read1_regex);
		
	    my $lane_total;# = keys %{$lane_data};
	    my $sample_db_id;
	    if ($writedb) {

			($lane_data) = modules::Pipeline::get_lane_info($read_directory .'/'.$sample_readdir,$read1_regex);
		
		    $lane_total = keys %{$lane_data};

		    my %sample_info = (
		    					sample_number=>$sample_number,
		    					total_lanes => $lane_total,
		    					sample_name => $sample_name,
		    					source_group_id=>$sourcegroup_db_id,
		    					sequence_type=>$seq_type,
		    					);
		    
		    if ($affected_flag) {
		    	$sample_info{sample_type} = 'related_affected';
		    	if (@affected_external_sample_names > $sample_count) {
		    		$sample_info{external_sample_name} = $affected_external_sample_names[$sample_count];
		    	} else {
		    		$sample_info{external_sample_name} = 'NONE';
		    	}
			if(defined $OPT{apf_req_id}){
				$sample_info{apf_request_id} = $apf_req_id_mapper{$affected_external_sample_names[$sample_count]};
			}
		    } else {
		    	$sample_info{sample_type} = 'related_unaffected';
		    	
		    	if (@unaffected_external_sample_names > $sample_count) {
		    		$sample_info{external_sample_name} = $unaffected_external_sample_names[$sample_count];
		    	} else {
		    		$sample_info{external_sample_name} = 'NONE';
		    	}
			if(defined $OPT{apf_req_id}){
				$sample_info{apf_request_id} = $apf_req_id_mapper{$unaffected_external_sample_names[$sample_count]};
			}
		    }

		    if($OPT{update_source} && !exists $new_samples_mapper{$sample_info{external_sample_name}}){
			next;
		    }
		    $sample_db_id =  modules::Adaptors::Sample->insert(\%sample_info);
		    print STDERR "\t\tCreate sample $sample_name with id $sample_db_id\n";   
	    
		    my %related_sample_info = (
		    						affected => $affected_flag,
		    						sample_id=>$sample_db_id
		    						);
		    
		    
		    if ($affected_flag) {
		    	
		    	
				$related_sample_info{relation} = $affected_relations[$sample_count];
			    if ($affected_relations[$sample_count] eq 'proband') {
			    	$related_sample_info{sex} = $proband_sex;
			    } elsif ($affected_relations[$sample_count] eq 'father' || $affected_relations[$sample_count] eq 'brother' || $affected_relations[$sample_count] eq 'grandfather' || $affected_relations[$sample_count] eq 'son' || $affected_relations[$sample_count] eq 'uncle' || $affected_relations[$sample_count] eq 'nephew') {
			    	$related_sample_info{sex} = 'male';
			    } elsif ($affected_relations[$sample_count] eq 'mother' || $affected_relations[$sample_count] eq 'sister'  || $affected_relations[$sample_count] eq 'grandmother' || $affected_relations[$sample_count] eq 'daughter' || $affected_relations[$sample_count] eq 'aunt' || $affected_relations[$sample_count] eq 'niece' ) {
			    	$related_sample_info{sex} = 'female';
			    } elsif ($affected_relations[$sample_count] eq 'other' || $affected_relations[$sample_count] eq 'sibling' || $affected_relations[$sample_count] eq 'parent' || $affected_relations[$sample_count] eq 'child' || $affected_relations[$sample_count] eq 'cousin' || $affected_relations[$sample_count] eq 'spouse') {
			    	$related_sample_info{sex} = 'unknown';
			    } else {
			    	modules::Exception->throw("ERROR: Can't determine gender for $sample_name");
			    }
		    } else {
		    	
		    	
		    	$related_sample_info{relation} = $unaffected_relations[$sample_count];
		    	if ($unaffected_relations[$sample_count] eq 'proband') {
					modules::Exception->throw("ERROR: Proband cannot be unaffected");
			    } elsif ($unaffected_relations[$sample_count] eq 'father' || $unaffected_relations[$sample_count] eq 'brother' || $unaffected_relations[$sample_count] eq 'grandfather' || $unaffected_relations[$sample_count] eq 'son' || $unaffected_relations[$sample_count] eq 'uncle' || $unaffected_relations[$sample_count] eq 'nephew' ) {
			    	$related_sample_info{sex} = 'male';
			    } elsif ($unaffected_relations[$sample_count] eq 'mother' || $unaffected_relations[$sample_count] eq 'sister'  || $unaffected_relations[$sample_count] eq 'grandmother' || $unaffected_relations[$sample_count] eq 'daughter' || $unaffected_relations[$sample_count] eq 'aunt' || $unaffected_relations[$sample_count] eq 'niece' ) {
			    	$related_sample_info{sex} = 'female';
			    } elsif ($unaffected_relations[$sample_count] eq 'other' || $unaffected_relations[$sample_count] eq 'sibling' || $unaffected_relations[$sample_count] eq 'parent' || $unaffected_relations[$sample_count] eq 'child'|| $unaffected_relations[$sample_count] eq 'cousin' || $unaffected_relations[$sample_count] eq 'spouse') {
			    	$related_sample_info{sex} = 'unknown';
			    } else {
			    	modules::Exception->throw("ERROR: Can't determine gender for $sample_name");
			    }
		    	
		    }
		     	
		    my $related_sample_db_obj =  modules::Adaptors::Human_Related_Sample->insert(\%related_sample_info);
		    my $related_sample_id = $related_sample_db_obj->id;
		    print STDERR "\t\t\tCreate related_sample $sample_name with id $related_sample_id\n";   
		     	
	    } else {
	    	#Get the existing sample
	    	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name,source_group_id=>$sourcegroup_db_id);
	    	print "\t\t\tSample $sample_name source_groupdb $sourcegroup_db_id\n";
	    	$sample_db_id = $sample_obj->id;
	    	
	    	$lane_total = $sample_obj->total_lanes;
	    }				
	    
	    
	    my @bams = ();			
	    my $lanes_outdir = $outdir.'/'.$sample_name.'_lanes';
	    my $runs_outdir = $outdir.'/'.$sample_name.'_runs';
	    
	    my $encoding;
		#Create the lanes for each sample
	    for ( my $lane_count = 0 ; $lane_count < $lane_total ; $lane_count++ ) {

	        my $lane_db_id;
			my $lane_bam;
			my %commands = ();

	    	my $lane_number = $lane_count;
	    	$lane_number++;
	    	my $lane_name = $sample_name.'_l'.$lane_number;
	        
	        #Get the read_file info
	        my $read1_compressed = 0;
	        my $read2_compressed = 0;
	        my $read1_suffix;
	        my $read2_suffix;
			my $read1_name = $lane_name.'_r1';
			my $read2_name = $lane_name.'_r2';
			
			my $read1_full_file;
			my $read2_full_file;
				
			#By default don't check the quality (assume phred33)
			my $quality = defined $OPT{skip_quality}?1:0;
								
	        if ($writedb) {					
				my $read1_actual = $lane_data->{$lane_count}{read1};
				my $read2_actual = $lane_data->{$lane_count}{read2};
				
				$read1_full_file = $read_directory .'/'. $sample_readdir. '/'. $lane_data->{$lane_count}{read1};
				$read2_full_file = $read_directory .'/' . $sample_readdir. '/'. $lane_data->{$lane_count}{read2};
				my $read1_symlink = $read_directory .'/' . $sample_readdir. '/'. $read1_name;
				my $read2_symlink = $read_directory .'/' . $sample_readdir. '/'. $read2_name;
				
				
				my $read_directory_full = $read_directory .'/'. $sample_readdir;
				
				
				#First check the read files exist
				if ( !-e $read1_full_file ) {
					modules::Exception->throw("File $read1_full_file doesn't exist");	
				}	
				if ( !-e $read2_full_file ) {
					modules::Exception->throw("File $read2_full_file doesn't exist");	
				}	
				
				#Create symlinks to the files if they don't exist
				if (!-e $read1_symlink) {
					system("cd $read_directory/$sample_readdir; ln -s $read1_actual $read1_name");
				}
				
				if (!-e $read2_symlink) {
					system("cd $read_directory/$sample_readdir; ln -s $read2_actual $read2_name");
				}
		
				#Check if data is compressed        
		        for my $suffix ( keys %suffices ) {
		            if ($lane_data->{$lane_count}{read1} =~ /$suffix$/) {
		            	$read1_compressed = 1;
		            	$read1_suffix = $suffix;
		            	#Add decompression if required
		            	if ($suffix eq 'bz2') {
							push @{$commands{bwa}},"$svndir/scripts/compress.pl -keep -suffix $suffix -threadNum $threadnum -files $read1_full_file";
						}
	
		            }
		            
		            if ($lane_data->{$lane_count}{read2} =~ /$suffix$/) {
		            	$read2_compressed = 1;
		            	$read2_suffix = $suffix;
		            	if ($suffix eq 'bz2') {
							push @{$commands{bwa}}, "$svndir/scripts/compress.pl -keep -suffix $suffix -threadNum $threadnum -files $read2_full_file";
						}
		            }
		        }
				
				#Create the aligning lane qsub files
	
				system("mkdir -p $lanes_outdir") if (!-d $lanes_outdir);
				  
		        my %lane_info = (
	        				lane_number=>$lane_number,
	        				lane_name=>$lane_name,
	        				sequencing_centre=>$sequencing_centre,
	        				sample_id=>$sample_db_id
	        				);
	        	$lane_db_id = modules::Adaptors::Lane->insert(\%lane_info);
	        	print STDERR "\t\t\t\tCreate lane $lane_name with id $lane_db_id\n";
	        	
	        	#Create the read_files
				my %read1_info = (
								file_name => $lane_data->{$lane_count}{read1},
								is_compressed=>$read1_compressed,
								compression_suffix=>$read1_suffix,
								read_file_number=>1,
								read_directory=>$read_directory_full,
								read_file_name=>$read1_name,
								lane_id=>$lane_db_id
								);
								
				my %read2_info = (
								file_name => $lane_data->{$lane_count}{read2},
								is_compressed=>$read2_compressed,
								compression_suffix=>$read2_suffix,
								read_file_number=>2,
								read_directory=>$read_directory_full,
								read_file_name=>$read2_name,
								lane_id=>$lane_db_id
								);        
			
				modules::Adaptors::Read_File->insert(\%read1_info);
				modules::Adaptors::Read_File->insert(\%read2_info);
			
				print STDERR "\t\t\t\t\tCreate two read_files $read1_name and $read2_name\n";
	        	
	        	
	        } else {
	        	my ($lane_obj) = modules::Adaptors::Lane->search(sample_id=>$sample_db_id,lane_name=>$lane_name);
				if (!$lane_obj) {
					modules::Exception->throw("ERROR: Can't find lane  $lane_name");
				}
	        	$lane_db_id = $lane_obj->id;
	        	
	        	my ($read_file_obj_r1) = modules::Adaptors::Read_File->search(lane_id=>$lane_db_id,read_file_number=>1);
				if (!$read_file_obj_r1) {
					modules::Exception->throw("ERROR: Can't find read file 1 for sample id $sample_db_id");
				}
				my ($read_file_obj_r2) = modules::Adaptors::Read_File->search(lane_id=>$lane_db_id,read_file_number=>2);
				if (!$read_file_obj_r2) {
					modules::Exception->throw("ERROR: Can't find read file 2 for sample id $sample_db_id");
				}	        	
	        	$read1_full_file = $read_directory .'/'. $sample_readdir. '/'. $read_file_obj_r1->file_name;
				$read2_full_file = $read_directory .'/'. $sample_readdir. '/'. $read_file_obj_r2->file_name;
	        }
	        
			if (keys %commands) {
				($encoding,$lane_bam) = $pipeline_sample->align_lane_qsub(-sample_name=>$sample_name,-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id,-commands=>\%commands);
			} else {
				($encoding,$lane_bam) = $pipeline_sample->align_lane_qsub(-sample_name=>$sample_name,-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id);
			}
			push @bams, $lane_bam;
			print "\n";
	    }
	    my $first_qsub = $sample_name .'*' . $first_lane_step . '.qsub';
	    $qsub_commands{"for f in $qsub_dir/lanes/$first_qsub; do qsub \$f; done"}++;
		system("mkdir -p $runs_outdir") if (!-d $runs_outdir);
	    

	   	#Now create the config xml for the runs
	    my $steps_xml_dir = $qsub_dir . '/runs/';
	    if (!-d $steps_xml_dir) {
	    	system("mkdir -p $steps_xml_dir");
	    }

	    my $steps_xml_out = $steps_xml_dir.$sample_name.'.xml';
	    #Assumes same encoding for a sample but likely ok
		$pipeline_sample->create_sample_xml(-sample_name=>$sample_name,-steps_xml_template=>$steps_xml,-xml_out=>$steps_xml_out,-encoding=>$encoding,-bams=>\@bams);

	   if ($parallel) {
                #Now create the run qsubs

                $pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$first_pipe_step,-end_step=>$last_align_step,-pipe_block=>1);

                #Generate the snv processing part of the pipeline run
                $pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$parvar_first,-end_step=>$parvar_first,-pipe_block=>2);
                $pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$parvar_merge,-end_step=>$last_pipe_step,-pipe_block=>3);

                if ($encoding eq 'phred33') {
					$pipeline_sample->create_parallel_snvcall_qsubs(-sample=>$sample_name,-source_group_name=>$source_group_name,-steps_xml=>$steps_xml_out);
				} else {
					$pipeline_sample->create_parallel_snvcall_qsubs(-sample=>$sample_name,-source_group_name=>$source_group_name,-steps_xml=>$steps_xml_out,-phred64=>1);
				}
	
 	    } else {

                #Now create the run qsubs; only a single block for each run here; not worth splitting 
                $pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$first_pipe_step,-end_step=>$last_pipe_step,-pipe_block=>1);

           }		
	    #$pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$first_pipe_step,-end_step=>$last_pipe_step,-pipe_block=>1);	  	    
	    
	}

}









