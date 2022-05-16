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
use modules::Adaptors::Human_Cancer_Sample;
use modules::Adaptors::Mouse_Cancer_Sample;
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
	   		"normal_external_names=s",
	   		"tumour_external_names=s",
	   		"normal_readdirs=s",
	   		"tumour_readdirs=s",
	   		"normal_sample_types=s",
	   		"tumour_sample_types=s",
	   		"sequencing_centre=s",
	   		"tumour_types=s",
	   		"tumour_tissues=s",
	   		"read1_pattern=s",
	   		"cluster_xml=s",
	   		"steps_xml=s",
	   		"skip_quality",
	   		"new_source_group",
	   		"new_source",
	   		"new_sample",
	   		"test",
	   		"only_qsubs",
	   		"source_group_number=i",
	   		"submit",
	   		"project_name=s",
	   		"mouse",
	   		"sequence_type=s",
	   		"no_mdss"
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{external_name} || !$OPT{normal_external_names} || !$OPT{tumour_external_names} || !$OPT{normal_readdirs} || !$OPT{tumour_readdirs} || !$OPT{normal_sample_types} || !$OPT{tumour_sample_types} || !$OPT{tumour_types});

	   
=pod

=head1 SYNOPSIS

add_cancer_sample.pl -normal_external_names external_name_normal -sequence_type sequence_experiment_type (exome,genome, or targeted; default=genome)-submit submit_jobs -new_source_group new_source_group_to_existing_patient(default=new_patient) -tumour_external_name external_name_tumour -external_name external_name(e.g.tb412) -read1_pattern read1_regex_pattern(default=R1) -steps_xml xml_steps_file(default=../conf/steps/xml) -cluster_xml xml_cluster_file(default=../conf/cluster.xml) -sequencing_center sequencing_centre(option=BRF,AGRF,RAM;default=AGRF) -normal_sample_types comma_delim_sample_type -tumour_sample_types comma_delim_sample_type -tumour_types comma_delim_tumour_types(options=primary,metatstatic)  -normal_readdirs comma_delim_fq_directories -tumour_readdirs comma_delim_fq_directories -skip_quality Skip_read_check_and_assume_phred33_qualities -only_qsubs only_create_qsubs -source_group_number source_group_number_for_only_qsubs_flag [options]

Required flags: -normal_external_names -tumour_external_names -external_name -normal_readdirs -tumour_readdirs -normal_sample_types -tumour_sample_types -tumour_types (-new_source || -new_source_group || -new_sample)

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

add_cancer_sample.pl -> Add patients, source_groups, and sample to db

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./add_cancer_sample.pl -tumour_readdirs MASCRI-434_27604_Primary,MASCRI-434_27604_Metastasis -tumour_types primary,metastatic -sequencing_centre AGRF -read1_pattern R1 -tumour_sample_types tumour_primary_patient,tumour_metastatic_patient -external_name MASCRI-434_27604 -normal_readdirs MASCRI-434_27604_Blood -normal_sample_types normal_patient -normal_external_names 102.100.100/7701 -tumour_external_names 102.100.100/7699,102.100.100/7700
./add_cancer_sample.pl -tumour_readdirs A15_metastatic_2,A15_primary -tumour_types metastatic,metastatic -sequencing_centre AGRF -read1_pattern R1 -tumour_sample_types tumour_metastatic_cell_line,tumour_metastatic_cell_line -external_name A15 -normal_readdirs A15_LCL_2 -normal_sample_types normal_cell_line -normal_external_names 102.100.100/7710 -tumour_external_names 102.100.100/7709,102.100.100/7709 -skip_quality -only_qsubs -source_group_number 2 -new_source_group -test
./add_cancer_sample.pl -tumour_tissues brain -tumour_readdirs tumour -tumour_types metastatic -sequencing_centre BRF -read1_pattern R1 -tumour_sample_types tumour_metastatic_patient -external_name SCC11 -normal_readdirs normal -normal_sample_types normal_patient -normal_external_names 102.100.100/7731 -tumour_external_names 102.102.100/7730 -skip_quality -only_qsubs -test -read1_pattern read1
./add_cancer_sample.pl -normal_external_names APFNS1890 -tumour_external_names APFNS1889 -external_name IGL01740_IGL01741 -normal_readdirs IGL01741 -tumour_readdirs IGL01740 -normal_sample_types mouse_cancer_normal -tumour_sample_types mouse_cancer_tumour -tumour_types metastatic -test -read1_pattern R1 -mouse -only_qsubs

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

if ($svndir =~ /trunk/) {
	modules::Exception->throw("ERROR: Do not use trunk for svn directory") unless $OPT{test};
}


#Check we have a single '-new' flag passed
if (!$OPT{new_sample} && !$OPT{new_source_group} && !$OPT{new_source} && !$OPT{only_qsubs}) {
	modules::Exception->throw("ERROR: Must pass either -new_sample, -new_source_group, or -new_source");
}

if ($OPT{new_sample} && $OPT{new_source_group}) {
	modules::Exception->throw("ERROR: Can only call single -new argument");
} elsif ($OPT{new_sample} && $OPT{new_source}) {
	modules::Exception->throw("ERROR: Can only call single -new argument");
} elsif ($OPT{new_source_group} && $OPT{new_source}) {
	modules::Exception->throw("ERROR: Can only call single -new argument");	
}

my $seq_type = defined $OPT{sequence_type}?$OPT{sequence_type}:'genome';

my $pipe_config = modules::Pipeline::get_pipe_conf();
my %possible_sample_types = map {$_ => 1} split(",",$pipe_config->read('common','sample_types'));

my $writedb = defined $OPT{only_qsubs}?0:1;

#Get the arguments
my $external_name = $OPT{external_name};
$external_name =~ s/ /_/g;
my $project_name = defined $OPT{project_name}?$OPT{project_name}:'Melanoma';
my $cancer_type = lc($project_name);
my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
if (!defined $project_obj) {
	modules::Exception->throw("ERROR: Can't find project name $project_name");
}

#Default is human_cancer
my $source_type = defined $OPT{mouse}?'mouse_cancer':'human_cancer';

#First create the patient entry

my $source_db_id;
my ($source_obj) = modules::Adaptors::Source->search('external_source_name'=>$external_name);

#Flag for whether it's new source, new sample group, or new sample
my $new_source = defined $OPT{new_source}?1:0;
my $new_source_group = defined $OPT{new_source_group}?1:0;
my $new_sample = defined $OPT{new_sample}?1:0;

if (defined $source_obj && $new_source) {
	modules::Exception->throw("ERROR: Trying to create new source but source exists with name $external_name");
}

if (!defined $source_obj && !$new_source) {
	modules::Exception->throw("ERROR: Trying to add to source $external_name but it doesn't exist");
}

my @normal_sample_types = split(",",$OPT{normal_sample_types});

if (@normal_sample_types > 1){
	modules::Exception->throw("ERROR: Can't have more than 1 normal sample");
}

my @tumour_sample_types = split(",",$OPT{tumour_sample_types});

#Check all the sample_types
for my $sample_type (@normal_sample_types, @tumour_sample_types) {
	if (!exists $possible_sample_types{$sample_type}) {
		modules::Exception->throw("ERROR: $sample_type is not a correct sample_type option");
	}    

	if ($sample_type !~ /cancer/) {
		modules::Exception->throw("ERROR: $sample_type is not a cancer sample_type");
	}
}

#Create the pipeline object
#Get the xml files to get variables we need first
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $steps_xml = defined $OPT{steps_xml}?$OPT{steps_xml}:"$svndir/conf/".$source_type.".steps.xml";

if ( !-e $steps_xml ) {
	modules::Exception->throw("File $steps_xml doesn't exist");	
}

my %sequencing_centres = map{ $_ => 1 } split(",",$pipe_config->read('common','sequencing_centres'));

my $sequencing_centre = defined $OPT{sequencing_centre}?$OPT{sequencing_centre}:'AGRF';
if (!exists $sequencing_centres{$sequencing_centre}) {
	my $seq_centre_str = join (',',keys %sequencing_centres);
	modules::Exception->throw("ERROR: Sequencing centre must be $seq_centre_str from pipe.xml");
}

my @normal_external_names = split(",",$OPT{normal_external_names});
my @tumour_external_names = split(",",$OPT{tumour_external_names});

my @normal_readdirs = split(",",$OPT{normal_readdirs});
my @tumour_readdirs = split(",",$OPT{tumour_readdirs});

my $normal_num = @normal_readdirs;
my $tumour_num = @tumour_readdirs;

my @tumour_types = split(",",$OPT{tumour_types});
my @tumour_tissues = ();

if (defined $OPT{tumour_tissues}) {
	@tumour_tissues = split(",",$OPT{tumour_tissues});
} else {
	#Default is skin
	@tumour_tissues = (('skin') x $tumour_num);
}

#Check all the argument numbers match
if (@tumour_external_names != @tumour_readdirs || @tumour_sample_types != @tumour_readdirs || @tumour_sample_types != @tumour_types || @tumour_sample_types != @tumour_tissues) {
	modules::Exception->throw("ERROR: Number of tumour arguments must match");
}

if (@normal_sample_types != @normal_readdirs) {
	modules::Exception->throw("ERROR: Number of normal arguments must match");
}




my $read_directory = $cluster_config->read($source_type,'base_directories','base_read_directory').'/'.$external_name;

my $read1_regex = defined $OPT{read1_pattern}?$OPT{read1_pattern}:'R1';

my $mdss = $cluster_config->read('common','mdss','mdss_flag');
my $sys_call = modules::SystemCall->new();

if ($mdss && $writedb) {
	#Create a backup read files in this case
	my $mdss_read_dir = $cluster_config->read('common','mdss','mdss_reads');
	
	#marcin:
	#my $mdss_command = "mdss put -r $read_directory $mdss_read_dir";	
	my $mdss_command = "$svndir/scripts/mdss_put.sh $read_directory $mdss_read_dir";
	$sys_call->run($mdss_command) unless $OPT{no_mdss};
}

#Check reads match if alread
if ( !-d $read_directory ) {
	modules::Exception->throw("File $read_directory doesn't exist");	
}

for my $readdir (@tumour_readdirs, @normal_readdirs) {
	my $output = `ls $read_directory/$readdir/*$read1_regex* 2>/dev/null`;
	if (!$output) {
		modules::Exception->throw("ERROR: No read files match  $read_directory/$readdir/*$read1_regex*");
	}
	
}


if ($writedb) {
	#normal case; writing out new objectss
	if ($new_source) {
		#normal case here; insert objects into database
		my %source_info = (
							external_source_name=>$external_name,
							project_id=>$project_obj->id,
							source_type=>$source_type
				   );
		$source_db_id = modules::Adaptors::Source->insert(\%source_info);
		($source_obj) = modules::Adaptors::Source->search('id'=>$source_db_id);
		print STDERR "Create source $external_name with id $source_db_id\n";
	} else {
		$source_db_id = $source_obj->id;
		print STDERR "Add to source $external_name with id $source_db_id\n";
	}
} else {
	$source_db_id = $source_obj->id;
}

if (!defined $source_obj) {
	modules::Exception->throw("ERROR: Can't retrieve source object for $external_name");
}

my $source_group_count;
if ($writedb) {
	if ($new_source) {
		#First source group
		$source_group_count = 1;
	} elsif ($new_source_group) {
		#Get the current highest source group and add one
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
	} elsif (defined $OPT{source_group_number}) {
		#updating qsubs we need the source_group_number if more than one source_group
		$source_group_count = $OPT{source_group_number}; 
	} else {
		#Just adding new sample; without source_group_number confirm there's only one source group that exists
		#Writing out qsubs; confirm there aren't multiple source_groups as the source_group_number would be needed
		my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
		if (@source_group_objs > 1) {
			modules::Exception->throw("ERROR: Multiple source groups for source $external_name");
		}
		$source_group_count = 1;
	}
} else {
	#Updating qsubs we need the source_group_number if more than one source_group
 	if (defined $OPT{source_group_number}) {
		$source_group_count = $OPT{source_group_number};
	} else {
		my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
		if (@source_group_objs > 1) {
			modules::Exception->throw("ERROR: Multiple source groups for source $external_name");
		}
		$source_group_count = 1;
	}
}

#Now create the database entry or get the id if it already exists
my $source_group_name = $external_name.'_sg'.$source_group_count;
my $sourcegroup_db_id;
my $sample_total;

if ($writedb) {
	if ($new_source || $new_source_group) {
		#Create new source_group
		#Next create the sample group
		$sample_total = @normal_sample_types + @tumour_sample_types;
		my %source_group_info = (
							total_samples=>$sample_total,
							source_group_number=>$source_group_count,
							source_group_name=>$source_group_name,
							source_id=>$source_db_id
							);
		$sourcegroup_db_id = modules::Adaptors::Source_Group->insert(\%source_group_info);
		print STDERR "\tCreate source_group $source_group_name with id $sourcegroup_db_id\n";		
	}  elsif ($new_sample) {
		my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$source_group_name,source_id=>$source_db_id);
		
		#Don't create source_group but update total_samples variable
		my $current_samples = $source_group_obj->total_samples;
		$sample_total = $current_samples + 1; #Account for sample we're adding
		$source_group_obj->total_samples($sample_total);
		$source_group_obj->update();
		$sourcegroup_db_id = $source_group_obj->id;
	}
	
} else {
	my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$source_group_name,source_id=>$source_db_id);
	if (!defined $source_group_obj) {
		modules::Exception->throw("No source_group object for $source_group_name and source id $source_db_id");
	}
	$sourcegroup_db_id = $source_group_obj->id;
}


my $first_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first'); 
my $last_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'last');
my $last_align_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'finish_align');
my $parvar_first = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_first');
my $parvar_merge = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_merge');
my $nonpar_start = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'start_non_parvar');
my $first_lane_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first_lane');

#Read the cluster specific variables 
my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory').'/'.$external_name;
system("mkdir $qsub_dir") if (!-d $qsub_dir);

my $outdir = $cluster_config->read($source_type,'base_directories','base_run_directory').'/'.$external_name;
system("mkdir -p $outdir") if (!-d $outdir);
my $snv_call_dir = $outdir .'/' .$source_group_name.'_snvcalls';
system("mkdir -p $snv_call_dir") if (!-d $snv_call_dir);

my $bam_link_dir = $outdir . '/bam_links';
system("mkdir -p $bam_link_dir") if (!-d $bam_link_dir);



my $threadnum = $cluster_config->read('common','qsub_vars','thread_num');
my $scheduler = $cluster_config->read('common','scheduler');
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
my @normal_samples  = ();
push @normal_samples, $source_group_name.'_normal1';

my @tumour_samples = ();




if ($new_source || $new_source_group || !$writedb) {
	for ( my $sample_count = 0 ; $sample_count < $tumour_num ; $sample_count++ ) {
		my $sample_number = $sample_count+1;
		my $sample_name = $source_group_name.'_tumour'.$sample_number;
		push @tumour_samples,$sample_name;
	}
	#First normal
	&Create_Samples('normal',\@normal_samples);
	#Next tumour
	&Create_Samples('tumour',\@tumour_samples);
} elsif ($new_sample) {
	#Just create the new sample
	my $new_sample_count = 0;
	#Get the next tumour count
	my @sample_objs = modules::Adaptors::Sample->search(source_group_id=>$sourcegroup_db_id);
	for my $sample_obj (@sample_objs) {
		if ($sample_obj->sample_name =~ /tumour(\d)/) {
			if ($1 > $new_sample_count) {
				$new_sample_count = $1;
			}
		}
	}
	$new_sample_count++;
	push @tumour_samples, $source_group_name.'_tumour'.$new_sample_count;
	&Create_Samples('tumour',\@tumour_samples);
}


#Create the snv_calling qsubs # germline calls are now obtained as subset of pair calls
$pipeline_sample->create_parallel_cancer_snvcall_qsubs(-normal_sample=>$normal_samples[0],-tumour_samples=>\@tumour_samples,-source_group_name=>$source_group_name,-call_normal=>0);

if (keys %qsub_commands) {
	print STDERR "Run the following commands...\n";
	for my $command ( keys %qsub_commands ) {
	    print "\t$command\n";
	    system("$command") if $OPT{submit};
	}
	print "\n";
}

#Does all the work creating samples, lanes, and read_files
sub Create_Samples () {
	my ($tumour_normal,$sample_names) = @_;	
	my $tumour_flag = $tumour_normal eq 'tumour'?1:0;


	for ( my $sample_count = 0 ; $sample_count < @{$sample_names} ; $sample_count++ ) {
	    my $cell_line = 0;
	    if ($tumour_flag) {
	    	if ($tumour_sample_types[$sample_count] =~ /cell_line/) {
		    	$cell_line = 1;
		    }
	    } else {
		    if ($normal_sample_types[$sample_count] =~ /cell_line/) {
		    	$cell_line = 1;
		    }
	    }
		my $sample_number = $sample_count+1;
	
		my $sample_name = $sample_names->[$sample_count];
	
		#Figure out how many lanes there are
		my $lane_data;
		
		#Get the local sample_directory for getting reads later
			
		my $read_directory_full;
		my $sample_readdir;
		
		if ($tumour_flag) {
			$sample_readdir = $tumour_readdirs[$sample_count];
			$read_directory_full = $read_directory .'/'. $sample_readdir;
			($lane_data) = modules::Pipeline::get_lane_info($read_directory_full,$read1_regex);
		} else {
			$sample_readdir = $normal_readdirs[$sample_count];
			$read_directory_full = $read_directory .'/'. $sample_readdir;
			
			($lane_data) = modules::Pipeline::get_lane_info($read_directory_full,$read1_regex);
		}
	    my $lane_total = keys %{$lane_data};
	    	    
	    	    
	    
	    my $sample_db_id;
	    my $cancer_sample_db_id;
	    if ($writedb) {
		    my %sample_info = (
		    					sample_number=>$sample_number,
		    					total_lanes => $lane_total,
		    					sample_name => $sample_name,
		    					source_group_id=>$sourcegroup_db_id,
		    					sequence_type=>$seq_type
		    					);
		    
		    
		    if ($tumour_flag) {
		    	$sample_info{external_sample_name} = $tumour_external_names[$sample_count];	    						    	
		    	$sample_info{sample_type} = $tumour_sample_types[$sample_count];
		    } else {
		    	$sample_info{external_sample_name} = $normal_external_names[$sample_count];
		    	$sample_info{sample_type} = $normal_sample_types[$sample_count];
		    }
		    
	    	$sample_db_id =  modules::Adaptors::Sample->insert(\%sample_info);
		    print STDERR "\t\tCreate sample $sample_name with id $sample_db_id\n";   
		    
		    my %cancer_sample_info = (
		    						tumour => $tumour_flag,
		    						cell_line => $cell_line,
		    						cancer_type => $cancer_type,
		    						sample_id=>$sample_db_id
		    						);
		    
		    if ($tumour_flag) {
				$cancer_sample_info{tissue}= $tumour_tissues[$sample_count];
				$cancer_sample_info{tumour_type} = $tumour_types[$sample_count];
		    } else {
		    	$cancer_sample_info{tissue}='blood';
		    }
		    
		    if ($source_type eq 'mouse_cancer') {
		    	$cancer_sample_db_id =  modules::Adaptors::Mouse_Cancer_Sample->insert(\%cancer_sample_info);
		    } else {
			    $cancer_sample_db_id =  modules::Adaptors::Human_Cancer_Sample->insert(\%cancer_sample_info);
		    } 	
		    print STDERR "\t\t\tCreate cancer_sample $sample_name with id $cancer_sample_db_id\n";   
		     	
	    } else {
	    	#Get the existing sample
	    	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name,source_group_id=>$sourcegroup_db_id);
	    	print "Sample $sample_name source_groupdb $sourcegroup_db_id\n";
	    	$sample_db_id = $sample_obj->id;
	    }				
	    
	    
  
	    
	    my @bams = ();			
	    my $lanes_outdir = $outdir.'/'.$sample_name.'_lanes';
	    my $runs_outdir = $outdir.'/'.$sample_name.'_runs';
	    
	    my $encoding;
		#Create the lanes for each sample
	    for ( my $lane_count = 0 ; $lane_count < $lane_total ; $lane_count++ ) {
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
			
			
			my $read1_actual = $lane_data->{$lane_count}{read1};
			my $read2_actual = $lane_data->{$lane_count}{read2};
			
			my $read1_full_file = $read_directory .'/'. $sample_readdir. '/'. $lane_data->{$lane_count}{read1};
			my $read2_full_file = $read_directory .'/' . $sample_readdir. '/'. $lane_data->{$lane_count}{read2};
			my $read1_symlink = $read_directory .'/' . $sample_readdir. '/'. $read1_name;
			my $read2_symlink = $read_directory .'/' . $sample_readdir. '/'. $read2_name;
			
			
			
			
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
	
			my %commands = ();
	
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
			my $lane_bam;
			system("mkdir -p $lanes_outdir") if (!-d $lanes_outdir);
			
			#By default don't check the quality (assume phred33)
			my $quality = defined $OPT{skip_quality}?1:0;
			
	        
	        
	        my $lane_db_id;
	        if ($writedb) {
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
	        	$lane_db_id = $lane_obj->id;
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
	    $qsub_commands{"for f in $first_qsub; do qsub \$f; done"}++;
		system("mkdir -p $runs_outdir") if (!-d $runs_outdir);
	    
	   	#Now create the config xml for the runs
	    my $steps_xml_dir = $qsub_dir . '/runs/';
	    if (!-d $steps_xml_dir) {
	    	system("mkdir -p $steps_xml_dir");
	    }
	    
	    my $steps_xml_out = $steps_xml_dir.$sample_name.'.xml';
	    #Assumes same encoding for a sample but likely ok
		$pipeline_sample->create_sample_xml(-sample_name=>$sample_name,-steps_xml_template=>$steps_xml,-xml_out=>$steps_xml_out,-encoding=>$encoding,-bams=>\@bams);

	    #Now create the run qsubs
		$pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$first_pipe_step,-end_step=>$last_align_step,-pipe_block=>1);	  	    
	    
	    #Generate the snv processing part of the pipeline run
	    if ($tumour_flag == 1) {
	    	if ($source_type eq 'mouse_cancer') {
				$pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$nonpar_start,-end_step=>$last_pipe_step,-pipe_block=>2);								    		
	    	} else {
			    $pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$parvar_first,-end_step=>$parvar_first,-pipe_block=>2);
			    $pipeline_sample->create_run_qsubs(-sample_name=>$sample_name,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$parvar_merge,-end_step=>$last_pipe_step,-pipe_block=>3);
	    	}
	    } 
	}
}









