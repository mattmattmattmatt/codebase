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
use modules::Adaptors::Human_Single_Sample;
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
	   		"readdir=s",
	   		"sequencing_centre=s",
	   		"read1_pattern=s",
	   		"cluster_xml=s",
	   		"steps_xml=s",
	   		"skip_quality",
	   		"test",
	   		"only_qsubs",
	   		"project_name=s",
	   		"unaffected",
	   		"submit",
	   		"sequence_type=s",
	   		"external_sample_name=s",
	   		"no_mdss",
	   		"no_parallel",
	   		"gatk",
	   		"apf_req_id=s",
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{external_name});

	   
=pod

=head1 SYNOPSIS

add_single_human_sample.pl -no_parallel run_in_serial(default=parallel) -sequence_type sequence_experiment_type (exome,genome, or targeted; default=exome) -project_name db_project_name(default='random human') -submit submit_jobs -external_name external_name(e.g.tb412) -external_sample_name external_sample_name(default=external_name) -unaffected sample_is_unaffected(default=affected) -read1_pattern read1_regex_pattern(default=R1) -steps_xml xml_steps_file(default=../conf/steps/xml) -cluster_xml xml_cluster_file(default=../conf/cluster.xml) -sequencing_center sequencing_centre(option=BRF,AGRF,RAM;default=AGRF) -sample_type sample_type(default=base_on_affected_flag) -readdir fq_directories(default=source_type_xml/external_name) -skip_quality Skip_read_check_and_assume_phred33_qualities -only_qsubs only_create_qsubs [options]

Required flags: -external_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

add_single_human_sample.pl -> Add single human patient not part of any cohort

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./add_single_human_sample.pl 

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

if ($svndir =~ /trunk/) {
	modules::Exception->throw("ERROR: Do not use trunk for svn directory") unless $OPT{test};
}

my $seq_type = defined $OPT{sequence_type}?$OPT{sequence_type}:'exome';
my $use_gatk = defined $OPT{gatk}?1:0;
my $source_type = 'human_single';
$source_type .= '_gatk' if $use_gatk;


#Determine whether it's affected or not and set sample_type
my $affected = defined $OPT{unaffected}?0:1;
my $sample_type;

if ($affected) {
	$sample_type = 'single_affected';
} else {
	$sample_type = 'single_unaffected';
}

my $parallel = 1;

my $report_conf = modules::Pipeline::get_report_conf();


my $writedb = defined $OPT{only_qsubs}?0:1;

#Get the arguments
my $external_name = $OPT{external_name};
$external_name =~ s/ /_/g;

my $external_sample_name = defined $OPT{external_sample_name}?$OPT{external_sample_name}:$external_name;

my $project_name = defined $OPT{project_name}?$OPT{project_name}:'Random human';
my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
if (!defined $project_obj) {
	modules::Exception->throw("ERROR: Can't find project name $project_name; must create in the db first if not already done");
}

#Create the pipeline object
#Get the xml files to get variables we need first
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $steps_xml = defined $OPT{steps_xml}?$OPT{steps_xml}:"$svndir/conf/".$source_type.".steps.xml";

if ( !-e $steps_xml ) {
	modules::Exception->throw("File $steps_xml doesn't exist");	
}

my $pipe_config = modules::Pipeline::get_pipe_conf();
my %sequencing_centres = map{ $_ => 1 } split(",",$pipe_config->read('common','sequencing_centres'));

my $sequencing_centre = defined $OPT{sequencing_centre}?$OPT{sequencing_centre}:'AGRF';
if (!exists $sequencing_centres{$sequencing_centre}) {
	my $seq_centre_str = join (',',keys %sequencing_centres);
	modules::Exception->throw("ERROR: Sequencing centre must be $seq_centre_str from pipe.xml");
}

my $read_base = $cluster_config->read($source_type,'base_directories','base_read_directory');
my $readdir = defined $OPT{readdir}?$OPT{readdir}:$external_name;
my $read_directory = "$read_base/$readdir";

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


if (!-d $read_directory) {
	modules::Exception->throw("ERROR: Can't open readdir $read_directory");
}

my $read1_regex = defined $OPT{read1_pattern}?$OPT{read1_pattern}:'R1';

my $output = `ls $read_directory/*$read1_regex* 2>/dev/null`;
if (!$output && $writedb) {
	modules::Exception->throw("ERROR: No read files match $read_directory/*$read1_regex*");
}

#First create the patient entry
my %source_info = (
					external_source_name=>$external_name,
					project_id=>$project_obj->id,
					source_type=>$source_type
				   );

my $source_db_id;
my ($source_obj)  = modules::Adaptors::Source->search('external_source_name'=>$external_name);

if ($writedb) {
	#normal case; writing out new objects
	if (defined $source_obj) {
		modules::Exception->throw("ERROR: A source already exists with this external_name ($external_name). Must delete source first");
	}
	#normal case here; insert objects into database
	$source_db_id = modules::Adaptors::Source->insert(\%source_info);
	($source_obj) = modules::Adaptors::Source->search('id'=>$source_db_id);
	print STDERR "Create source $external_name with id $source_db_id\n";
} else {
	($source_obj) = modules::Adaptors::Source->search('external_source_name'=>$external_name);
	$source_db_id = $source_obj->id;
}

if (!defined $source_obj) {
	modules::Exception->throw("ERROR: Can't retrieve source object for $external_name");
}

my $source_group_count;
if ($writedb) {
	$source_group_count = 1;
} else {
	#Writing out qsubs in patient; confirm there aren't multiple source_groups
	my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
	if (@source_group_objs > 1) {
		modules::Exception->throw("ERROR: Multiple source groups for source $external_name");
	}
	$source_group_count = 1;
}

my $first_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first'); 
my $last_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'last');
my $last_align_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'finish_align');
my $parvar_first = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_first');
my $parvar_merge = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_merge');
my $nonpar_start = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'start_non_parvar');
my $first_lane_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first_lane');

#Next create the sample group
my $sample_total = 1;

#This is the first source_group
my $source_group_name = $external_name.'_sg'.$source_group_count;


my $sourcegroup_db_id;
if ($writedb) {
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


# Make directories for output, for SNV calls and BAM file links
my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory').'/'.$external_name;
system("mkdir -p $qsub_dir") if (!-d $qsub_dir);

my $outdir = $cluster_config->read($source_type,'base_directories','base_run_directory').'/'.$external_name;
system("mkdir -p $outdir") if (!-d $outdir);

my $snv_call_dir = $outdir .'/' .$source_group_name.'_snvcalls';
system("mkdir -p $snv_call_dir") if (!-d $snv_call_dir);

my $bam_link_dir = $outdir . '/bam_links';
system("mkdir -p $bam_link_dir") if (!-d $bam_link_dir);

# Read the cluster specific variables
my $threadnum = $cluster_config->read('common','qsub_vars','thread_num');
my $scheduler = $cluster_config->read('common','scheduler');
my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
my $pipeline_sample = modules::PipelineSample->new(-source_name=>$external_name,-cluster_obj=>$cluster_obj,-source_type=>$source_type);

#Acceptable suffices
my %suffices = (
    gz=>1,
    bz2=>1
    );

my %qsub_commands = ();

# Create samples
# - builds all the specific qsub files for running this sample

&Create_Single_Sample(); # MONSTER subroutine that does _everything_.  With global scope.  Many things changed by this call.

my $sample_name = $source_group_name.'_humansingle1';

if ($parallel) {
    # Really important, as a separate set of qsubs made for parallel SNV calling
	#$pipeline_sample->create_parallel_snvcall_qsubs(-sample            => $sample_name,
#							-source_group_name => $source_group_name);
}

if (keys %qsub_commands) {
	print STDERR "Run the following commands...\n";
	for my $command ( keys %qsub_commands ) {
	    print "\t$command\n";
	    system("$command") if $OPT{submit};
	}
	print "\n";
}

#Does all the work creating samples, lanes, and read_files

sub Create_Single_Sample {
    my $sample_number = 1; # this code has heritage in the add_cohort_sample.pl, which does this subroutine step several times.
	
    my $sample_name = $source_group_name.'_humansingle'.$sample_number; # GLOBAL: $source_group_name

    #Get the local sample_directory for getting reads later
    my $lane_data; 	

    #Figure out how many lanes there are
    my $lane_total;# = keys %{$lane_data};
    
    my $sample_db_id;
    my $human_single_sample_db_id;

    # Choose what to do depending if we are running in qsub only mode, or are attached to db
    if ($writedb) { # GLOBAL : $writedb

	($lane_data) = modules::Pipeline::get_lane_info($read_directory,$read1_regex); # GLOBAL: $read_directory, $read1_regex
	
	$lane_total = keys %{$lane_data};

	my %sample_info = (
	    sample_number        => $sample_number,
	    total_lanes          => $lane_total,
	    sample_name          => $sample_name,
	    source_group_id      => $sourcegroup_db_id, # GLOBAL: $sourcegroup_db_id
	    external_sample_name => $external_sample_name, # GLOBAL : $external_sample_name
	    sample_type          => $sample_type, # GLOBAL: $sample_type
	    sequence_type        => $seq_type, # GLOBAL : $seq_type
	    );
    #Optional apf request id field, applies to samples from APF only
	if(defined $OPT{apf_req_id}){
		$sample_info{apf_request_id} = $OPT{apf_req_id};
	}    
	$sample_db_id =  modules::Adaptors::Sample->insert(\%sample_info);
	print STDERR "\t\tCreate sample $sample_name with id $sample_db_id\n";   
	    
	my %human_single_sample_info = (
	    affected  => $affected, # GLOBAL : $affected
	    sample_id => $sample_db_id,
	    );
	    	
	# -- add this sample to human_single_samples table
	$human_single_sample_db_id =  modules::Adaptors::Human_Single_Sample->insert(\%human_single_sample_info);
	print STDERR "\t\t\tCreate single_sample $sample_name with id $human_single_sample_db_id\n";   
		     	
    } else {
    	# or if this sample already exists, get the existing sample
    	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name,source_group_id=>$sourcegroup_db_id);
    	print "Sample $sample_name source_groupdb $sourcegroup_db_id\n";
    	$sample_db_id = $sample_obj->id;
	$lane_total = $sample_obj->total_lanes;
    }				
    
    my @bams = ();			
    my $lanes_outdir = $outdir.'/'.$sample_name.'_lanes';
    my $runs_outdir = $outdir.'/'.$sample_name.'_runs';
    
    my $encoding;

    # Create the lanes for each sample 
    # - this includes counting the number of input read pair files to determine the number of lanes, and for each of these creating the qsub files for alignment execution
    for ( my $lane_count = 0 ; $lane_count < $lane_total ; $lane_count++ ) {

	my $lane_db_id;
	my $lane_bam;
    	my %commands = ();
	
	my $lane_number = $lane_count + 1;
    	my $lane_name = $sample_name.'_l'.$lane_number;
        
        #Get the read_file info - get actual file names for each lane, the full path and generate a symlink name for each
        my $read1_compressed = 0;
        my $read2_compressed = 0;
        my $read1_suffix;
        my $read2_suffix;
	my $read1_name = $lane_name.'_r1';
	my $read2_name = $lane_name.'_r2';
			
	my $read1_full_file;
	my $read2_full_file;

	# -- by default don't check the quality (assume phred33my $quality = defined $OPT{skip_quality}?1:0;
	my $quality = defined $OPT{skip_quality}?1:0;
	                
	if($writedb) {	

		my $read1_actual = $lane_data->{$lane_count}{read1};
		my $read2_actual = $lane_data->{$lane_count}{read2};
		$read1_full_file = $read_directory . '/'. $lane_data->{$lane_count}{read1};
		$read2_full_file = $read_directory . '/'. $lane_data->{$lane_count}{read2};
		my $read1_symlink = $read_directory . '/'. $read1_name;
		my $read2_symlink = $read_directory . '/'. $read2_name;
				
		#First check the read files exist
		if ( !-e $read1_full_file ) {
		    modules::Exception->throw("File $read1_full_file doesn't exist");	
		}	
		if ( !-e $read2_full_file ) {
		    modules::Exception->throw("File $read2_full_file doesn't exist");	
		}	
		
		#Create symlinks to the files if they don't exist
		if (!-e $read1_symlink) {
		    system("cd $read_directory; ln -s $read1_actual $read1_name");
		}
			
		if (!-e $read2_symlink) {
		    system("cd $read_directory; ln -s $read2_actual $read2_name");
		}
	
	
		#Check if data is compressed and add decompression commands if required       
	        for my $suffix ( keys %suffices ) {
	
	            if ($lane_data->{$lane_count}{read1} =~ /$suffix$/) {
	            	$read1_compressed = 1;
	            	$read1_suffix = $suffix;
	
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
			
			       
		# -- insert information for each lane into db, it this hasn't already been done
        
	    	my %lane_info = (
			lane_number       => $lane_number,
			lane_name         => $lane_name,
			sequencing_centre => $sequencing_centre,
			sample_id         => $sample_db_id,
			);
	    
	   	 $lane_db_id = modules::Adaptors::Lane->insert(\%lane_info);
	   	 print STDERR "\t\t\t\tCreate lane $lane_name with id $lane_db_id\n";
        		
		 #Create the information for inserting entries into the read_files table in the database
		 my %read1_info = (
			file_name          => $lane_data->{$lane_count}{read1},
			is_compressed      => $read1_compressed,
			compression_suffix => $read1_suffix,
			read_file_number   => 1,
			read_directory     => $read_directory,
			read_file_name     => $read1_name,
			lane_id            => $lane_db_id,
		 );
		
	    	my %read2_info = (
			file_name          => $lane_data->{$lane_count}{read2},
			is_compressed      => $read2_compressed,
			compression_suffix => $read2_suffix,
			read_file_number   => 2,
			read_directory     => $read_directory,
			read_file_name     => $read2_name,
			lane_id            => $lane_db_id,
		);        
		
	    	modules::Adaptors::Read_File->insert(\%read1_info);
	    	modules::Adaptors::Read_File->insert(\%read2_info);
		
	    	print STDERR "\t\t\t\t\tCreate two read_files $read1_name and $read2_name\n";
        		
        } else { # if already in db, re-derive this information from the db
	    	my ($lane_obj) = modules::Adaptors::Lane->search(sample_id => $sample_db_id,lane_name => $lane_name);
		if (!$lane_obj) {
			modules::Exception->throw("ERROR: Can't find land $lane_name");
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

		$read1_full_file = $read_directory .'/'. $read_file_obj_r1->file_name;
		$read2_full_file = $read_directory .'/'. $read_file_obj_r2->file_name;
	}
        
	my %align_lane_qsub_args 
	    = (-sample_name  => $sample_name,
	       -skip_quality => $quality, 
	       -outdir       => $lanes_outdir,
	       -read1        => $read1_full_file,
	       -read2        => $read2_full_file,
	       -lane_name    => $lane_name,
	       -lane_id      => $lane_db_id,
	    );

	$align_lane_qsub_args{-commands} = \%commands 
	    if (keys %commands); # presently, these commands pertain to decompressing reads, and need to be appended here if this is needed

	($encoding,
	 $lane_bam) = $pipeline_sample->align_lane_qsub(%align_lane_qsub_args);

	push @bams, $lane_bam;
	print STDERR "\n";
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
    $pipeline_sample->create_sample_xml(-sample_name        => $sample_name,
					-steps_xml_template => $steps_xml,
					-xml_out            => $steps_xml_out,
					-encoding           => $encoding,
					-bams               => \@bams);

    if ($parallel) {
	if ($encoding eq 'phred33') {
		$pipeline_sample->create_parallel_snvcall_qsubs(-sample            => $sample_name,
                                                        -source_group_name => $source_group_name,
                                                        -steps_xml   => $steps_xml_out);	
	} else {
		$pipeline_sample->create_parallel_snvcall_qsubs(-sample            => $sample_name,
                                                        -source_group_name => $source_group_name,
														-steps_xml   => $steps_xml_out,
														-phred64=>1);     
	}
	#Now create the run qsubs - pipe_block 1 - ending after align step
	$pipeline_sample->create_run_qsubs(-sample_name => $sample_name,
					   -qsub_dir    => $steps_xml_dir,
					   -steps_xml   => $steps_xml_out,
					   -start_step  => $first_pipe_step,
					   -end_step    => $last_align_step,
					   -pipe_block  => 1);
	  	    
	#Generate the snv processing part of the pipeline run - pipe_block 2 - starting at first parallel step and ending at this same step
	$pipeline_sample->create_run_qsubs(-sample_name => $sample_name,
					   -qsub_dir    => $steps_xml_dir,
					   -steps_xml   => $steps_xml_out,
					   -start_step  => $parvar_first,
					   -end_step    => $parvar_first,
					   -pipe_block  => 2);

	# And get the rest done - pipe_block 3 - starting at the merging of the results of the parallel steps and going straight on to the end
	$pipeline_sample->create_run_qsubs(-sample_name => $sample_name,
					   -qsub_dir    => $steps_xml_dir,
					   -steps_xml   => $steps_xml_out,
					   -start_step  => $parvar_merge,
					   -end_step    => $last_pipe_step,
					   -pipe_block  => 3);		
    } else {

	#Now create the run qsubs; only a single block for each run here; not worth splitting 
	$pipeline_sample->create_run_qsubs(-sample_name => $sample_name,
					   -qsub_dir    => $steps_xml_dir,
					   -steps_xml   => $steps_xml_out,
					   -start_step  => $first_pipe_step,
					   -end_step    => $last_pipe_step,
					   -pipe_block  => 1);	  	    
    }

} #end subroutine Create_Single_Sample











