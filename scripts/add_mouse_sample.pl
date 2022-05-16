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
use modules::Adaptors::Mouse_Single_Sample;
use modules::Adaptors::Lane;
use modules::Adaptors::Read_File;
use modules::Adaptors::Run;
use Pod::Usage;
use Cwd;

use vars qw(%OPT);


GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"sample_name=s",
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
	   		"new_source",
	   		"new_source_group",
	   		"new_sample",
	   		"source_group_number=i",
	   		"sample_type=s",
	   		"strain=s",
	   		"submit",
	   		"sequence_type=s",
	   		"no_mdss"
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{external_name} || !$OPT{sample_name});

	   
=pod

=head1 SYNOPSIS

add_mouse_sample.pl -sample_type sample_type(default='G1_internal') -sequence_type sequence_experiment_type (exome,genome, or targeted; default=exome) -sample_name external_sample_name -submit submit_commands -strain strain(default=B6) -project_name db_project_name(default='ENU mouse') -external_name external_name(e.g.BRF_mouse1) -new_sample new_sample_to_existing_source_group -new_source_group new_source_group_to_existing_patient(default=new mouse) -sample_group_number sample_group_number_when_called -read1_pattern read1_regex_pattern(default=R1) -steps_xml xml_steps_file(default=../conf/steps/xml) -cluster_xml xml_cluster_file(default=../conf/cluster.xml) -sequencing_center sequencing_centre(option=BRF,AGRF,RAM;default=AGRF) -sample_type -new_source new_source(no_bioreps_in_db) -readdir fq_directories(default=source_type_xml/external_name) -skip_quality Skip_read_check_and_assume_phred33_qualities -only_qsubs only_create_qsubs [options]

Required flags: -external_name -sample_name (-new_source || -new_source_group || -new_sample || -only_qsubs)

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

add_mouse_sample.pl -> Add single mouse as new source, add a new source_group, or add a new source to an existing source group

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./add_mouse_sample.pl -new_source -external_name IGL00506_group -sample_name IGL00506 
./add_mouse_sample.pl -new_sample -external_name IGL00506_group -sample_name IGL00507 

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

my $seq_type = defined $OPT{sequence_type}?$OPT{sequence_type}:'exome';

my $source_type = 'mouse_single';

my $sample_type = defined $OPT{sample_type}?$OPT{sample_type}:'G1_internal';

my $report_conf = modules::Pipeline::get_report_conf();

my $writedb = defined $OPT{only_qsubs}?0:1;

my $strain = defined $OPT{strain}?$OPT{strain}:'B6';

#Get the arguments
my $external_name = $OPT{external_name};

$external_name =~ s/ /_/g;


my $sample_name = $OPT{sample_name};

my $project_name = defined $OPT{project_name}?$OPT{project_name}:'ENU mouse';
my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
if (!defined $project_obj) {
	modules::Exception->throw("ERROR: Can't find project name $project_name; must create in the db first if not already done");
}

my $source_db_id;
my ($source_obj)  = modules::Adaptors::Source->search('external_source_name'=>$external_name);

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
my $readdir = defined $OPT{readdir}?$OPT{readdir}:$sample_name;
my $read_directory = "$read_base/$readdir";

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

if (!-d $read_directory) {
	modules::Exception->throw("ERROR: Can't open readdir $read_directory");
}
my $read1_regex = defined $OPT{read1_pattern}?$OPT{read1_pattern}:'R1';

my $output = `ls $read_directory/*$read1_regex* 2>/dev/null`;
if (!$output) {
	modules::Exception->throw("ERROR: No read files match $read_directory/*$read1_regex*");
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

#First get the source_group_count
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
		#Create new source_group and single sample in either case
		$sample_total = 1;
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
		$sample_total = $current_samples + 1;
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
my $first_lane_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'first_lane');

#Read the cluster specific variables 
my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory').'/'.$external_name;
system("mkdir $qsub_dir") if (!-d $qsub_dir);
my $outdir = $cluster_config->read($source_type,'base_directories','base_run_directory').'/'.$external_name;
system("mkdir -p $outdir") if (!-d $outdir);

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
if ($writedb) {
	&Create_Single_Sample($sample_total);
} else {
	#Just update all the sample qsubs
	&Update_Qsubs();
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
	my ($sample_number) = @_;
	
	my $sample_name_local = $source_group_name.'_mousesingle'.$sample_number;

	#Get the local sample_directory for getting reads later
	my ($lane_data) = modules::Pipeline::get_lane_info($read_directory,$read1_regex);	

	#Figure out how many lanes there are
    my $lane_total = keys %{$lane_data};
    
    my $sample_db_id;
    my $mouse_single_sample_db_id;
    
	if ($writedb) {    	    
	    my %sample_info = (
	    					sample_number=>$sample_number,
	    					total_lanes => $lane_total,
	    					sample_name => $sample_name_local,
	    					source_group_id=>$sourcegroup_db_id,
	    					external_sample_name=>$sample_name,
	    					sample_type=>$sample_type,
	    					sequence_type=>$seq_type
	    					);
	    
	    $sample_db_id =  modules::Adaptors::Sample->insert(\%sample_info);
	    print STDERR "\t\tCreate sample $sample_name_local with id $sample_db_id\n";   
	    my $internal = $sample_type =~ /internal/?1:0;
	    my %mouse_single_sample_info = (
	    								sample_id=>$sample_db_id,
	    								internal=>$internal,
	    								strain=>$strain
	    								);
	    								
		$mouse_single_sample_db_id =  modules::Adaptors::Mouse_Single_Sample->insert(\%mouse_single_sample_info);
		print STDERR "\t\t\tCreate single_sample $sample_name_local with id $mouse_single_sample_db_id\n";   
			     	
	} else {
    	#Get the existing sample
    	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name,source_group_id=>$sourcegroup_db_id);
    	print "Sample $sample_name source_groupdb $sourcegroup_db_id\n";
    	$sample_db_id = $sample_obj->id;
    }	
	    
	    
  
    
    my @bams = ();			
    my $lanes_outdir = $outdir.'/'.$sample_name_local.'_lanes';
    my $runs_outdir = $outdir.'/'.$sample_name_local.'_runs';
    
    my $encoding;
	#Create the lanes for each sample
    for ( my $lane_count = 0 ; $lane_count < $lane_total ; $lane_count++ ) {
    	my $lane_number = $lane_count;
    	$lane_number++;
    	my $lane_name = $sample_name_local.'_l'.$lane_number;
        
        #Get the read_file info
        my $read1_compressed = 0;
        my $read2_compressed = 0;
        my $read1_suffix;
        my $read2_suffix;
		my $read1_name = $lane_name.'_r1';
		my $read2_name = $lane_name.'_r2';
		
		
		my $read1_actual = $lane_data->{$lane_count}{read1};
		my $read2_actual = $lane_data->{$lane_count}{read2};
		
		my $read1_full_file = $read_directory . '/'. $lane_data->{$lane_count}{read1};
		my $read2_full_file = $read_directory . '/'. $lane_data->{$lane_count}{read2};
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
								read_directory=>$read_directory,
								read_file_name=>$read1_name,
								lane_id=>$lane_db_id
								);
								
			my %read2_info = (
								file_name => $lane_data->{$lane_count}{read2},
								is_compressed=>$read2_compressed,
								compression_suffix=>$read2_suffix,
								read_file_number=>2,
								read_directory=>$read_directory,
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
			($encoding,$lane_bam) = $pipeline_sample->align_lane_qsub(-sample_name=>$sample_name_local,-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id,-commands=>\%commands);
		} else {
			($encoding,$lane_bam) = $pipeline_sample->align_lane_qsub(-sample_name=>$sample_name_local,-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id);
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
    
    my $steps_xml_out = $steps_xml_dir.$sample_name_local.'.xml';
    #Assumes same encoding for a sample but likely ok
	$pipeline_sample->create_sample_xml(-sample_name=>$sample_name_local,-steps_xml_template=>$steps_xml,-xml_out=>$steps_xml_out,-encoding=>$encoding,-bams=>\@bams);

    #Now create the run qsubs; only a single block for each run here; not worth splitting 
	$pipeline_sample->create_run_qsubs(-sample_name=>$sample_name_local,-qsub_dir=>$steps_xml_dir,-steps_xml=>$steps_xml_out,-start_step=>$first_pipe_step,-end_step=>$last_pipe_step,-pipe_block=>1);	  	    
    	
}












