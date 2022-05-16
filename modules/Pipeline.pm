package modules::Pipeline;

#modules contains generic functions for running the pipeline; no objects required

use strict;
use modules::QualityEncoding;
use modules::Exception;
use Data::Dumper;
use modules::Cluster;
use modules::ConfigXML;
use modules::Adaptors::Syscall;
use modules::Adaptors::Source;
use modules::Adaptors::Pipeline_Step;
use modules::Adaptors::Sample;
use modules::Adaptors::Human_Cancer_Sample;
use modules::Adaptors::Mouse_Cancer_Sample;
use modules::Adaptors::Human_Related_Sample;
use modules::Adaptors::Project;
use modules::Adaptors::Source_Group;


sub new  
{
	my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

    return $self;
}

#Get source id from sample_name or source_group_name
sub get_source_name
{
	my ($sample_name) = @_;
	if ($sample_name =~ /^(.+)_sg[0-9]_\w+$/) {
		return $1;
	} elsif ($sample_name =~ /^(.+)_sg[0-9]$/) {
		return $1;
	} else {
		modules::Exception->throw("ERROR: Can't get source name from $sample_name");
	}
}

#Get the sample_type from the run_id or the sample_name
sub get_sample_type
{	
    my %args = @_;
	
	my $sample_obj;
	
	if (exists $args{-run_id}) {
		my ($run_obj) = modules::Adaptors::Run->search(id=>$args{-run_id});
	
		if (!defined $run_obj) {
			modules::Exception->throw("ERROR: Can't retrieve run db obj for run id $args{-run_id}");
		}
		
		($sample_obj) = modules::Adaptors::Sample->search(id=>$run_obj->sample_id);
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-run_id}");
		}
	
		
	} elsif (exists $args{-sample_name}) {
				
		($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$args{-sample_name});
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-sample_name}");
		}
		
	}
	return $sample_obj->sample_type;
}


#Get the sequence_type from the run_id or the sample_name
sub get_sequence_type
{	
    my %args = @_;
	
	my $sample_obj;
	
	if (exists $args{-run_id}) {
		my ($run_obj) = modules::Adaptors::Run->search(id=>$args{-run_id});
	
		if (!defined $run_obj) {
			modules::Exception->throw("ERROR: Can't retrieve run db obj for run id $args{-run_id}");
		}
		
		($sample_obj) = modules::Adaptors::Sample->search(id=>$run_obj->sample_id);
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-run_id}");
		}
	
		
	} elsif (exists $args{-sample_name}) {
				
		($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$args{-sample_name});
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-sample_name}");
		}
		
	}
	return $sample_obj->sequence_type;
}

sub get_sample_name
{
	my %args = @_;
	
	my $sample_obj;
	
	if (exists $args{-run_id}) {
		my ($run_obj) = modules::Adaptors::Run->search(id=>$args{-run_id});
	
		if (!defined $run_obj) {
			modules::Exception->throw("ERROR: Can't retrieve run db obj for run id $args{-run_id}");
		}
		
		($sample_obj) = modules::Adaptors::Sample->search(id=>$run_obj->sample_id);
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-run_id}");
		}
	
		
	} else {
		modules::Exception->throw("ERROR: require argument -run_id");
	}
	return $sample_obj->sample_name;
}


#Get source type from run, sample_name, or source_group_name
sub get_source_type
{
    my %args = @_;
	
	my $source_name;
	
	if (exists $args{-run_id}) {
		my ($run_obj) = modules::Adaptors::Run->search(id=>$args{-run_id});
	
		if (!defined $run_obj) {
			modules::Exception->throw("ERROR: Can't retrieve run db obj for run id $args{-run_id}");
		}
		
		my ($sample_obj) = modules::Adaptors::Sample->search(id=>$run_obj->sample_id);
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-run_id}");
		}
	
		($source_name) = $sample_obj->sample_name =~ /^(.+)_sg[0-9]_\w+$/;
		
	} elsif (exists $args{-sample_name}) {
		($source_name) = $args{-sample_name} =~ /^(.+)_sg[0-9]_\w+$/;
	} elsif (exists $args{-source_group_name}) {
		($source_name) = $args{-source_group_name} =~ /^(.+)_sg[0-9]$/;
	} elsif (exists $args{-source_name}) {
		$source_name = $args{-source_name};
	} else {
		modules::Exception->throw("ERROR: Problem with arguments: require -run_id, -source_name, -sample_name, or -source_group_name");
	}
		
	my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source_name);
	
	if (!defined $source_obj) {
		modules::Exception->throw("ERROR: Can't get db source obj for source name $source_name");
	}
	return $source_obj->source_type;
	
}

#Get source id from sample_name or source_group_name
sub get_source_group_name
{
	my ($sample_name) = @_;
	if ($sample_name =~ /^(.+_sg[0-9])_\w+$/) {
		return $1;
	} elsif ($sample_name =~ /^(.+_sg[0-9])$/) {
		return $1;
	} else {
		modules::Exception->throw("ERROR: Can't get source name from $sample_name");
	}	
}

#Get the latest stepname for a run
sub get_latest_step {
	my ($run_id) = @_;
	my $last_step_run = 'no_steps_recorded';
	my @pipeline_steps_obj = modules::Adaptors::Pipeline_Step->search_all();
								
	for my $pipeline_step_obj (@pipeline_steps_obj) {
		my ($pipeline_step_run_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$pipeline_step_obj->id);
		if (defined $pipeline_step_run_obj) {
			$last_step_run = $pipeline_step_obj->name;
		}
	}
	return $last_step_run;
}

#Get the next step to run
sub get_next_step {
	my ($run_id,$xml,$pipe_config,$source_type) = @_;
	my $next_step_run;
	my %steps_runs = ();
	my @pipeline_steps_obj = modules::Adaptors::Pipeline_Step->search_all();
								
	for my $pipeline_step_obj (@pipeline_steps_obj) {
		my ($pipeline_step_run_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$pipeline_step_obj->id);
		if (defined $pipeline_step_run_obj) {
			$steps_runs{$pipeline_step_obj->name}++;
		}
	}
	
	my $sample_xml = modules::ConfigXML->new($xml);
	my $steps = $sample_xml->read('steps_order', 'step');
	my @steps;
	#Get the steps to run	
	if (ref($steps) eq 'ARRAY'){ # Cope with steps being a single step or an array of steps
	    @steps = @$steps;
	} else {
	    @steps = ($steps);
	}
	
	
	#Find the first occurence that hasn't been run from the xml
	my $last_step = $steps[0];
	for my $step (@steps) {
		my $by_chr  = 0;
		if ($sample_xml->read('steps', 'step', $last_step, 'by_chr')) {
			$by_chr = 1;
		}
		#If it's by_chr need to return the by_chr step; otherwise return the next step
		if (!exists $steps_runs{$step}) {
			if ($by_chr) {
				my @chr = split(" ",$pipe_config->read($source_type,'annotation_version','chr'));
				my ($current_pipeline_step_obj) = modules::Adaptors::Pipeline_Step->search(name=>$last_step);
				my @pipeline_step_run_obj = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$current_pipeline_step_obj->id);
				
				if(@pipeline_step_run_obj != @chr){
					return $last_step;
				}else{
					return $step;
				}
			} else {
				return $step;
			}		
		}
		$last_step = $step;
	}
}

#Get the commands from a qsub file
sub get_commands {
	my ($qsub_file) = @_;
	
	if ( !-e $qsub_file ) {
		modules::Exception->throw("File $qsub_file doesn't exist");	
	}
	
	my @commands;
	
	open(FILE,"$qsub_file") || modules::Exception->throw("Can't open file $qsub_file\n");
	my $command_flag = 0;
	while (<FILE>) {
		chomp;
		push @commands, $_ if $command_flag;
		if (/source/) {
			$command_flag = 1;
		}
	}
	return \@commands;
}

#Move the run files to subdirectories from workdir
sub move_run_files {
	my ($workdir) = @_;
   	my $overlap_files = $workdir .'/*sg*';
   	my $vep_files = $workdir .'/vep*';
   	my $overlapdir = $workdir.'/overlap';
   	my $overlap_flag = `ls $overlap_files 2>/dev/null`;
   	my $vep_flag = `ls $vep_files 2>/dev/null`;
   	system("mv $overlap_files $overlapdir") if $overlap_flag;
   	system("mv $vep_files $overlapdir") if $vep_flag;
}


#Move the run files back to workdir from subdirectories
sub move_run_files_back {
	my ($workdir) = @_;
   	my $overlap_files = $workdir .'/overlap/*sg*'; #Remainder goes to overlap
   	my $vep_files = $workdir .'/overlap/vep*';
   	my $overlap_flag = `ls $overlap_files 2>/dev/null`;
   	my $vep_flag = `ls $vep_files 2>/dev/null`;
   	system("mv $overlap_files $workdir") if $overlap_flag;
   	system("mv $vep_files $workdir") if $vep_flag;
}

#Get lane info from a directory containing the full set of fastq input data or aligned BAM files
sub get_lane_info {
	my ($read_directory,
	    $read1_regex) = @_;
	
	my %lane_info;

	# Try to spoof the read2 regex
	my $read2_regex;
	if ($read1_regex =~ /\w/) {
	    #Replace the last occurrence of '1' with 2
	    ($read2_regex = $read1_regex) =~  s/(.*)1/${1}2/;
	} else {
	    $read1_regex = 'R1';
	    $read2_regex = 'R2';
	}
	
	if ($read1_regex eq $read2_regex) {
		modules::Exception->throw("ERROR: Read regexes match $read1_regex $read2_regex");
	}
	
	# Read input directory, check that read files are paired nicely and count the number of distinct pair.  Each pair is a lane.
	opendir(DIR,"$read_directory") || modules::Exception->throw("ERROR: Can't open readdir $read_directory $!\n");
	my @files = readdir DIR;
	my @reads1_all = sort grep {/$read1_regex/} @files;
	my @reads2_all = sort grep {/$read2_regex/} @files;
	closedir(DIR);
	
	if (!@reads1_all) {
		modules::Exception->throw("ERROR: No reads match regex $read1_regex");
	} elsif (!@reads2_all) {
		modules::Exception->throw("ERROR: No reads match regex $read2_regex");
	} elsif (@reads1_all != @reads2_all) {
		print "ERROR with read numbers\n";
		print Dumper \@reads1_all;
		print Dumper \@reads2_all;
		modules::Exception->throw("ERROR: Read numbers don't match");
	}
	
	for ( my $count = 0 ; $count < @reads1_all; $count++ ) {
	    $lane_info{$count}{read1} = $reads1_all[$count];
	    $lane_info{$count}{read2} = $reads2_all[$count];
	}
	return (\%lane_info);
}

#Insert the step_command and bin_vers entry if required
sub insert_step_command {
	my ($stepname,$command,$runid,$chr) = @_;
	
	
	
	#Now create the pipeline_step_run entry
	my ($stepobj) = modules::Adaptors::Pipeline_Step->search('name'=>$stepname);
	if (! defined $stepobj) {
		modules::Exception->throw("Step $stepname doesn't exist in the database");
	}
	my $stepid = $stepobj->id;
	#Insert the step_command entry
    my %step_command = (
    					'pipeline_step_id' => $stepid,
    					'run_id' => $runid
    					);
    my $pipeline_step_run_id = modules::Adaptors::Pipeline_Step_Run->insert(\%step_command);		
    
    #First create the syscall entry
	my %syscall_entry = (
						'command' => $command,
						'level' => 'sample',
						'pipeline_steps_run_id'=>$pipeline_step_run_id,
						'analysis_step_name'=>$stepname
						);
	
	#If it's a chromosome specific command					
	if (defined $chr) {
		$syscall_entry{'chr'} = $chr;
	}	
					
	my $syscall_id = modules::Adaptors::Syscall->insert(\%syscall_entry);	
	
    return $pipeline_step_run_id;
}

#Gets the conf dir for writing out xmls
sub get_run_dir {
	my %args = @_;
	
	my @required_args = (
			             '-sample_name',
			             '-run_id'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $sample_name = $args{-sample_name};
    my $run_id = $args{-run_id};
    my $cluster_config = modules::Pipeline::get_cluster_conf();
	my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);
	my $source_name = modules::Pipeline::get_source_name($sample_name);
	my $rundir_base = $cluster_config->read($source_type,'base_directories','base_run_directory');
	my $final_rundir = $rundir_base . '/' . $source_name . '/' .$sample_name .'_runs/'. $sample_name .'_'.$run_id;
	if ( !-d $final_rundir ) {
		modules::Exception->throw("Directory $final_rundir doesn't exist");	
	}
	return $final_rundir;
}

sub get_project_name {
	my %args = @_;
	
	
	my $source_obj;
	
	if (exists $args{-run_id}) {
		my ($run_obj) = modules::Adaptors::Run->search(id=>$args{-run_id});
	
		if (!defined $run_obj) {
			modules::Exception->throw("ERROR: Can't retrieve run db obj for run id $args{-run_id}");
		}
		
		my ($sample_obj) = modules::Adaptors::Sample->search(id=>$run_obj->sample_id);
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-run_id}");
		}
	
		my ($source_group_obj) = modules::Adaptors::Source_Group->search(id=>$sample_obj->source_group_id);
		
		if (!defined $source_group_obj) {
			modules::Exception->throw("ERROR: Can't retrieve source_group db obj for run id $args{-run_id}");
		}
		
		($source_obj)  = modules::Adaptors::Source->search(id=>$source_group_obj->source_id);
		
	} elsif (exists $args{-sample_name}) {
		my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$args{-sample_name});
		
		if (!defined $sample_obj) {
			modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-sample_name}");
		}
		
		my ($source_group_obj) = modules::Adaptors::Source_Group->search(id=>$sample_obj->source_group_id);
		
		if (!defined $source_group_obj) {
			modules::Exception->throw("ERROR: Can't retrieve source_group db obj for run id $args{-sample_name}");
		}
		
		($source_obj)  = modules::Adaptors::Source->search(id=>$source_group_obj->source_id);
		
	} elsif (exists $args{-source_group_name}) {
		my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$args{-source_group_name});
		
		if (!defined $source_group_obj) {
			modules::Exception->throw("ERROR: Can't retrieve source_group db obj for run id $args{-sample_name}");
		}
		
		($source_obj)  = modules::Adaptors::Source->search(id=>$source_group_obj->source_id);
		
		
	} elsif (exists $args{-source_name}) {
		($source_obj)  = modules::Adaptors::Source->search(external_source_name=>$args{-source_name});
		
	}
	
	if (!defined $source_obj) {
		modules::Exception->throw("ERROR: Can't retrieve source db obj");
	}
	my ($project_obj) = modules::Adaptors::Project->search(id=>$source_obj->project_id);
	
	return $project_obj->project_name;
	
	
}

sub get_cluster_conf {
   	my $self = shift;
    my %args = @_;	
	my $svndir;
	if (!exists $ENV{'SVNDIR'}) {
		modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
	} else {
		$svndir = $ENV{'SVNDIR'};
	}

	#Get the xml files and create the pipeline object
        my $cluster_xml;

        if(defined $args{-cluster_xml}){ # this allows usage of custom cluster_xml file
                $cluster_xml = $args{-cluster_xml};
        } else {
                $cluster_xml = "$svndir/conf/cluster.xml";
        }
 
	if ( !-e $cluster_xml ) {
		modules::Exception->throw("File $cluster_xml doesn't exist");	
	}
	
	my $cluster_config = modules::ConfigXML->new($cluster_xml);
	
	return $cluster_config;
}

sub get_pipe_conf {
	my $svndir;
	if (!exists $ENV{'SVNDIR'}) {
		modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
	} else {
		$svndir = $ENV{'SVNDIR'};
	}
	
	#Get the xml files and create the pipeline object
	my $pipe_xml = "$svndir/conf/pipe.xml";

	if ( !-e $pipe_xml ) {
		modules::Exception->throw("File $pipe_xml doesn't exist");	
	}
	
	my $pipe_config = modules::ConfigXML->new($pipe_xml);
	
	return $pipe_config;
}

sub get_sample_conf {
	my $svndir;
	if (!exists $ENV{'SVNDIR'}) {
		modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
	} else {
		$svndir = $ENV{'SVNDIR'};
	}
	
	#Get the xml files and create the pipeline object
	my $sample_file = "$svndir/conf/sample_info.csv";

	if ( !-e $sample_file ) {
		modules::Exception->throw("File $sample_file doesn't exist");	
	}	
	
	return $sample_file;
}

sub get_report_conf {
   	my $self = shift;
    my %args = @_;
	my $svndir;
	if (!exists $ENV{'SVNDIR'}) {
		modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
	} else {
		$svndir = $ENV{'SVNDIR'};
	}
	
	#Get the xml files and create the pipeline object
        my $report_xml;

        if(defined $args{-report_xml}){ # this allows usage of custom report_xml file
                $report_xml = $args{-report_xml};
        } else {
                $report_xml = "$svndir/conf/report.xml";
        }

        if ( !-e $report_xml ) {
                modules::Exception->throw("File $report_xml doesn't exist");
        }

        my $report_config = modules::ConfigXML->new($report_xml);

        return $report_config;	
}

#Check whether a sample is whole genome and not cancer as these regular whole genome cases are handled differently by the database
sub get_whole_genome {
	#modules::Pipeline::get_whole_genome(-sample_name=>$sample_name)
    my %args = @_;
	
    my @required_args = (
			             '-sample_name'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $sample_name = $args{-sample_name};
	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);
	if (!$sample_obj) {
		modules::Exception->throw("ERROR: Cannot retrieve sample object for $sample_name");
	}
	
	my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_name);
	my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);
	
	if ($sequence_type eq 'genome') {
		return 1;	
	} else {
		return 0;
	}
	
	
}

#Check whether a sample is a tumour from sample_name
sub get_tumour_flag {
    my %args = @_;
	
    my @required_args = (
			             '-sample_name'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $sample_name = $args{-sample_name};
	my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);
	if (!$sample_obj) {
		modules::Exception->throw("ERROR: Cannot retrieve sample object for $sample_name");
	}
	
	my ($human_cancer_obj) = modules::Adaptors::Human_Cancer_Sample->search(sample_id=>$sample_obj->id);
	my ($mouse_cancer_obj) = modules::Adaptors::Mouse_Cancer_Sample->search(sample_id=>$sample_obj->id);
	
	if ($human_cancer_obj && $human_cancer_obj->tumour) {
		return 1;	
	} elsif ($mouse_cancer_obj && $mouse_cancer_obj->tumour) {
		return 1;	
	}else {
		return 0;
	}
	
	
}

#Get cancer type from run_id or sample_name
sub get_cancer_type {
	my %args = @_;
	
	
	my $sample_obj;
	
	if (exists $args{-run_id}) {
		my ($run_obj) = modules::Adaptors::Run->search(id=>$args{-run_id});
	
		if (!defined $run_obj) {
			modules::Exception->throw("ERROR: Can't retrieve run db obj for run id $args{-run_id}");
		}
		($sample_obj) = modules::Adaptors::Sample->search(id=>$run_obj->sample_id);
	} elsif (exists $args{-sample_name}) {
		($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$args{-sample_name});
	}
	
	if (!defined $sample_obj) {
		modules::Exception->throw("ERROR: Can't retrieve sample db obj for run id $args{-sample_name}");
	}

	my ($human_cancer_obj) = modules::Adaptors::Human_Cancer_Sample->search(sample_id=>$sample_obj->id);
	my ($mouse_cancer_obj) = modules::Adaptors::Mouse_Cancer_Sample->search(sample_id=>$sample_obj->id);
	
	if ($human_cancer_obj) {
		return $human_cancer_obj->cancer_type;	
	} elsif ($mouse_cancer_obj) {
		return $mouse_cancer_obj->cancer_type;	
	}else {
		my $sample_name = $sample_obj->sample_name;
		modules::Exception->throw("ERROR: Can't get cancer type info for sample $sample_name");
	}
	
}


#Get the special pipeline step from the xml label
sub get_xml_step {
    my %args = @_;
	
	my $pipe_conf = get_pipe_conf(); 
	
    my @required_args = (
			             '-source_type',
			             '-xml_name'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $source_type = $args{-source_type};
	my $xml_name = $args{-xml_name};

	if ($pipe_conf->exists($source_type,'pipe_steps',$xml_name)) {
		#If the source type has a non default value
		return $pipe_conf->read($source_type,'pipe_steps',$xml_name);
	} elsif ($pipe_conf->exists('common','pipe_steps',$xml_name)) {
		return $pipe_conf->read('common','pipe_steps',$xml_name);
	} else {
		modules::Exception->throw("ERROR: No xml label match source_type $source_type");
	}
}


return 1;
