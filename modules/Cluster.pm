package modules::Cluster;

use strict;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;
use modules::Utils;

sub new {
	my ($class, @args) = @_;
    my %args = @args;
    my $self = bless {}, $class;
	my @required_args = (
			             -scheduler,
			             -svndir
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
    }

	my $scheduler = $args{-scheduler};
	my $svndir = $self->svndir($args{-svndir});

    if ($scheduler eq 'PBS') {
		$self->{pbs} = 1;	
    } elsif ($scheduler eq 'SGE') {
    	$self->{sge} = 1;
    } else {
		modules::Exception->throw("Only PBS and SGE support at this time");
    }

	return $self;
}

#Subroutine to create the qsub file
sub create_single_qsub {
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -qsub_dir,
			             -qsub_file,
			             -commands
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $qsub_dir = $args{-qsub_dir};
	if ( !-d $qsub_dir ) {
		system("umask 0007; mkdir -p $qsub_dir");
	}
	my $qsub_file = $args{-qsub_file};
	my $qsub_out = $qsub_dir . '/' . $qsub_file;
	my ($date_stamp) = split('_',modules::Utils->GetTime());
	my $qsub_error = $qsub_out . '.err.'.$date_stamp;
	my $qsub_output = $qsub_out . '.out.'.$date_stamp;
	my $qsub_next = defined $args{-qsub_next_file}?$args{-qsub_next_file}:0;
	
	my $bash = defined $args{-bash}?$args{-bash}:"#!/bin/bash";
	#Add shebang if not included
	if ($bash !~ /^#/) {
		$bash = "#!".$bash;
	}
	my $queue = defined $args{-queue}?$args{-queue}:"normal";
	my $project = defined $args{-project}?$args{-project}:"u86";
	
	my $walltime = defined $args{-walltime}?$args{-walltime}:"24:00:00";
	if ($walltime !~ /\d+:\d\d:\d\d/) {
		modules::Exception->throw("ERROR: Walltime arguments must be in the format [H]H:MM:SS\n");
	}
	my $jobfs = defined $args{-jobfs}?$args{-jobfs}:"10GB";
	my $mem = defined $args{-mem}?$args{-mem}:"4GB";
	$mem = uc($mem);
	
	if ($mem !~ /GB$/ && $mem !~ /MB$/) {
		modules::Exception->throw("ERROR: Memory must end in GB or MB\n");
	}
		
	my $cpus = defined $args{-cpus}?$args{-cpus}:1;
	if ($cpus !~ /^\d+$/) {
		modules::Exception->throw("ERROR: cpus must be an integer\n");
	}
	
	my $modules = defined $args{-modules}?$args{-modules}:"";
	my $commands = $args{-commands};
	if (ref($commands) ne 'ARRAY') {
		modules::Exception->throw("ERROR: commands args must be array ref");
	} 
	my $command_str = join("\n",@{$commands});

	my $source_line = 'source '. $self->svndir(). '/conf/export_env.txt '.$self->svndir();

	my $resource_line;

	open(QSUB,">$qsub_out") || modules::Exception->throw("ERROR: Cannot open qsub directory $!");
	print QSUB "$bash\n";
	my $qsub_out_single = $qsub_out . '.single';
	open(SINGLE,">$qsub_out_single") || modules::Exception->throw("ERROR: Cannot open qsub directory $!");
	print SINGLE "$bash\n";


	#Now create the qsub file
	#Sean Li added -W umask=0007 and -W group_list=u86 for err and out files from PBS to be read- and writeable by group members
	if ($self->{pbs}) {
		$resource_line = "walltime=$walltime,mem=$mem,jobfs=$jobfs,ncpus=$cpus";
		print QSUB "#PBS -P $project\n";
		print QSUB "#PBS -q $queue\n";
		print QSUB "#PBS -l $resource_line\n";
		print QSUB "#PBS -l other=gdata2\n"; 
		print QSUB "#PBS -W umask=0007\n"; 
		print QSUB "#PBS -W group_list=u86\n"; 
		print QSUB "#PBS -o $qsub_output\n";
		print QSUB "#PBS -e $qsub_error\n";
		print SINGLE "#PBS -P $project\n";
		print SINGLE "#PBS -q $queue\n";
		print SINGLE "#PBS -l $resource_line\n";
		print SINGLE "#PBS -l other=gdata2\n"; 
		print SINGLE "#PBS -W umask=0007\n"; 
		print SINGLE "#PBS -W group_list=u86\n"; 
		print SINGLE "#PBS -o $qsub_output\n";
		print SINGLE "#PBS -e $qsub_error\n";
		if ($modules =~ /module/) {
			if ($modules =~ /java/) {
				print QSUB "#PBS -l other=hyperthread\n";
				print SINGLE "#PBS -l other=hyperthread\n";
			}
			print QSUB "$modules\n";
			print SINGLE "$modules\n";
		}
		print QSUB "$source_line\n";
		print SINGLE "$source_line\n";
		print QSUB "$command_str\n";
		print SINGLE "$command_str\n";
		if ($qsub_next) {
			my $qsub_next_out = $qsub_dir . '/' . $qsub_next;
			print QSUB "qsub $qsub_next_out\n"; #only do this for pipeline code
		}			
	} elsif ($self->{sge}) {
		$resource_line = "h_rt=$walltime,virtual_free=$mem\n#\$ -pe threads $cpus";
		print QSUB "#\$ -cwd\n";
		print QSUB "#\$ -S $bash\n";
		print QSUB "#\$ -q $queue\n";
		print QSUB "#\$ -l $resource_line\n";
		print QSUB "$source_line\n";
		print QSUB "$command_str\n";			
	}
	close QSUB;
	close SINGLE;
}

#Subroutine to create the bash wrapper file with dependent job submission
sub create_qsub_wrapper {
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -qsub_dir,
			             -qsub_file,
			             -jobs
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $qsub_dir = $args{-qsub_dir};
	if ( !-d $qsub_dir ) {
		system("mkdir -p $qsub_dir");
	}
	my $qsub_file = $args{-qsub_file};
	my $qsub_out = $qsub_dir . '/' . $qsub_file;
	my $jobs = $args{-jobs};
	if (ref($jobs) ne 'ARRAY') {
		modules::Exception->throw("ERROR: jobs args must be array ref");
	} 
	my @full_jobs = map {$_ = $qsub_dir .'/'.$_} @{$jobs};
	my $first_job = shift @full_jobs;
	my $dependent_jobs = join(" ",@full_jobs);
		
	
	my $bash = defined $args{-bash}?$args{-bash}:"#!/bin/bash";
	#Add shebang if not included
	if ($bash !~ /^#/) {
		$bash = "#!".$bash;
	}
	
	#Now create the qsub wrapper file
	open(WRAPPER,">$qsub_out") || modules::Exception->throw("ERROR: Cannot open qsub directory $!");
	print WRAPPER "$bash\n";
	print WRAPPER "jobid=`qsub $first_job`\n";
	print WRAPPER "regex=\"([0-9]+)\"\n";
	print WRAPPER "jobstr=''\n";
	print WRAPPER "if [[ \$jobid =~ \$regex ]]\nthen\n";
	print WRAPPER "\tfor job in $dependent_jobs\n";	
	if ($self->{pbs}) {
		print WRAPPER "\t\n\tdo\n\t\tif [[ \$jobid =~ \$regex ]]\n\t\tthen\n\t\t\tif [ -z \$jobstr ]\n\t\t\tthen\n\t\t\t\tjobstr=\${BASH_REMATCH[1]}\n\t\t\telse\n\t\t\t\tjobstr+=\":\"\${BASH_REMATCH[1]}\n\t\t\tfi\n\t\t\tjobid=`qsub -W depend=afterok:\$jobstr \$job`\n\t\telse\n\t\t\techo 'ERROR: Later submission did not return jobid' \n\t\tfi\n\tdone\nelse\n\techo 'ERROR: First submission did not return jobid'\nfi\n";
	} elsif ($self->{sge}) {
		print WRAPPER "\t\n\tdo\n\t\tif [[ \$jobid =~ \$regex ]]\n\t\tthen\n\t\t\tjobid=`qsub -hold_jid \${BASH_REMATCH[1]} \$job`\n\t\telse\n\t\t\techo 'ERROR: Later submission did not return jobid' \n\t\tfi\n\tdone\nelse\n\techo 'ERROR: First submission did not return jobid'\nfi\n";
	}
	close WRAPPER;
	
}

#Set svndir
sub svndir
{
    my ($self, $svndir) = @_;

    if (defined $svndir) {
		$self->{'svndir'} = $svndir;
    } elsif (! defined $self->{'svndir'}) {
		modules::Exception->throw("svndir not set");
    }

    return $self->{'svndir'};
}

sub get_job_list {
	#Get a list of running jobs
	my ($self) = @_;
	
	my %jobs = ();
	
	if ($self->{pbs}) {
		my $qstat_out = `qstat -f 2>&1`;
		return 0 unless $qstat_out;
		my @qstat_lines = split("\n",$qstat_out);
		my $job_id;
		my $job_state;
		my $job_name;
		my $sample_found = 0;
		my $sample_name = '';
		my $new_job = 0;

		#marcin modified so it records cohort summaries; it's really crappy code, should be re-written
    for my $qstat_line(@qstat_lines) 
    {
	    if($qstat_line =~ /Job Id: (\S+)/) 
	    {
	    	$job_id = $1;
	    	$new_job = 1;
	    	$sample_found = 0;
	    	$job_state = '';
			}
			elsif($qstat_line =~ /job_state = (\S+)/) 
			{
	    	$job_state = $1;
	    }
	    elsif($qstat_line =~ /([^\/\.\s]+_sg\d+_[^\/\.\s]+)/)
	    {
	    	#Matches vc_cohort2_sg1_affected1, vc_cohort2_sg1_unaffected2,vc60310_sg1_humansingle1, SCC10-0092_32147_sg1_humansingle1, etc
	     	$sample_found = 1;
				$sample_name = $1;
	   	}
	    elsif($qstat_line =~ /([^\/\.\s]+_cohort).qsub/) #marcin: that's the cohort summary, sample name will be <sample_ID>_cohort, eg: APOSLE_cohort93_cohort
	    {
	    	#Matches vc_cohort2_cohort.qsub, SCC10-0092_32147_cohort.qsub, ...
	     	$sample_found = 1;
				$sample_name = $1;
	   	}
	   	elsif($qstat_line =~ /(MOUSE_cancer\d+_sg\d+_[a-z]+\d+)/) 
	   	{
	    	$sample_found = 1;
				$sample_name = $1;
	    }
	    
	    if ($sample_found && $new_job) 
	    {
	    	$jobs{$sample_name} = $job_id;
	    }
		}
	}
	elsif($self->{sge})
	{
		#TODO: Write for sge
	}
	$self->{jobs} = \%jobs;
}

#Checks whether merge_vcf job is running already or not; required for resume_run to avoid submitting same job twice
sub check_job_running
{
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -sample_name
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $sample_name = $args{-sample_name};
	if (exists $self->{jobs}->{$sample_name}) {
		return $self->{jobs}->{$sample_name};
	} else {
		return 0;	
	}
	
	
	
}

#Checks whether merge_vcf job is running already or not; required for resume_run to avoid submitting same job twice
sub check_lanejob_running
{
	my ($self, @args) = @_;
	
	my %args = @args;

	my @required_args = (
			             -lane_name
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}	
	}
	
	my $lane_name = $args{-lane_name};
	if (exists $self->{jobs}->{$lane_name}) {
		return $self->{jobs}->{$lane_name};
	} else {
		return 0;	
	}
	
	
	
}


return 1;
