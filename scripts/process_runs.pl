#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::Adaptors::Source;
use modules::Adaptors::Align_Lane;
use modules::Adaptors::Lane;
use modules::Adaptors::Syscall;
use modules::Adaptors::Lane_Run;
use modules::Adaptors::Sample;
use modules::Adaptors::Release_File;
use modules::ConfigXML;
use modules::Pipeline;
use File::Basename;
use Pod::Usage;
use Cwd;
use Cwd 'abs_path';
use vars qw(%OPT);

GetOptions(\%OPT, 
		   	"help|h",
		   	"man|m",
		   	"submit",
		   	"restart"
		   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

process_run.pl -submit submit_jobs -restart restart_incomplete_jobs_without_existing_jobs_running_or_queued

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

process_run.pl -> Start/resume pipeline runs when certain conditions are met

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

commence_run.pl 

=cut  

#Keep the svndir 
my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $restart = defined $OPT{restart}?1:0;



#Get the xml files and create the pipeline object
my $cluster_xml = "$svndir/conf/cluster.xml";

if ( !-e $cluster_xml ) {
	modules::Exception->throw("File $cluster_xml doesn't exist");	
}

my $cluster_config = modules::ConfigXML->new($cluster_xml);

my @commands = ();

#First retrieve all the samples
my @samples_obj = modules::Adaptors::Sample->search_all();
my $unaligned_lanes = 0;
my $first_pipe_step = modules::Pipeline::get_xml_step(-source_type=>'common',-xml_name=>'first'); 

my $scheduler = $cluster_config->read('common','scheduler');
my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
$cluster_obj->get_job_list();
my %sample_missing_lanes = ();


print "Runs ready to start....\n";

for my $sample_obj ( @samples_obj ) {
    my $sample_name = $sample_obj->sample_name;
    
    
    my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);
	
    #Get the total lanes variable
    my $total_lanes = $sample_obj->total_lanes;

	#Get the number of lanes db entries and ensure they match    
    my @lanes_objs = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);
    
    if (@lanes_objs != $total_lanes) {
    	my $lane_obj_count = @lanes_objs;
    	modules::Exception->throw("ERROR: Sample $sample_name has 'total_lanes' of $total_lanes and $lane_obj_count db lanes entries; These should be the same");
    }
    
    #Check if each lane has a lane_align entry
    my $align_lane_total = 0;
    my @unaligned_lanes = ();
    for my $lane_obj ( @lanes_objs ) {
        my $lane_name = $lane_obj->lane_name;
        my @align_lane = modules::Adaptors::Align_Lane->search(lane_id=>$lane_obj->id);
        if (@align_lane > 1) {
        	modules::Exception->throw("ERROR: Lane number $lane_name has more than one align_lane entry");
        } elsif (@align_lane == 1) {
        	$align_lane_total++;
        } else {
        	push @unaligned_lanes, $lane_name;
        	my $job_id = $cluster_obj->check_lanejob_running(-lane_name=>$lane_name);
						
			if (!$job_id) {
				print STDERR "ERROR: Lane $lane_name isn't aligned and no jobs are queued\n";
				next;
			} 
        }
    }
    
    if ($align_lane_total == $total_lanes) {
	    #Check if all these lanes are part of existing run
	    #We need to trigger new run if there exists at least one lane not used in previous run (either completely new or new lane added)  
	    my $run_lane_total = 0;
	    for my $lane_obj ( @lanes_objs ) {
			my @lane_run_objs = modules::Adaptors::Lane_Run->search(lane_id=>$lane_obj->id);
			if (@lane_run_objs) {
				$run_lane_total++;
			}
		}
		
		#If we need to create a new run
		if ($run_lane_total != $total_lanes) {
			my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
						
			if ($job_id) {
				#print STDERR "Sample $sample_name is ready to create run but $job_id is queued\n";
				next;
			} else {
				print STDERR "Begin new run for sample $sample_name\n";
			}
			
			my $qsub_base = $cluster_config->read($source_type,'base_directories','base_qsub_directory');
			my ($patient_name) = modules::Pipeline::get_source_name($sample_name);
			my $qsub_dir = $qsub_base . '/' . $patient_name . '/runs/';
			my $qsub_first = $qsub_dir. $sample_name . '.pipe1.'.$first_pipe_step.'.qsub';
			if (!-e $qsub_first) {
				modules::Exception->throw("ERROR: No qsub file for $sample_name (file = $qsub_first)");
			}
			push @commands, "qsub $qsub_first";
			
		} else {
			#print STDERR "Skip sample $sample_name: Run with these lanes already exists\n";
		}
	
	
    	
    } else {
    	my $lane_count = @unaligned_lanes;
    	my $lane_str = join(", ",@unaligned_lanes);
    	print STDERR "Skip sample $sample_name: There are $lane_count lanes unaligned: ($lane_str)\n";
    	$unaligned_lanes += $lane_count;
    	$sample_missing_lanes{$sample_name}++;
    }
    
	
	    
    
}

print "\nTotal unaligned lanes: $unaligned_lanes\n" unless $unaligned_lanes == 0;


print "\nRuns to progress...\n";

my %step_ids = ();

my @step_objs = modules::Adaptors::Pipeline_Step->search_all();

for my $step_obj (@step_objs) {
	my $step = $step_obj->name;
	$step_ids{$step} = $step_obj->id;
}

#Get the xml files and create the pipeline object
my $pipe_config = modules::Pipeline::get_pipe_conf();



#First retrieve all the sample groups
my @source_group_objs = modules::Adaptors::Source_Group->search_all();

#Keep track on normal samples that don't have bams yet; these are required for tumour/related snv calling
my %missing_normals = ();
my %missing_unaffected = ();

					
for my $source_group_obj ( @source_group_objs) {
	my $sample_total_field = $source_group_obj->total_samples;
	my $source_group_name = $source_group_obj->source_group_name;
	
	my $source_type = modules::Pipeline::get_source_type(-source_group_name=>$source_group_name);
	
	
	
	#Get the steps for tracking progress
	my $paired_variants = $pipe_config->read($source_type,'paired_variants');
	my $parallel_vcf = $pipe_config->read($source_type,'parallel_vcf');

	my $last_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'last');
	my $last_align_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'finish_align');
	my $parvar_first = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_first');
	my $parvar_record = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_record');
	my $parvar_merge = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'parvar_merge');
	my $nonpar_start = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'start_non_parvar');
	
	
	if (!$paired_variants && !$parallel_vcf) {
		#Simplest case where pipeline is basically one block; human_single, human_related and mouse_single source group
		my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);	
		
		
	
		
		for my $sample_obj ( @samples ) {
			my $sample_name = $sample_obj->sample_name;
			
			next if exists $sample_missing_lanes{$sample_name};
			
			my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
			my ($source_name) = modules::Pipeline::get_source_name($sample_name);
			my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/runs';
			
			
			#Get the latest run if any
	    	my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_obj->id);
	    
	    	#Check if last_step has been run
	    	if (defined $run_obj) {
	    		my $run_id = $run_obj->id;
	    		my @last_step = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_pipe_step});	    
		    	if (!@last_step) {
		    		#Haven't completed the latest step
		    		my $last_step_run = modules::Pipeline::get_latest_step($run_id);
					print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
					
					#Here we get the sample xml in case we need to resubmit jobs    
				    if ($restart && !$job_id) {
						my $sample_xml = $qsub_dir .'/'. $sample_name . '.xml';
						my $next_step = modules::Pipeline::get_next_step($run_id,$sample_xml,$pipe_config,$source_type,$source_type);
						my $qsub_file = $sample_name . '.pipe3.'.$next_step .'.qsub';
						push @commands, "qsub $qsub_dir/$qsub_file";
				    }
					
		    	}
			} else {
					
				my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
						
				if ($job_id) {
					print STDERR "Sample $sample_name is ready to create run but $job_id is queued\n";
				} else {
					print STDERR "Sample group $source_group_name sample $sample_name has no runs created..\n";
				}
			} 
		}
	} elsif ($parallel_vcf && !$paired_variants) {
		my %chromosomes = map {$_ => 1} split(" ",$pipe_config->read($source_type,'annotation_version','chr'));
		my $chr_count = keys %chromosomes;
		my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);	
		for my $sample_obj ( @samples ) {
			my $sample_name = $sample_obj->sample_name;
			next if exists $sample_missing_lanes{$sample_name};
			#Get the latest run if any
	    	my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_obj->id);
	    	my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
	    	
	    	my $db_total_lanes = $sample_obj->total_lanes;
	    	my ($source_name) = modules::Pipeline::get_source_name($sample_name);
	    	my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/runs';
	    	
	    
	    	#Check if last_step has been run
	    	if (defined $run_obj) {
	    		my $run_id = $run_obj->id;
	    		my @bamstats_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_align_step});	
		    	if (@bamstats_steps == 1) {
		    		#Now check if submit_snvs step has been run
					my @callsnv_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$parvar_first});
					
					#If snvs not submitted yet
					if (!@callsnv_steps) {
						
						if ($job_id) {
							print STDERR "Sample $sample_name is ready to submit snvs but $job_id is queued\n";
						} else {
							print STDERR "Sample $sample_name can submit snvs with $db_total_lanes lanes..\n";
							my $qsub_file = $sample_name . '.pipe2.'.$parvar_first.'.qsub'; 
							my $qsub_snvcalls = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/snvcalls';
							push @commands, "qsub $qsub_dir/$qsub_file";	
						}
					} else {
						#Here jobs have been submitted so check if there is a complete run
						my @mergevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$parvar_merge});
						
						
						if (@mergevcf_steps) {
							#Get the release file and check the total lane numbers
							my @total_lanes = modules::Adaptors::Release_File->search_total_lanes($run_id,$step_ids{$parvar_merge});
							if (! @total_lanes) {
								modules::Exception->throw("ERROR: Sample group $source_group_name sample $sample_name has run $parvar_merge but no release file");
							} elsif (@total_lanes > 1) {
								modules::Exception->throw("ERROR: Sample group $source_group_name sample $sample_name has multiple release files for $parvar_merge");
							}
							
							my $release_file_lane_total = $total_lanes[0]->total_lanes;
												
							#If the lane number is different then submit jobs 
							if ($release_file_lane_total != $db_total_lanes) {
								my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
								
								if ($job_id) {
									print STDERR "Sample $sample_name is ready to submit snvs but $job_id is queued\n";
								} else {
									print STDERR "Sample $sample_name can submit snvs as it has different lane number ($db_total_lanes) than existing run ($release_file_lane_total) for run $run_id\n";
									my $qsub_file = $sample_name . '.pipe2.'.$parvar_first.'.qsub'; 
									my $qsub_snvcalls = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/snvcalls';
									push @commands, "qsub $qsub_dir/$qsub_file";	
								}
							} else {
								#Check the last step is run
								my ($last_step_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_pipe_step});
								if (defined $last_step_obj) {
									#print STDERR "\tSample $sample_name with $db_total_lanes lanes has finished step and is complete..\n";					
								} else {
									#Last step not complete; find the latest step run
									my $last_step_run = modules::Pipeline::get_latest_step($run_id);
									print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
									
									#Here we get the sample xml in case we need to resubmit jobs    
								    if ($restart && !$job_id) {
										my $sample_xml = $qsub_dir .'/'. $sample_name . '.xml';
										my $next_step = modules::Pipeline::get_next_step($run_id,$sample_xml,$pipe_config,$source_type,$source_type);
										my $qsub_file = $sample_name . '.pipe3.'.$next_step .'.qsub';
										push @commands, "qsub $qsub_dir/$qsub_file";
								    }
									
								}
							}
						} else {
							#Check how many (if any) single_vcf jobs are done
							my @singlevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$parvar_record});
							#Get the syscalls for the single_vcfs
							my %chr_done = ();
							
							for my $singlevcf_step_obj ( @singlevcf_steps ) {
							    my ($syscall_obj) = modules::Adaptors::Syscall->search('pipeline_steps_run_id'=>$singlevcf_step_obj->id);
							    my $pipe_step_id = $singlevcf_step_obj->id;
							    
							    if (!defined $syscall_obj) {
	    							#Delete orphaned pipeline_steps_runs entries (no syscall)
	    							modules::Adaptors::Pipeline_Step_Run->search(id=>$pipe_step_id)->delete_all;
	    						} else {
							    	$chr_done{$syscall_obj->chr} = 1 if defined $syscall_obj;
	    						}
							}
							
							if (keys %chr_done == $chr_count) {
								#If run for each chr we start running pipe3 from merge_vcf to the end unless job already started
								my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/runs';
								my $qsub_file = $sample_name . '.pipe3.'.$parvar_merge.'.qsub';
								 
								#Check the job hasn't been started yet
								my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
								
								if ($job_id) {
									print STDERR "Sample $sample_name has completed all 24 $parvar_record but $job_id has started running..\n";
								} else {
									print STDERR "Sample $sample_name has completed all 24 $parvar_record so begin pipe3 ..\n";											
									push @commands, "qsub $qsub_dir/$qsub_file";
								}
								
							} else {
								my @missing_chr = ();
								for my $chr (keys %chromosomes) {
									if (!exists $chr_done{$chr}) {
										push @missing_chr, $chr;
									}
								}
								my $chr_missing = join(",",sort @missing_chr);
								my $missing_count = @missing_chr;
								print STDERR "Sample $sample_name has submitted jobs but $missing_count chromosome are not complete yet ($chr_missing)..\n";	
																
							}	
						} 
					}
		    	} elsif (@bamstats_steps == 0) {
		    		my $last_step_run = modules::Pipeline::get_latest_step($run_id);
					print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
					#Here we get the sample xml in case we need to resubmit jobs    
				    if ($restart && !$job_id) {
						my $sample_xml = $qsub_dir .'/'. $sample_name . '.xml';
						my $next_step = modules::Pipeline::get_next_step($run_id,$sample_xml,$pipe_config,$source_type);
						my $qsub_file = $sample_name . '.pipe1.'.$next_step .'.qsub';
						push @commands, "qsub $qsub_dir/$qsub_file";
				    }
		    	} else {
		    		modules::Exception->throw("ERROR: Sample $sample_name has run $last_align_step multiple times");
		    	}
			} else {
					
				my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
						
				if ($job_id) {
					print STDERR "Sample $sample_name is ready to create run but $job_id is queued\n";
				} else {
					print STDERR "Sample group $source_group_name sample $sample_name has no runs created..\n";
				}
			} 
		}
	} elsif ($paired_variants && $parallel_vcf) {
		#Tumour/normal pairs -> parallel and paired
		
		my %chromosomes = map {$_ => 1} split(" ",$pipe_config->read($source_type,'annotation_version','chr'));
		
		
		my @normal_samples = my @tumour_samples = ();
		
		my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);
		
		for my $sample_obj ( @samples ) {
		    if (modules::Pipeline::get_tumour_flag(-sample_name=>$sample_obj->sample_name)) {
		    	push @tumour_samples, $sample_obj;
		    } else {
		    	push @normal_samples, $sample_obj;
		    }
		}
		
		if (!@normal_samples) {
			modules::Exception->throw("ERROR: Sample group $source_group_name has no normal sample");
		} elsif (@normal_samples > 1) {
			modules::Exception->throw("ERROR: Sample group $source_group_name has >1 normal samples");
		} 
			    
			    
		#my $normal_sample_id = $normal_samples[0]->id;
		
		my $normal_sample_name = $normal_samples[0]->sample_name;
		my $normal_total_lanes = $normal_samples[0]->total_lanes;
		my ($source_name) = modules::Pipeline::get_source_name($normal_sample_name);
	
		#Now check if the tumour sample is ready for snv calling 
	
		#The +1 is for the normal sample
		my $sample_total_db = 1 + @tumour_samples; 
	
	
		#Check the sample numbers agree
		if ($sample_total_field != $sample_total_db) {
			modules::Exception->throw("ERROR: Sample group $source_group_name expects $sample_total_field samples and there are $sample_total_db samples in the database\n");
		}
		
		my $sample_ready_count = 0;
		my $complete = 0;
		
		
		#now check the runs are complete
		for my $sample ( @normal_samples,@tumour_samples ) {
			my $sample_id = $sample->id;
		    my $sample_name = $sample->sample_name;
		    next if exists $sample_missing_lanes{$sample_name};
		    
		    #Check the job hasn't been started yet
			my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
		    my ($source_name) = modules::Pipeline::get_source_name($sample_name);
	    	my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/runs';
			my $db_total_lanes;
			my $tumour_flag = modules::Pipeline::get_tumour_flag(-sample_name=>$sample_name);
			#Tumour uses normal lanes as well
			if ($tumour_flag) {
				$db_total_lanes = $normal_total_lanes  + $sample->total_lanes;
			} else {
				$db_total_lanes = $normal_total_lanes;
			}
			
		    
		    #Get the latest run if any
		    my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_id);
		    
		    #Check if bam_stats has been run
		    if (defined $run_obj) {
		    	my $run_id = $run_obj->id;
		    	my @bamstats_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_align_step});	
		    	if (@bamstats_steps == 1) {
		    		#Don't proceed past here with normal samples...
		    		if (!$tumour_flag) {
			    		$complete++;
			    		next;				
		    		}
		    		
		    		if ($source_type eq 'mouse_cancer') {
		    			#Simpler case as don't do parallel vcf
		    			if (exists $missing_normals{$source_group_name}) {
							print "\n";
							next;
						}
						
						my ($last_step_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_pipe_step});
						if (defined $last_step_obj) {
							next;
						}
						
						#Ready to resume at call variants
			    		my ($source_name) = modules::Pipeline::get_source_name($sample_name);

						
			    		
						my $qsub_file = $sample_name . '.pipe2.'.$nonpar_start.'.qsub';
						my $last_step_run = modules::Pipeline::get_latest_step($run_id);
                        if ($last_step_run ne $last_align_step) {
                        	#Here it in process in the last block
                        	print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
                        	#Here we get the sample xml in case we need to resubmit jobs    
						    if ($restart && !$job_id) {
								my $sample_xml = $qsub_dir .'/'. $sample_name . '.xml';
								my $next_step = modules::Pipeline::get_next_step($run_id,$sample_xml,$pipe_config,$source_type);
								my $qsub_file = $sample_name . '.pipe3.'.$next_step .'.qsub';
								push @commands, "qsub $qsub_dir/$qsub_file";
						    }
                        	next;
                        }
									 
						
						if ($job_id) {
							print STDERR "Sample $sample_name is ready to begin pipe2 but $job_id has started running..\n";
						} else {
							print STDERR "Sample $sample_name is ready to begin pipe2 ..\n";											
							push @commands, "qsub $qsub_dir/$qsub_file";
						}
						
		    			next;
		    		}
		    		
		    		
					#Now check if submit_snvs step has been run
					my @callsnv_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$parvar_first});
					
					#If snvs not submitted yet
					if (!@callsnv_steps) {
						#Tumour snv calling requires the normal bam file so skip if it doesn't exist
						if (exists $missing_normals{$source_group_name}) {
							print "\n";
							next;
						}
						my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
								
						if ($job_id) {
							print STDERR "Sample $sample_name is ready to submit snvs but $job_id is queued\n";
						} else {
							print STDERR "Sample $sample_name can submit snvs with $db_total_lanes lanes..\n";
							my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/runs';
							my $qsub_file = $sample_name . '.pipe2.'.$parvar_first.'.qsub'; 
							my $qsub_snvcalls = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/snvcalls';
							push @commands, "qsub $qsub_dir/$qsub_file";
						}	
					} else {
						#Here jobs have been submitted so check if there is a complete run
						my @mergevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$parvar_merge});
						
						
						if (@mergevcf_steps) {
							#Get the release file and check the total lane numbers
							my @total_lanes = modules::Adaptors::Release_File->search_total_lanes($run_id,$step_ids{$parvar_merge});
							if (! @total_lanes) {
								modules::Exception->throw("ERROR: Sample group $source_group_name sample $sample_name has run $parvar_merge but no release file");
							} elsif (@total_lanes > 1) {
								modules::Exception->throw("ERROR: Sample group $source_group_name sample $sample_name has multiple release files for $parvar_merge");
							}
							
							my $release_file_lane_total = $total_lanes[0]->total_lanes;
												
							#If the lane number is different then submit jobs 
							if ($release_file_lane_total != $db_total_lanes) {
								
								if ($job_id) {
									print STDERR "Sample $sample_name is ready to submit snvs but $job_id is queued\n";
								} else {
									print STDERR "Sample $sample_name can submit snvs as it has different lane number ($db_total_lanes) than existing run ($release_file_lane_total) for run $run_id\n";
									my $qsub_file = $sample_name . '.pipe2.'.$parvar_first.'.qsub'; 
									my $qsub_snvcalls = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/snvcalls';
									push @commands, "qsub $qsub_dir/$qsub_file";	
										
								}
							} else {
								#Check the last step is run
								my ($last_step_obj) = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_pipe_step});
								if (defined $last_step_obj) {
									$complete++;
									#print STDERR "\tSample $sample_name with $db_total_lanes lanes has finished step and is complete..\n";					
								} else {
									#Last step not complete; find the latest step run
									my $last_step_run = modules::Pipeline::get_latest_step($run_id);
									print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
									 if ($restart && !$job_id) {
										my $sample_xml = $qsub_dir .'/'. $sample_name . '.xml';
										my $next_step = modules::Pipeline::get_next_step($run_id,$sample_xml,$pipe_config,$source_type);
										my $qsub_file = $sample_name . '.pipe3.'.$next_step .'.qsub';
										push @commands, "qsub $qsub_dir/$qsub_file";
								    }
								}
							}
						} else {
							#Check how many (if any) single_vcf jobs are done
							my @singlevcf_steps = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$parvar_record});
							
							#Get the syscalls for the single_vcfs
							my %chr_done = ();
							
							for my $singlevcf_step_obj ( @singlevcf_steps ) {
							    my ($syscall_obj) = modules::Adaptors::Syscall->search('pipeline_steps_run_id'=>$singlevcf_step_obj->id);
							    my $pipe_step_id = $singlevcf_step_obj->id;
							    
							    if (!defined $syscall_obj) {
	    							#Delete orphaned pipeline_steps_runs entries (no syscall)
	    							modules::Adaptors::Pipeline_Step_Run->search(id=>$pipe_step_id)->delete_all;
	    						} else {
							    	$chr_done{$syscall_obj->chr} = 1 if defined $syscall_obj;
	    						}
							}
							
							if (keys %chr_done == 24) {
								#If run for each chr we start running pipe3 from merge_vcf to the end unless job already started
								my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/runs';
								my $qsub_file = $sample_name . '.pipe3.'.$parvar_merge.'.qsub';
								 
								#Check the job hasn't been started yet
								my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
								
								if ($job_id) {
									print STDERR "Sample $sample_name has completed all 24 $parvar_record but $job_id has started running..\n";
								} else {
									print STDERR "Sample $sample_name has completed all 24 $parvar_record so begin pipe3 ..\n";											
									push @commands, "qsub $qsub_dir/$qsub_file";
								}
								
							} else {
								my @missing_chr = ();
								for my $chr (keys %chromosomes) {
									if (!exists $chr_done{$chr}) {
										push @missing_chr, $chr;
									}
								}
								my $chr_missing = join(",",sort @missing_chr);
								my $missing_count = @missing_chr;
								print STDERR "Sample $sample_name has submitted jobs but $missing_count chromosome are not complete yet ($chr_missing)..\n";	
																
							}	
						} 
					}
		    	} elsif (@bamstats_steps == 0) {
		    		$missing_normals{$source_group_name}++ if $tumour_flag == 0;
		    		my $last_step_run = modules::Pipeline::get_latest_step($run_id);
					print STDERR "Sample $sample_name has only run up to step $last_step_run...\n";
					if ($restart && !$job_id) {
						my $sample_xml = $qsub_dir .'/'. $sample_name . '.xml';
						my $next_step = modules::Pipeline::get_next_step($run_id,$sample_xml,$pipe_config,$source_type);
						my $qsub_file = $sample_name . '.pipe1.'.$next_step .'.qsub';
						push @commands, "qsub $qsub_dir/$qsub_file";
				    }
		    	} else {
		    		modules::Exception->throw("ERROR: Sample $sample_name has run $last_align_step multiple times");
		    	}
		    } else {
		    	$missing_normals{$source_group_name}++ if $tumour_flag == 0;
		    	#Check the job hasn't been started yet
				my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample_name);
								
				if ($job_id) {
					print STDERR "Sample $sample_name is ready to create a run but $job_id is queued..\n";
				} else {
		    		print STDERR "Sample group $source_group_name sample $sample_name has no runs created..\n";
				}
		    }
		}
		my $total_samples = @normal_samples + @tumour_samples;
	    print "\n" unless $total_samples == $complete;
	} else {
		modules::Exception->throw("ERROR: No handling for sample $source_group_name\n");
	}
}


print STDERR "\nRun the following commands:\n\n" if @commands;
for my $command ( @commands ) {
   print STDERR "$command\n";
   system("$command") if $OPT{submit};
}











