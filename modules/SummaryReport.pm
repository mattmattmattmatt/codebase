package modules::SummaryReport;

use strict;
use Data::Dumper;
use modules::ConfigXML;
use modules::VariantXML;
use modules::Adaptors::Pipeline_Step;
use modules::Utils;
use modules::Pipeline;
use modules::Exception;
use modules::Adaptors::Sample;
use modules::Adaptors::Lane;
use modules::Adaptors::Group_Summary;
use modules::Adaptors::Project_Summary;
use modules::Adaptors::Release_File;
use modules::Adaptors::Human_Related_Sample;
use modules::Adaptors::Human_Single_Sample;
use modules::Adaptors::Human_Cancer_Sample;
use modules::SystemCall;
use File::Basename;

sub new {
	my ($class, @args) = @_;
	
	my $self = bless {}, $class;

    my %args = @args;

	

    my @required_args = (
    					-sources,
    					-snv_gene_col_name,
    					-indel_gene_col_name,
    					-sv_gene_col_name,
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
   
   	$self->{sources} = $args{-sources};
   	$self->{snv_gene_col_name} = $args{-snv_gene_col_name};
   	$self->{indel_gene_col_name} = $args{-indel_gene_col_name};
   	$self->{sv_gene_col_name} = $args{-sv_gene_col_name};
   
   	if (exists $args{-production} && $args{-production}) {
   		$self->{production} = 1;
   	} else {
   		$self->{production} = 0;
   	}
   	
   	if (exists $args{-fail} && $args{-fail}) {
   		$self->{fail} = 1;
   	} else {
   		$self->{fail} = 0;
   	}
   	
   	if (exists $args{-family} && $args{-family}) {
   		$self->{family} = 1;
   	} else {
   		$self->{family} = 0;
   	}
   	
   	if (exists $args{-report_type}) {
   		$self->{report_type} = $args{-report_type};
   	} else {
   		$self->{report_type} = 'all_db';
   	}
   	
   	my $filter_output = 0;
   	
   	if (exists $args{-chrom}) {
   		$filter_output = 1;
   		$self->{filters}{chr} = $args{-chrom};
   	}
   	
   	if (exists $args{-coord}) {
   		$filter_output = 1;
   		my ($chr,$start,$end) = $args{-coord} =~ /([0-9XYM]+):(\d+)\-(\d+)/;
   		$self->{filters}{chr} = $chr; #overwrite $args{-chr} if also passed in
   		$self->{filters}{start} = $start;
   		$self->{filters}{end} = $end;
   	}
   	
   	if (exists $args{-gene_list}) {
   		$filter_output = 1;
   		$self->{filters}{gene_list} = $args{-gene_list};
   	}
   	
   	if (exists $args{-pileups} && $args{-pileups}) {
   		$self->{pileups} = 1;
   	} else {
   		$self->{pileups} = 0;
   	}

   	if (exists $args{-denovo}) {
   		$filter_output = 1;
   		$self->{filters}{denovo} = 1;
   	} 
   	
   	if (exists $args{-comhet}) {
   		$filter_output = 1;
   		$self->{filters}{comhet} = 1;
   	}
   	
   	if (exists $args{-phase_block_size}) {
   		$filter_output = 1;
   		$self->{filters}{phase_block_size} = $args{-phase_block_size};
   	} 
   	
   	if (exists $args{-phase_var_num}) {
   		$filter_output = 1;
   		$self->{filters}{phase_var_num} = $args{-phase_var_num};
   	} 
   	
   	if (exists $args{-inheritance}) {
   		$filter_output = 1;
   		$self->{filters}{inheritance} = $args{-inheritance};
   	} 
   	
   	if (exists $args{-max_allele_freq}) {
   		$filter_output = 1;
   		$self->{filters}{max_allele_freq} = $args{-max_allele_freq};
   	} 
   	
   	if (exists $args{-min_num_aff}) {
   		$filter_output = 1;
   		$self->{filters}{min_num_aff} = $args{-min_num_aff};
   	} 
   	
   	if (exists $args{-comhet}) {
   		$filter_output = 1;
   		$self->{filters}{comhet} = 1;
   	}
   	
   	if (exists $args{-sift}) {
   		$filter_output = 1;
   		$self->{filters}{sift} = $args{-sift};
   	} 
   	
   	if (exists $args{-polyphen}) {
   		$filter_output = 1;
   		$self->{filters}{polyphen} = $args{-polyphen};
   	} 
   	
   	if (exists $args{-min_num_aff}) {
   		$filter_output = 1;
   		$self->{filters}{min_num_aff} = $args{-min_num_aff};
   	}
   	
   	if ($filter_output) {
   		$self->{filter_output} = 1;
   	}

	if(exists $args{-no_chr}){
		$self->{no_chr}=1;
	}else{
		$self->{no_chr}=0;
	}
   	
   	if(exists $args{-sv}){
		$self->{sv} = 1;
	}else{
		$self->{sv} = 0;
	}
   	
   	if(exists $args{-score} && $args{-score}) {
   		$self->{score} = 1;
   	} else {
   		$self->{score} = 0;
   	}
   	
   	if (exists $args{-force}) {
   		$self->{force} = 1;
   	} else {
   		$self->{force} = 0;
   	}
   	
   	if(exists $args{-report_xml}) {
   		$self->{report_xml} = $args{-report_xml};
   	}
   	
   	print Dumper $self;
   	
    return $self;
}


#Loads all the sample info into the object
sub set_sample_info {
	my ($self,$sample_included) = @_;
	
	my %included_samples = ();
	if ($sample_included ne 'all') {
		%included_samples = %{$sample_included};
	}
	
	my $samples_to_process = 0;
	
	#Get all the step mappings from the database
	my %step_ids = ();
	
	my @step_objs = modules::Adaptors::Pipeline_Step->search_all();
	
	for my $step_obj (@step_objs) {
		my $step = $step_obj->name;
	}
	
	my $project_sample_count;
	
	#Get project sample count here to avoid recalculating
	if ($self->{report_type} eq 'project') {
		for my $source (@{$self->{sources}}) {
			my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source);
			my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
			
			for my $source_group_obj (@source_group_objs) {
				my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);
				$project_sample_count += @samples;
			}
		}
	}
	
	my $clus_conf = modules::Pipeline::get_cluster_conf();
	my ($date_stamp) = split('_',modules::Utils->GetTime());
	
	#Get the report info from the release files in the database
	for my $source (@{$self->{sources}}) {
		
		my $source_type = modules::Pipeline::get_source_type(-source_name=>$source);
		
		#Now check all the samples have run to completion
		my $last_pipe_step = modules::Pipeline::get_xml_step(-source_type=>$source_type,-xml_name=>'last');
		my $rundir = $clus_conf->read($source_type,'base_directories','base_run_directory');
		my $resultdir = $clus_conf->read($source_type,'base_directories','base_results_directory');
		
		my $basedir = $rundir . '/' . $source;
		
		if (!-d $basedir) {
			modules::Exception->throw("ERROR: base directory $basedir for source $source doesn't exist");
		}
		
		my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source);
		if (! defined $source_obj ) {
			#next;
			modules::Exception->throw("ERROR: Can't retrieve source with name $source");
		}
		my @source_group_objs = modules::Adaptors::Source_Group->search(source_id=>$source_obj->id);
		
		if (! @source_group_objs ) {
			modules::Exception->throw("ERROR: Can't retrieve source_group with source id $source_obj->id");
		}
		
		for my $source_group_obj (@source_group_objs) {
			my $project_name;
			my $project_name_nospace;
			my $sg_name; 
			my $sample_summary_number;
			


			if ($self->{report_type} eq 'project') {
				$project_name = modules::Pipeline::get_project_name(-source_name=>$source);
				$self->{project_name} = $project_name;
				($project_name_nospace = $project_name) =~ s/ /_/g;
				$sg_name = $project_name_nospace;
			} elsif ($self->{report_type} eq 'sample') {
				$sg_name = 'combined_samplelist';				
				#Check if there are any matching samples in the source
				my $match_sample = 0;
				my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);	
				for my $sample_obj ( @samples ) {
					my $sample_name = $sample_obj->sample_name;
					# "include" field per sample rather than using $match_sample
					$self->{source_info}{$sg_name}{samples}{$sample_name}{include} = 0;
				
					#Check if in -sample_list file
#					$match_sample = 1 if exists $included_samples{$sample_name};
					$self->{source_info}{$sg_name}{samples}{$sample_name}{include} = 1 if exists $included_samples{$sample_name};
					next unless $self->{source_info}{$sg_name}{samples}{$sample_name}{include};
				}

				$sample_summary_number= keys %included_samples;
				#next unless $match_sample;
			} elsif ($self->{report_type} eq 'all_db') {
				$sg_name = 'alldb';
			} else {
				$sg_name = $source_group_obj->source_group_name;
				$sample_summary_number = $source_group_obj->total_samples;
			}
			
			my $outdir;
			my $coverage_file;
			
			if ($self->{production}) {
				my $prod_entry_found = 0;
				if ($self->{report_type} eq 'project') {
					my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
					my (undef,$base_results) = fileparse($resultdir);
					$outdir = $base_results . '/project_summaries/' . $project_name_nospace; 
					
					mkdir($outdir) if !-d $outdir;
					
					my @project_summary_objs = modules::Adaptors::Project_Summary->search(project_id=>$project_obj->id);
					if (@project_summary_objs) {
						#Here there are recorded summaries
						for my $proj_summary_obj (@project_summary_objs) {
							my $proj_summary_obj_count = $proj_summary_obj->total_samples;
							
							$coverage_file = $outdir . '/' . $project_name_nospace ."_" . $project_sample_count .'_samples_cover_' .$date_stamp.'.tsv';

							#Skip check if we've run earlier iterations with less samples
							next unless $proj_summary_obj_count == $project_sample_count;
							
							if (!-d $outdir) {
								modules::Exception->throw("ERROR: Summary project obj exists but output directory $outdir doesn't");
							}
							
							if (-e $coverage_file) {
								print "Skip project $project_name for cover file $coverage_file\n";
								$prod_entry_found = 1;
							}
						}
					} else {
						$coverage_file = $outdir . '/' . $project_name_nospace ."_" . $project_sample_count .'_samples_cover_' .$date_stamp.'.tsv';
					}
				} elsif ($self->{family}) { #either report_type = source or -all_related
					$outdir = $resultdir . '/' . $source . '/source_group_summaries'; #Standard production location
					#check the entry doesn't exist
					my @summary_group_objs = modules::Adaptors::Group_Summary->search(source_group_id=>$source_group_obj->id);
					if (@summary_group_objs) {
						#Here there are recorded summaries
						for my $summary_group_obj (@summary_group_objs) {
							if (!-d $outdir) {
								modules::Exception->throw("ERROR: Summary group obj exists but output directory $outdir doesn't");
							}
							$sample_summary_number = $summary_group_obj->total_samples;
							$coverage_file = $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_cover_' .$date_stamp.'.tsv';	
							if (-e $coverage_file) {
								print "Skip source $source source_group $source_group_obj for $sg_name (id = $sample_summary_number) samples $summary_group_obj\n";
								$prod_entry_found = 1;
							}
						}
					} else {
						$coverage_file = $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_cover_' .$date_stamp.'.tsv';
					}
				} else {
					#Everything else in cwd
					$outdir = './';
				}
				
				if ($prod_entry_found && !$self->{sv}) {
					next;
				}
				mkdir($outdir) if !-d $outdir;
				
				#Check if analysis has been done
				if ( -e $coverage_file && $self->{report_type} ne 'project' && !$self->{sv}) { #Don't do coverage for sv so skip this bail out
					modules::Exception->warning("Skip for coverage_file $coverage_file as it already exists");	
				}
			
				
			} else {
				$outdir = './'; #non production write to current directory							
				
				if ($self->{report_type} eq 'project') {
					$coverage_file = $outdir . '/' . $project_name_nospace ."_" . $project_sample_count .'_samples_cover_' .$date_stamp.'.tsv';
				} elsif ($self->{report_type} eq 'all_db') {
					$coverage_file = $outdir . '/' . $sg_name .'_cover_' .$date_stamp.'.tsv';
				} else {
					$coverage_file = $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_cover_' .$date_stamp.'.tsv';
				}
			}
			
			
			
			
			if ($self->{family}) {
				#For cohorts check all samples within the cohort are complete
				my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);	
				
				my $complete_sample_count = 0;
				
				for my $sample_obj ( @samples ) {
					my $sample_name = $sample_obj->sample_name;
					
					if (keys %included_samples) {
						#Check if in -sample_list file
						next unless exists $included_samples{$sample_name};
					}
					
					
					#Get the latest run if any
		    		my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_obj->id);
		    
		    		#Check if last_step has been run
		    		if (defined $run_obj) {
		    			my $run_id = $run_obj->id;
		    			my @last_step = modules::Adaptors::Pipeline_Step_Run->search(run_id=>$run_id,pipeline_step_id=>$step_ids{$last_pipe_step});	    
			    		if (!@last_step) {
			    			next;
			    		}
		    		} else {
		    			next;
		    		}
		    		$complete_sample_count++;
				}			
				
				if ($complete_sample_count != $sample_summary_number) {
					print "Skip summary for $sg_name; not all samples are complete ($complete_sample_count != $sample_summary_number)\n";
					#exit;
				}
			
			} 
			
			#If we make it here then we need to run the summary
			$self->{source_info}{$sg_name}{dir_base} =  $outdir;
		
			
				
			if ($self->{report_type} eq 'project') {
				open(COVER,">>$coverage_file") || modules::Exception->throw("Can't open file to write $coverage_file\n");
				if ($self->{sv}) {
					$self->{source_info}{$sg_name}{sv_file} =   $outdir . '/' . $project_name_nospace ."_" . $project_sample_count.'_samples_sv_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{sv_xml} =   $outdir . '/' . $project_name_nospace ."_" . $project_sample_count.'_samples_sv_' .$date_stamp.'.xml';
				} else {
					$self->{source_info}{$sg_name}{snv_file} =   $outdir . '/' . $project_name_nospace ."_" . $project_sample_count.'_samples_snv_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{snv_xml} =   $outdir . '/' . $project_name_nospace ."_" . $project_sample_count.'_samples_snv_' .$date_stamp.'.xml';
					$self->{source_info}{$sg_name}{indel_file} =   $outdir . '/' . $project_name_nospace ."_" . $project_sample_count.'_samples_indel_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{indel_xml} =   $outdir . '/' . $project_name_nospace ."_" . $project_sample_count.'_samples_indel_' .$date_stamp.'.xml';
				}
				$self->{project_info}{sample_number}++;	
			} elsif ($self->{report_type} eq 'all_db') {
				open(COVER,">>$coverage_file") || modules::Exception->throw("Can't open file to write $coverage_file\n");
				if ($self->{sv}) {
					$self->{source_info}{$sg_name}{sv_file} =   $outdir . '/' . $sg_name . '_sv_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{sv_xml} =   $outdir . '/' . $sg_name . '_sv_' .$date_stamp.'.xml';
				} else {
					$self->{source_info}{$sg_name}{snv_file} =   $outdir . '/' . $sg_name . '_snv_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{snv_xml} =   $outdir . '/' . $sg_name . '_snv_' .$date_stamp.'.xml';
					$self->{source_info}{$sg_name}{indel_file} =   $outdir . '/' . $sg_name . '_indel_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{indel_xml} =   $outdir . '/' . $sg_name . '_indel_' .$date_stamp.'.xml';
				}
				$self->{alldb_info}{sample_number}++;	
			} else {
				open(COVER,">$coverage_file") || modules::Exception->throw("Can't open file to write $coverage_file\n");
				if ($self->{sv}) {
					$self->{source_info}{$sg_name}{sv_file} =   $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_sv_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{sv_xml} =   $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_sv_' .$date_stamp.'.xml';
				} else {
					$self->{source_info}{$sg_name}{snv_file} =   $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_snv_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{snv_xml} =   $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_snv_' .$date_stamp.'.xml';
					$self->{source_info}{$sg_name}{indel_file} =   $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_indel_' .$date_stamp.'.tsv';
					$self->{source_info}{$sg_name}{indel_xml} =   $outdir . '/' . $sg_name ."_" . $sample_summary_number .'_samples_indel_' .$date_stamp.'.xml';
				}
				$self->{source_info}{$sg_name}{sample_number} = $sample_summary_number;
			}
			
			my @sample_objs = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);
			if (! @sample_objs ) {
				modules::Exception->throw("ERROR: Can't retrieve sample with source_group id $source_group_obj->id");
			}
			for my $sample_obj (@sample_objs) {
				
				my $sample_name = $sample_obj->sample_name;
				
				## skip not included samples if report type is 'sample'
				if ($self->{report_type} eq 'sample' && !$self->{source_info}{$sg_name}{samples}{$sample_name}{include}){
                	next;
              	}
				
				#Get the sequencing centre from a lane obj
				my ($lane_obj) = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);
				my $seq_centre_name = $lane_obj->sequencing_centre;
				
	
				my @release_file_objs = modules::Adaptors::Release_File->search_release_files_sample($sample_obj->id);
				if (! @release_file_objs ) {
					print "Skip sample $sample_name; No release files\n";
					next;
				}
	
	
				my $human_sample_obj;
				my $affected;
							
				$samples_to_process = 1; #set the flag that there are sample to process
				#Set the relation and parent_count fields
				if ($source_type =~ /human_single/) {
					($human_sample_obj) = modules::Adaptors::Human_Single_Sample->search(sample_id=>$sample_obj->id);
					if (!defined $human_sample_obj) {
						modules::Exception->throw("ERROR: Can't get human sample object for $sample_name");
					}
					$self->{source_info}{$sg_name}{related} = 0;
					$self->{source_info}{$sg_name}{samples}{$sample_name}{relation} = 0;
					$self->{source_info}{$sg_name}{parent_count} = 0;
					$affected = $human_sample_obj->affected?1:0;
				} elsif ($source_type =~ /human_related/) {
					($human_sample_obj) = modules::Adaptors::Human_Related_Sample->search(sample_id=>$sample_obj->id);
					if (!defined $human_sample_obj) {
						modules::Exception->throw("ERROR: Can't get human sample object for $sample_name");
					}					
					$self->{source_info}{$sg_name}{related} = 1;
					$self->{source_info}{$sg_name}{samples}{$sample_name}{relation} = $human_sample_obj->relation;
					if ($human_sample_obj->relation eq 'mother') {
						$self->{source_info}{$sg_name}{parent_count}++;
						$self->{source_info}{$sg_name}{mother} = 1;
					} elsif ($human_sample_obj->relation eq 'father') {
						$self->{source_info}{$sg_name}{parent_count}++;
						$self->{source_info}{$sg_name}{father} = 1;
					}	
					$affected = $human_sample_obj->affected?1:0;
					
				} elsif ($source_type =~ /human_cancer/){
					($human_sample_obj) = modules::Adaptors::Human_Cancer_Sample->search(sample_id=>$sample_obj->id);
					if (!defined $human_sample_obj) {
						modules::Exception->throw("ERROR: Can't get human sample object for $sample_name");
					}
					if(!$human_sample_obj->tumour){
						next;
					}
					$self->{source_info}{$sg_name}{related} = 0;
					$self->{source_info}{$sg_name}{samples}{$sample_name}{relation} = 0;
					$self->{source_info}{$sg_name}{parent_count} = 0;					
					$affected = 1;
				}
				
				if ($affected) {
					$self->{source_info}{$sg_name}{aff_count}++;
					$self->{source_info}{$sg_name}{affected} = 1;
				} else {
					$self->{source_info}{$sg_name}{unaff_count}++;
					$self->{source_info}{$sg_name}{unaffected} = 1; #Flag for whether to report disease inheritance
				}
				
				$self->{source_info}{$sg_name}{samples}{$sample_name}{affected} = $affected;
				
				for my $release_file_obj (@release_file_objs) {
					my ($sample_and_run) = split('\.',$release_file_obj->file_name);
					my $local_rundir = $rundir  . '/' . $source . '/' . $sample_name . '_runs/' . $sample_and_run;
					
					if ($release_file_obj->file_name =~ /combined/) {
						my $full_file = $resultdir . '/' . $source . '/' . $release_file_obj->file_name;
						if ( !-e $full_file && $full_file =~ /sv/) {
							modules::Exception->throw("File $full_file doesn't exist for $sample_name");	
						}
						$self->{source_info}{$sg_name}{samples}{$sample_name}{sv_file} =  $full_file;	
					} elsif ($release_file_obj->file_name =~ /summary.tsv$/) {
						$self->{source_info}{$sg_name}{samples}{$sample_name}{type} = $sample_obj->sample_type;
						my $full_file = $resultdir . '/' . $source . '/' . $release_file_obj->file_name;
						if ( !-e $full_file && ($full_file =~ /snv/ || $full_file =~ /indel/) ) {
							modules::Exception->throw("File $full_file doesn't exist for $sample_name");	
						}
						if ($full_file =~ /snv/) {
							$self->{source_info}{$sg_name}{samples}{$sample_name}{snv_file} =  $full_file;	
						} if ($full_file =~ /indel/) {
							$self->{source_info}{$sg_name}{samples}{$sample_name}{indel_file} =  $full_file;	
						}
					} elsif ($release_file_obj->file_name =~  /bam$/) {
						next if $self->{sv};
						my $full_bam_file = $local_rundir .'/bam/' . $release_file_obj->file_name;
						$full_bam_file =~ s/merge_bam.out/gatk/;
#						print "bam file found $full_bam_file for sample $sample_name\n";
						$self->{source_info}{$sg_name}{samples}{$sample_name}{bam_file} =  $full_bam_file;	
						if (!-e $full_bam_file) {
							modules::Exception->throw("ERROR: Can't access bam file $full_bam_file for $sample_name");
						}
#						next;
#						my $full_bam_file = $basedir . '/bam_links/' . $sample_name . '.bam';
#						$self->{source_info}{$sg_name}{samples}{$sample_name}{bam_file} =  $full_bam_file;	
#						
#						if (!-e $full_bam_file) {
#							modules::Exception->throw("ERROR: Can't access bam file $full_bam_file for $sample_name");
#						}
						
						#Get the bam_stats file for coverage stats
						my $read_report_file = $local_rundir . '/summary/' . $sample_and_run . '.readReport.summary.txt';
						if ( !-e $read_report_file ) {
							modules::Exception->throw("File $read_report_file doesn't exist");	
						}
						open(FILE,$read_report_file) || modules::Exception->throw("ERROR: Can't open $read_report_file");
						while (<FILE>) {
							if (/^aligned: (\d+)/) {
								
								if ($sample_obj->sequence_type eq 'genome') {
									my $coverage = sprintf ("%.2f",100 * $1 / 2897310462);
									print COVER join("\t",
													$sample_name,
													$seq_centre_name,
													$coverage . 'X'
												) . "\n";							
								}
								
							} elsif (/^exon_aligned: (\d+)/) {
								if ($sample_obj->sequence_type eq 'exome') {
									my $coverage = sprintf ("%.2f",100 * $1 / 33438284);
									print COVER join("\t",
												$sample_name,
												$seq_centre_name,
												$coverage . 'X'
											) . "\n";							
								}
							}
						}
						close FILE;		 			
						
					} 
				}
			}
		close COVER;
		}
	}
#	print Dumper $self;
	return $samples_to_process;
}

#Set the master headers from the first report files for an affected for each cohort
sub set_master_headers {
	my ($self) = @_;
	my %snv_headers = ();
	my %indel_headers = ();
	my %sv_headers = ();
	
	my %sample_info = %{$self->{source_info}};
	#print Dumper \%sample_info;
	#Get the master headers index to set for all files; use affected here as they have a few different headers
	for my $sg_name ( keys %sample_info ) {
		for my $sample ( keys %{$sample_info{$sg_name}{samples}} ) {
			
			## skip not included samples if report type is 'sample'
			if ($self->{report_type} eq 'sample' && !$self->{source_info}{$sg_name}{samples}{$sample}{include}){
				next;
			}
		
			if ($self->{sv}) {
				if (!exists $sample_info{$sg_name}{samples}{$sample}{sv_file}) {
					modules::Exception->throw("ERROR: sv file for $sample doesn't exist");
				}
				
				
				#next unless $sample_info{$sg_name}{samples}{$sample}{affected} == 1;
				open(SV_REPORT,$sample_info{$sg_name}{samples}{$sample}{sv_file}) || modules::Exception->throw("ERROR: Can't open snv file $sg_name $sample $sample_info{$sg_name}{samples}{$sample}{sv_file}");	
				
				while (<SV_REPORT>) {
					chomp;
					if (/^SV/) {
						my @sv_report_headers = split("\t");
						
						#Set the master ordering from the first file we see
						for (my $i = 0; $i < scalar @sv_report_headers; $i++){
		    				$sv_headers{$sv_report_headers[$i]} = $i;
						}
						
					
					}
				}
				
			} else {
				if (!exists $sample_info{$sg_name}{samples}{$sample}{snv_file}) {
					modules::Exception->throw("ERROR: snv file for $sample doesn't exist");
				}
				
				#next unless $sample_info{$sg_name}{samples}{$sample}{affected} == 1;
				open(SNV_REPORT,$sample_info{$sg_name}{samples}{$sample}{snv_file}) || modules::Exception->throw("ERROR: Can't open snv file $sg_name $sample $sample_info{$sg_name}{samples}{$sample}{snv_file}");	
				
				while (<SNV_REPORT>) {
					chomp;
					if (/^chr/) {
						my @snv_report_headers = split("\t");
						
						#Set the master ordering from the first file we see
						for (my $i = 0; $i < scalar @snv_report_headers; $i++){
		    				$snv_headers{$snv_report_headers[$i]} = $i;
						}
						
					
					}
				}
				
				if (!exists $sample_info{$sg_name}{samples}{$sample}{indel_file}) {
					modules::Exception->throw("ERROR: indel file for $sample doesn't exist");
				}
				
				open(INDEL_REPORT,$sample_info{$sg_name}{samples}{$sample}{indel_file}) || modules::Exception->throw("ERROR: Can't open indel file $sample_info{$sg_name}{samples}{$sample}{indel_file}");
				while (<INDEL_REPORT>) {
					chomp;
					if (/^chr/) {
						my @indel_report_headers = split("\t");
						
						for (my $i = 0; $i < scalar @indel_report_headers; $i++){
		    				$indel_headers{$indel_report_headers[$i]} = $i;
						}
					}
				}
				
			}
		
		
		}
		if ($self->{sv}) {
			$self->{source_info}{$sg_name}{sv_headers} = \%sv_headers;
		} else {
			$self->{source_info}{$sg_name}{snv_headers} = \%snv_headers;
			$self->{source_info}{$sg_name}{indel_headers} = \%indel_headers;
			
		}
	}
	
	#Contains the pileup info
	my @extra_headers = ();
	push @extra_headers, 'samples';
	
	#Now set the additional headers
	if ($self->{family} && !$self->{sv}) {
		my $mother_in_any = 0;
		my $father_in_any = 0;
		for my $sg_name ( keys %sample_info ) {
			if ($self->{source_info}{$sg_name}{affected}) {
				push @extra_headers, 'number_affected_variant/total_affected';
			}
			if ($self->{source_info}{$sg_name}{unaffected}) {
				push @extra_headers, 'number_unaffected_variant/total_unaffected';	
			}
			
			if ($self->_has_two_parents($sg_name)) {
				push @extra_headers, 'mendelian_inheritance','disease_inheritance','gene','total_variants_in_gene','unique_coord_gene_variants','%_coding_bases_variant','mother_allele','father_allele','parent_allele_common_to_affected','affected_allele_block','definite_compound_het','possible_compound_het','rare_definite_compound_het','rare_possible_compound_het';
			} elsif ($self->_has_mother($sg_name)) {
				push @extra_headers,  'mendelian_inheritance','disease_inheritance','gene','total_variants_in_gene','unique_coord_gene_variants','%_coding_bases_variant','mother_allele','parent_allele_common_to_affected','affected_allele_block','definite_compound_het','possible_compound_het','rare_definite_compound_het','rare_possible_compound_het';
			} elsif ($self->_has_father($sg_name)) {
				push @extra_headers,'mendelian_inheritance','disease_inheritance','gene','total_variants_in_gene','unique_coord_gene_variants','%_coding_bases_variant','father_allele','parent_allele_common_to_affected','affected_allele_block','definite_compound_het','possible_compound_het','rare_definite_compound_het','rare_possible_compound_het';			
			} else {
				push @extra_headers, 'disease_inheritance','gene','total_variants_in_gene','unique_coord_gene_variants','%_coding_bases_variant';
			}
			$self->{source_info}{$sg_name}{extra_headers} = \@extra_headers;
			my @sample_headers = ();
			#here we include the relation info in the sample header			
			for my $sample (sort keys %{$self->{source_info}{$sg_name}{samples}}) {
				push @sample_headers, $sample.'_alleles ('.$sample_info{$sg_name}{samples}{$sample}{relation}.')' if $self->{pileups};
			}
			$self->{source_info}{$sg_name}{sample_headers} = \@sample_headers;
			
		}
		
	} else {
		for my $sg_name ( keys %sample_info ) {
			if ($self->{source_info}{$sg_name}{affected}) {
				push @extra_headers, 'number_affected_variant/total_affected';
			}
			if ($self->{source_info}{$sg_name}{unaffected}) {
				push @extra_headers, 'number_unaffected_variant/total_unaffected';	
			}
			push @extra_headers,'gene','total_variants_in_gene','unique_coord_gene_variants','%_coding_bases_variant';
			$self->{source_info}{$sg_name}{extra_headers} = \@extra_headers;
			
			my @sample_headers = ();
			for my $sample (sort keys %{$self->{source_info}{$sg_name}{samples}}) {
				push @sample_headers, $sample.'_alleles' if $self->{pileups};
			}
			$self->{source_info}{$sg_name}{sample_headers} = \@sample_headers;
		}	
	}
	
	
	
	
	
}

#Parse the snv and indel reports to load this data
sub parse_reports {
	my ($self) = @_;
	
	my %sample_info = %{$self->{source_info}};
	
	my $found_line = 0;
	
	
	if ($self->{sv}) {
		my $pipe_conf = modules::Pipeline::get_pipe_conf();
		my @sv_types = split(",",$pipe_conf->read('common','sv_types'));
		
		my $svndir = $ENV{'SVNDIR'};
		if (!-d "$svndir/utils") {
			modules::Exception->throw("ERROR: Can't open dir $svndir\n");
		}
		
		my $overlap_bin = "$svndir/utils/overlap_files.pl";
		my $sys_call = modules::SystemCall->new();
		my %sv_header_lookup = ();
		my %sv_index_lookup = ();
		my @sv_report_headers = ();
		my $sv_gene_column_number;
		my $sv_type_column_number;
		my $fail_flag_sv = 0;
		
		
		#Here we need to pick the sample to compare all other in group to and then do individual overlaps
		for my $sg_name ( keys %sample_info ) {
			my $source_type = modules::Pipeline::get_source_type(-source_group_name=>$sg_name);
			
			my $max_distance = $pipe_conf->read($source_type,'cutoffs','sv_max_distance');
			my $proband = $self->_get_proband_name($sg_name);
			my ($sv_xml,$outdir) = fileparse($self->{source_info}{$sg_name}{sv_xml});
			my $tmpdir = $outdir . '/tmp';
			
			if (-d $tmpdir) {
				$sys_call->run("rm -Rf $tmpdir"); 
			}
			$sys_call->run("mkdir $tmpdir");

			my $sv_proband_file = $sample_info{$sg_name}{samples}{$proband}{sv_file};
	
			for my $sv_type (@sv_types) {
				my $command;
				#fussy grepping for tab
				my $regex = "'".'\\t'.${sv_type}.'\\t'."'";
				if ($sv_type eq 'tra') {
					#here we need the second chr
					$command = "grep -P $regex $sv_proband_file | awk -F'\t' '{print ".'$6,$7,$8,$9,$12":"$5'."}' | sort | uniq > $tmpdir/$proband.$sv_type";
				} else {
					$command = "grep -P $regex $sv_proband_file | awk -F'\t' '{print ".'$6,$7,$9,$12":"$5'."}' | sort | uniq > $tmpdir/$proband.$sv_type";
				}
				`$command`;
			}
			
			#Parse the proband sv file first to get the data
			open(PROBAND,"$sv_proband_file") || modules::Exception->throw("Can't open file $sv_proband_file\n");
			
			my ($sv_caller_count,$total_sv_calls,$exon_overlap,$gene_overlap,$sv_unique_key,$priority);
			
			my $line_count = 0;
			
			while (<PROBAND>) {
				chomp;
				next unless /\S/;
				
				#build lookups
				if (/^SV/) {
					@sv_report_headers = split("\t");
					for (my $i = 0; $i < scalar @sv_report_headers; $i++){
						if ($sv_report_headers[$i] eq $self->{sv_gene_col_name}) {
							$sv_gene_column_number = $i;
						}
						if ($sv_report_headers[$i] eq 'sv_type') {
							$sv_type_column_number = $i;
						}
						
						#10 -> ref_allele for example
		    			$sv_header_lookup{$i} = $sv_report_headers[$i];
		    			#ref_allele -> 10 for example
						$sv_index_lookup{$sv_report_headers[$i]} = $;
					}
					next;
				}
				
				if (/---/) {
					if (/^FAIL/ || /^LOW_PRIORITY/) {
						$priority = 'LOW';
						$fail_flag_sv = 1;
						next;
					} elsif (/^HIGH/) {
						$priority = 'HIGH';
						next;
					} elsif (/^MEDIUM/) {
						$priority = 'MEDIUM';
						next;
					}
					next;
				}
	
				if ($fail_flag_sv) {
					next unless $self->{fail};
				}
	
				my @values = split("\t");
				
				
								
				if ($values[0] =~ /\d/) {
					#New entry here
					($sv_caller_count,$total_sv_calls,$exon_overlap,$gene_overlap,$sv_unique_key) = @values;
					$line_count = 0;
				} else {
					#Other line for existing entry
					$sv_unique_key = $values[4];
					$line_count++;
				}
				
				my $gene_name;
				my $sv_type = $values[$sv_type_column_number];
				
				#handle these separately as they as sometimes empty due to grouping
				$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_unique_key}{$priority}{$line_count}{"0:$sv_header_lookup{0}"} = $sv_caller_count;
				$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_unique_key}{$priority}{$line_count}{"1:$sv_header_lookup{1}"} = $total_sv_calls;
				$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_unique_key}{$priority}{$line_count}{"2:$sv_header_lookup{2}"} = $exon_overlap;
				$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_unique_key}{$priority}{$line_count}{"3:$sv_header_lookup{3}"} = $gene_overlap;				
				
				for (my $i = 0; $i < scalar @values; $i++){
					next if $i <= 3; #These are captured above
					my $header_name = $sv_header_lookup{$i};
					if ($header_name =~ /hgnc/) {
						$gene_name = $values[$i];
						
						$self->{source_info}{$sg_name}{gene_lookup}{$sv_unique_key} = $values[$i];
						
						if($gene_name =~ /No gene overlap/){
							$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_total}{$sv_unique_key}='N/A';
							$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{non_uniq_total}='N/A';
	                      	if ($sample_info{$sg_name}{samples}{$proband}{affected} == 1) {
	                        	$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_aff}{$sv_unique_key}='N/A';
	                     	} else {
	                       		$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_unaff}{$sv_unique_key}='N/A';
	                    	}
	
						} else{
                            $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_total}{$sv_unique_key}++;
                            $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{non_uniq_total}++;
							if ($sample_info{$sg_name}{samples}{$proband}{affected} == 1) {
								#$self->{source_info}{gene_counts}{$sg_name}{$values[$i]}{total_aff}++;
								$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_aff}{$sv_unique_key}++;
							} else {
								#$self->{source_info}{gene_counts}{$sg_name}{$values[$i]}{total_unaff}++;
								$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_unaff}{$sv_unique_key}++;
							}
						}
					} 
					
					
					#Set the value
					my %sv_headers = %{$self->{source_info}{$sg_name}{sv_headers}};
					my $index = $sv_headers{$header_name}; #Get the index from the first file seen as that represents master order
					
					if (!defined $index) {
						print "$proband $header_name\n";
					}
					$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_unique_key}{$priority}{$line_count}{"$index:$sv_header_lookup{$i}"} = $values[$i];
				}
				
				
			}
			
			
			#Keep track of which sample entries match the proband entries for later parsing of files
			my %proband_pairs = ();
			
			#Now generate the overlaps between each proband/sample pair -> we only focus on SVs found in proband as overlapping is too complex otherwise
			for my $sample ( sort keys %{$sample_info{$sg_name}{samples}} ) {
				if ($sample eq $proband) {
					next;
				} #Skip the proband as it is the baseline and already done
				
				my $sv_sample_file = $sample_info{$sg_name}{samples}{$sample}{sv_file};
				
				for my $sv_type (@sv_types) {
					my $command;
					#fussy grepping for tab
					my $regex = "'".'\\t'.${sv_type}.'\\t'."'";
					if ($sv_type eq 'tra') {
						#here we need the second chr
						$command = "grep -P $regex $sv_sample_file | awk -F'\t' '{print ".'$6,$7,$8,$9,$12":"$5'."}' | sort | uniq > $tmpdir/$sample.$sv_type";
					} else {
						$command = "grep -P $regex $sv_sample_file | awk -F'\t' '{print ".'$6,$7,$9,$12":"$5'."}' | sort | uniq > $tmpdir/$sample.$sv_type";
					}

					`$command`;
					
					if ($sv_type eq 'tra') {
						open(FILE1,"$tmpdir/$proband.$sv_type") || modules::Exception->throw("Can't open file $proband.$sv_type\n");
						open(FILE2,"$tmpdir/$sample.$sv_type") || modules::Exception->throw("Can't open file $sample.$sv_type\n");
						
						my %trans = ();
						
						while (<FILE1>) {
							chomp;
							my ($chr1,$coord1,$chr2,$coord2,$tmp) = split(" ");
							my ($id,$sv_unique_key) = split(':',$tmp);
							my $chr_pair;
							my $key;
							if ($chr1 lt $chr2) {
								$chr_pair = $chr1.':'.$chr2;
								$key = "$coord1:$coord2:$id:$sv_unique_key";
							} else {
								$chr_pair = $chr2.':'.$chr1;
								$key = "$coord2:$coord1:$id:$sv_unique_key";
							}
							$trans{$chr_pair}{proband}{$key}++;
							$self->{source_info}{$sg_name}{var_counts}{sv}{$sv_unique_key}{aff} = 1;
						}
						
						while (<FILE2>) {
							chomp;
							my ($chr1,$coord1,$chr2,$coord2,$tmp) = split(" ");
							my ($id,$sv_unique_key) = split(':',$tmp);
							my $chr_pair;
							my $key;
							if ($chr1 lt $chr2) {
								$chr_pair = $chr1.':'.$chr2;
								$key = "$coord1:$coord2:$id:$sv_unique_key";
							} else {
								$chr_pair = $chr2.':'.$chr1;
								$key = "$coord2:$coord1:$id:$sv_unique_key";
							}
							$trans{$chr_pair}{sample}{$key}++;
						}
						
						#Now check if the proband/sample event should be merged
						for my $chr_pair(keys %trans) {
							
							
							#If the chr_pair is found in both samples
							if (keys %{$trans{$chr_pair}} == 2) {
								
								for my $proband_str (keys %{$trans{$chr_pair}{proband}}) {
									my ($pro_coord1,$pro_coord2,$pro_id,$pro_event) = split(':',$proband_str);
									for my $sample_str (keys %{$trans{$chr_pair}{sample}}) {
										my ($sample_coord1,$sample_coord2,$sample_id,$sample_event) = split(':',$sample_str);
										
										if (abs($pro_coord1-$sample_coord1) < 500 && abs($pro_coord2-$sample_coord2) < 500) {
											#Here the events should be merged -> use proband coords for reference
											
											#Add the sample to match lookup table for later parsing
											$proband_pairs{$sample}{$sample_event}{$pro_event}++;
											#print "here $chr_pair P $proband_str S $sample_str EVENT $pro_event\n";
											#Record whether sample is affected or not
											if ($self->{source_info}{$sg_name}{samples}{$sample}{affected}) {
												$self->{source_info}{$sg_name}{var_counts}{sv}{$pro_event}{aff}++;
											} 
										} 
										
									}
								}
								
								
							} 
						}
						#print Dumper \%proband_pairs;
					}  else {
						#Dels/dups/invs get normal handling
						my $out = "$tmpdir/${proband}_${sample}.$sv_type";
						$sys_call->run("$overlap_bin -just_overlap -max_distance $max_distance -ref $tmpdir/$proband.$sv_type -coord $tmpdir/$sample.$sv_type -all -just_overlap > $out");
						open(OUT,"$out") || modules::Exception->throw("Can't open file $out\n");
						while (<OUT>) {
							chomp;
							my @fields = split("\t",$_);
							my (undef,$sv_unique_key) = split(':',$fields[3]);
							
							#Set the affected counter for the proband
							$self->{source_info}{$sg_name}{var_counts}{sv}{$sv_unique_key}{aff} = 1;
							
							if ($fields[4] ne '1') {
								
								my @matches = split('\^\^\^',$fields[4]);
								
								
								#Match for proband -> record sample info
								for my $match (@matches) {
									my (undef,$sample_event) = split(':',$match);
									$proband_pairs{$sample}{$sample_event}{$sv_unique_key}++;
								}

								#Record whether sample is affected or not
								if ($self->{source_info}{$sg_name}{samples}{$sample}{affected}) {
									$self->{source_info}{$sg_name}{var_counts}{sv}{$sv_unique_key}{aff}++;
								} 
							}	
						}
						
						close OUT;
						
						
					}
					
				}
			}
			
			#print Dumper \%proband_pairs;
			
			
			#Now use overlap info to get the lines from the non-proband samples
			for my $sample ( keys %proband_pairs ) {
				$found_line = 1;
				
				my $sv_sample_file = $sample_info{$sg_name}{samples}{$sample}{sv_file};
				
				#Now parse the sample sv file to read in lines
				open(SAMPLE,"$sv_sample_file") || modules::Exception->throw("Can't open file $sv_sample_file\n");
				while (<SAMPLE>) {
					chomp;
					next unless /\S/;
					my $sample_event;
					my @values = split("\t");
					if ($values[0] =~ /\d/) {
						($sv_caller_count,$total_sv_calls,$exon_overlap,$gene_overlap,$sample_event) = @values;
					} else {
						$sample_event = $values[4];
					}
				
					
				
					if (exists $proband_pairs{$sample}{$sample_event}) {
						for my $proband_group (keys %{$proband_pairs{$sample}{$sample_event}}) {
							#Only record cases that overlap with proband cases
							my $sample_line;
							if ($values[0] =~ /\d/) {
								$sample_line = join("\t",
													'',
													'',
													'',
													'',
													@values
													);
							} else {
								my @sample_line = ();
								push @sample_line, $sv_caller_count, $total_sv_calls, $exon_overlap, $gene_overlap;
								for ( my $var = 4 ; $var < @values ; $var++ ) {
								    push @sample_line, $values[$var];
								}
																	

								$sample_line = join("\t",
													'',
													'',
													'',
													'',
													@sample_line
													)	
							}
							#Assumes same sv_type column index
							push @{$self->{source_info}{$sg_name}{var_data}{sv}{$values[$sv_type_column_number]}{$proband_group}{samples}{$sample}}, $sample_line;
						}
						
					} else {
						#For now don't bother bother with non-proband entries
					}
								
				}
			    	
			}
		}
		#print Dumper $self;		
		
	} else {
		for my $sg_name ( keys %sample_info ) {
	
			for my $sample ( sort keys %{$sample_info{$sg_name}{samples}} ) {
				## skip not included samples if report type is 'sample'
				if ($self->{report_type} eq 'sample' && !$self->{source_info}{$sg_name}{samples}{$sample}{include}){
					next;
				}                        

				print STDERR "Sample $sample reports\n";

				#SNV/INDEL reports
				open(SNV_REPORT,$sample_info{$sg_name}{samples}{$sample}{snv_file}) || modules::Exception->throw("ERROR: Can't open snv file $sg_name $sample $sample_info{$sg_name}{samples}{$sample}{snv_file}");	
				print "reading snv file $sample_info{$sg_name}{samples}{$sample}{snv_file}\n";
				my $fail_flag_snv = 0;
				my %snv_header_lookup = ();
				my %snv_index_lookup = ();
				my @snv_report_headers = ();
				
				
				my $snv_gene_column_number;
				while (<SNV_REPORT>) {
					chomp;
					if (/^chr/) {
						@snv_report_headers = split("\t");
						for (my $i = 0; $i < scalar @snv_report_headers; $i++){
							if ($snv_report_headers[$i] eq $self->{snv_gene_col_name}) {
								$snv_gene_column_number = $i;
							}
							#10 -> ref_allele for example
		    				$snv_header_lookup{$i} = $snv_report_headers[$i];
		    				#ref_allele -> 10 for example
							$snv_index_lookup{$snv_report_headers[$i]} = $;
						}
						next;
					}
					
					if (/---/) {
						if (/FAIL/ || /LOW_PRIORITY/) {
							$fail_flag_snv = 1;
							next;
						}
						if (/PASS/ || /NOVEL/) {
							next;
						}
						next;
					}
					
		
					if ($fail_flag_snv) {
						next unless $self->{fail};
					}
		
					my @values = split("\t");
					my ($chr,$coord,$ref,undef,$var) = @values;
					my $start = my $end = $coord;
	#				if (/COMBINED:(\d+)-(\d+)/) {
	#					$start = $1;
	#					$end = $2;
	#				} 
					
					if ($self->{filter_output}) {
						if ($self->{filters}{gene_list}) {
							#Make sure the gene is in the list
							my $gene_name = $values[$snv_gene_column_number];
							next unless exists $self->{filters}{gene_list}{$values[$snv_gene_column_number]};
						}
						
						
						#Running on single chromosome
						if ($self->{filters}{chr}) {
							next unless $self->{filters}{chr} eq $chr;
						}
						
						if ($self->{filters}{start}) {
							next unless ($self->{filters}{start} <= $coord && $self->{filters}{end} >= $coord); 
						}
					}
					
					
					$found_line = 1;
					
					my $snv_unique_key = "$chr:$start:$end:$var";
					#$self->{source_info}{$sg_name}{snv}{$snv_unique_key}{total}++;
					
					if ($sample_info{$sg_name}{samples}{$sample}{affected} == 1) {
						$self->{source_info}{$sg_name}{var_counts}{snv}{$snv_unique_key}{aff}++;
					} else {
						$self->{source_info}{$sg_name}{var_counts}{snv}{$snv_unique_key}{unaff}++;
					}
					
					my $skip_line = 0;
					my $gene_name;
					
					
					for (my $i = 1; $i < scalar @values; $i++){
						next if $snv_header_lookup{$i} eq 'chr' || $snv_header_lookup{$i} eq 'coord';
						my $header_name = $snv_header_lookup{$i};
						if ($header_name =~ /hgnc/) {
							$gene_name = $values[$i];
							
							$self->{source_info}{$sg_name}{gene_lookup}{$snv_unique_key} = $values[$i];
							
							if($gene_name =~ /No gene overlap/){
								$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_total}{$snv_unique_key}='N/A';
								$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{non_uniq_total}='N/A';
		                      	if ($sample_info{$sg_name}{samples}{$sample}{affected} == 1) {
		                        	$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_aff}{$snv_unique_key}='N/A';
		                     	} else {
		                       		$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_unaff}{$snv_unique_key}='N/A';
		                    	}
		
							} else{
	                            $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_total}{$snv_unique_key}++;
	                            $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{non_uniq_total}++;
								if ($sample_info{$sg_name}{samples}{$sample}{affected} == 1) {
									#$self->{source_info}{gene_counts}{$sg_name}{$values[$i]}{total_aff}++;
									$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_aff}{$snv_unique_key}++;
								} else {
									#$self->{source_info}{gene_counts}{$sg_name}{$values[$i]}{total_unaff}++;
									$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_unaff}{$snv_unique_key}++;
								}
							}
						} elsif (exists $self->{filter_output}) {
							#Check for optional filtering steps
							if (exists $self->{filters}{polyphen} && $header_name eq 'polyphen_prediction') {
								$skip_line = 1 unless  $self->{filters}{polyphen} eq $values[$i];	
							}
							if (exists $self->{filters}{sift} && $header_name eq 'sift_prediction') {
								$skip_line = 1 unless  $self->{filters}{sift} eq $values[$i];	
							}
							if (exists $self->{filters}{max_allele_freq} && $header_name eq 'dbsnp_var_allele_freq') {
								if ($values[$i] =~ /\d/) {
									$skip_line = 1 unless $self->{filters}{max_allele_freq} >= $values[$i];
								} else {
									$skip_line = 1;
								}	
							}
						} 
						
						if ($header_name =~ /aa_length/) {
							$self->{source_info}{$sg_name}{gene_length}{$gene_name} = $values[$i]; 
						}
						
						
						#Set the value
						my %snv_headers = %{$self->{source_info}{$sg_name}{snv_headers}};
						my $index = $snv_headers{$header_name}; #Get the index from the first file seen as that represents master order
						
						if (!defined $index) {
							print "$sample $header_name\n";
						}
						$self->{source_info}{$sg_name}{var_data}{snv}{$snv_unique_key}{data}{$sample}{"$index:$snv_header_lookup{$i}"} = $values[$i];
						
					}
					
					if (!$skip_line) {
						$self->{source_info}{$sg_name}{pileup_coord}{"$chr:$coord"}++;
					} else {
						#Delete the entry
						delete $self->{source_info}{$sg_name}{var_data}{snv}{$snv_unique_key};
					}	
					
				}
				close SNV_REPORT;
				
				
				
				my %indel_header_lookup = ();		
				my @indel_report_headers = ();
				my $fail_flag_indel;
				my $indel_gene_column_number;
				
				if (exists $self->{filters}{polyphen} || exists $self->{filters}{sift}) {
					next; #Don't generate indel reports if filtering on polyphen or sift
				}
							
				
				open(INDEL_REPORT,$sample_info{$sg_name}{samples}{$sample}{indel_file}) || modules::Exception->throw("ERROR: Can't open indel file $sample_info{$sg_name}{samples}{$sample}{indel_file}");
				print "reading indel file $sample_info{$sg_name}{samples}{$sample}{indel_file}\n";
				while (<INDEL_REPORT>) {
					chomp;
					if (/^chr/) {
						@indel_report_headers = split("\t");
						for (my $i = 0; $i < scalar @indel_report_headers; $i++){
							if ($indel_report_headers[$i] eq $self->{indel_gene_col_name}) {
								$indel_gene_column_number = $i;
							}
		    				$indel_header_lookup{$i} = $indel_report_headers[$i];
						}
						next;
					}
					
					if (/---/) {
						if (/FAIL/ || /LOW_PRIORITY/) {
							$fail_flag_indel = 1;
							next;
						}
						if (/PASS/ || /NOVEL/) {
							next;
						}
						next;
					}
		
		
					if ($fail_flag_indel) {
						next unless $self->{fail};
					}
		
					my @values = split("\t");
					my ($chr,$start_coord,$end_coord,undef,undef,$type,undef,$var_base) = @values;
					
					if ($self->{filter_output}) {
						if ($self->{filters}{gene_list}) {
							#Make sure the gene is in the list
							next unless exists $self->{filters}{gene_list}{$values[$indel_gene_column_number]};
						}
						
						#Running on single chromosome
						if ($self->{filters}{chr}) {
							next unless $self->{filters}{chr} eq $chr;
						}
						
						if ($self->{filters}{start}) {
							next unless ($self->{filters}{start} <= $start_coord && $self->{filters}{end} >= $end_coord); 
						}
						
					}
					$found_line = 1;
					
					my $indel_unique_key = "$chr:$start_coord:$end_coord:$var_base";
					#$self->{source_info}{$sg_name}{var_counts}{$type}{$indel_unique_key}{total}++;
					
					if ($sample_info{$sg_name}{samples}{$sample}{affected} == 1) {
						$self->{source_info}{$sg_name}{var_counts}{$type}{$indel_unique_key}{aff}++;
					} else {
						$self->{source_info}{$sg_name}{var_counts}{$type}{$indel_unique_key}{unaff}++;
					}
					
					
					my $skip_line = 0;
					my $gene_name;
					
					for (my $i = 0; $i < scalar @values; $i++){
						next if $indel_header_lookup{$i} eq 'chr' || $indel_header_lookup{$i} eq 'start_coord' || $indel_header_lookup{$i} eq 'end_coord';
						my $header_name = $indel_header_lookup{$i};
						if ($header_name =~ /hgnc/) {
							$gene_name = $values[$i];
							$self->{source_info}{$sg_name}{gene_lookup}{$indel_unique_key} = $values[$i];
	
							if($gene_name =~ /No gene overlap/){
	                                                        $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_total}{$indel_unique_key}='N/A';
	                                                        $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{non_uniq_total}='N/A';
	                                                        if ($sample_info{$sg_name}{samples}{$sample}{affected} == 1) {
	                                                                $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_aff}{$indel_unique_key}='N/A';
	                                                        } else {
	                                                                $self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_unaff}{$indel_unique_key}='N/A';
	                                                        }
							
							}else{
								$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_total}{$indel_unique_key}++;
								$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{non_uniq_total}++;
								if ($sample_info{$sg_name}{samples}{$sample}{affected} == 1) {
									#$self->{source_info}{gene_counts}{$sg_name}{$values[$i]}{total_aff}++;
									$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_aff}{$indel_unique_key}++;
								} else {
									#$self->{source_info}{gene_counts}{$sg_name}{$values[$i]}{total_unaff}++;
									$self->{source_info}{$sg_name}{gene_counts}{$values[$i]}{uniq_unaff}{$indel_unique_key}++;
								}
							}
						} elsif (exists $self->{filter_output}) {
							#Check for optional filtering steps
							
							if (exists $self->{filters}{max_allele_freq} && $header_name eq 'dbsnp_var_allele_freq') {
								if ($values[$i] =~ /\d/) {
									$skip_line = 1 unless  $self->{filters}{max_allele_freq} >= $values[$i];
								} else {
									$skip_line = 1;;
								}	
							}
						} 
						
						if ($header_name =~ /aa_length/) {
							$self->{source_info}{$sg_name}{gene_length}{$gene_name} = $values[$i];
						}
						
						my %indel_headers = %{$self->{source_info}{$sg_name}{indel_headers}};
						my $index = $indel_headers{$header_name}; #Get the index from the first file seen as that represents master order
						$self->{source_info}{$sg_name}{var_data}{$type}{$indel_unique_key}{data}{$sample}{"$index:$indel_header_lookup{$i}"} = $values[$i];
					}
					if (!$skip_line) {
						my $pileup_start = $start_coord-1; #Account for how indels are reported
						if ($type eq 'DEL') {
							$self->{source_info}{$sg_name}{pileup_coord}{"$chr:$pileup_start"}++;
						} else {
							$self->{source_info}{$sg_name}{pileup_coord}{"$chr:$start_coord"}++;
						}
						
					} else {
						delete $self->{source_info}{$sg_name}{var_data}{$type}{$indel_unique_key};				
					}
				}
				close INDEL_REPORT;
				
			}
			
		}
	}
	return $found_line;
}


#subroutine creates the pileup files for getting allele strings later
sub generate_pileups {
	my ($self) = @_;
	
	my %sample_info = %{$self->{source_info}};
	#Pileup for snvs
	my $pipe_config = modules::Pipeline::get_pipe_conf();
	my $clus_conf = modules::Pipeline::get_cluster_conf();
	my $samtools_bin = $pipe_config->read('human_related_gatk','binaries','samtools','binary');

	my $ref = $clus_conf->read('human_related_gatk','svn','fasta_file');

	my $sys_call = modules::SystemCall->new();
	for my $sg_name ( keys %sample_info ) {
	
		my $pileup_coord_file = $self->{source_info}{$sg_name}{dir_base} . '/'. $sg_name .'.coord';
		
		open (COORD,">$pileup_coord_file") || modules::Exception->throw("Can't open file $pileup_coord_file");
		for my $coord_string (sort keys %{$self->{source_info}{$sg_name}{pileup_coord}}) {
			my ($chr,$start_coord) = split(":",$coord_string);
			
			my $full_chr = $chr;	#the hs37d5 genome has no 'chr' in front, just chr number!
                              
			print COORD "$full_chr\t$start_coord\n";
		}
		close COORD;
		
		for my $sample ( sort keys %{$sample_info{$sg_name}{samples}} ) {
			my $pileup_file = $self->{source_info}{$sg_name}{dir_base} . '/' . $sample .'.pileup';
			my $bam = $self->{source_info}{$sg_name}{samples}{$sample}{bam_file};
			my $mpileup_snv_command = "$samtools_bin mpileup -A -E -f $ref -l $pileup_coord_file $bam  > $pileup_file";
			$sys_call->run($mpileup_snv_command) unless -e $pileup_file; #Don't run longish command if not needed
			open(PILEUP,"$pileup_file") || modules::Exception->throw("Can't open pileup file $pileup_file");
			while (<PILEUP>) {
				$_ =~ s/^chr//;
				my @fields = split("\t");
				my ($pileup_string,$zyg) = modules::Utils->pileup_string($fields[4]);
				#print "$fields[0] $fields[1] $sample $fields[4] $pileup_string $zyg\n";
				$self->{source_info}{$sg_name}{pileup_str}{$fields[0].':'.$fields[1]}{$sample}{pileup} = $pileup_string;
				$self->{source_info}{$sg_name}{pileup_str}{$fields[0].':'.$fields[1]}{$sample}{zyg} = $zyg;
			}
			close PILEUP;
		
		}	
		
		my %pileup_strings = %{$self->{source_info}{$sg_name}{pileup_str}};
		
		#Set the 'no data' entries
		for my $pileup_str (keys %pileup_strings) {
			for my $sample ( keys %{$self->{source_info}{$sg_name}{samples}} ) {
				if (!exists $pileup_strings{$pileup_str}{$sample}) {
					$self->{source_info}{$sg_name}{pileup_str}{$pileup_str}{$sample}{pileup} = 'No data';
					$self->{source_info}{$sg_name}{pileup_str}{$pileup_str}{$sample}{zyg} = 'No data';
				}
			}
		}
	}
}	


#subroutine gets the raw parent alleles and sets the mendelian flag	
sub get_parent_alleles {
	
	my ($self) = @_;
	#print Dumper $self;	
	
	for my $sg_name (keys %{$self->{source_info}}) {
	
		my %relations = ();
		
		
		
		#Find out mother and father included in cohort
		for my $sample (keys %{$self->{source_info}{$sg_name}{samples}}) {
			$relations{$sample.':'.$self->{source_info}{$sg_name}{samples}{$sample}{relation}}++;
		}
		
		
		
		for my $var_type (keys %{$self->{source_info}{$sg_name}{var_data}}) {
			for my $var_key (keys %{$self->{source_info}{$sg_name}{var_data}{$var_type}}) {
			#print "VARKEY $var_key\n";	
				my $gmaf_snv_str = my $gmaf_indel_str;
				my $gmaf_header_name = 'gmaf_1000_genomes';
				
				my $filter_dbsnp_snv_str = my $filter_dbsnp_indel_str;
				my $filter_dbsnp_snv_header_name = 'filter_dbsnp_snv'; #NOVEL,RARE_ALLELE,NO_FREQ,RARE_REF,FAIL
				my $filter_dbsnp_indel_header_name = 'filter_dbsnp_indel';
				
				my $gmaf_snv_index = $self->{source_info}{$sg_name}{snv_headers}{$gmaf_header_name};
				my $gmaf_indel_index = $self->{source_info}{$sg_name}{indel_headers}{$gmaf_header_name};
				
				my $filter_dbsnp_snv_index = $self->{source_info}{$sg_name}{snv_headers}{$filter_dbsnp_snv_header_name};
				my $filter_dbsnp_indel_index = $self->{source_info}{$sg_name}{indel_headers}{$filter_dbsnp_indel_header_name};
				
				#Get the gmaf data first
				for my $sample (keys %{$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}}) {
					if (exists $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$gmaf_snv_index:$gmaf_header_name"}) {
						$gmaf_snv_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$gmaf_snv_index:$gmaf_header_name"};
					}
					if (exists $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$gmaf_indel_index:$gmaf_header_name"}) {
						$gmaf_indel_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$gmaf_indel_index:$gmaf_header_name"};
					}	
					if (exists $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$filter_dbsnp_snv_index:$filter_dbsnp_snv_header_name"}) {
						$filter_dbsnp_snv_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$filter_dbsnp_snv_index:$filter_dbsnp_snv_header_name"};
					}
					if (exists $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$filter_dbsnp_indel_index:$filter_dbsnp_indel_header_name"}) {
						$filter_dbsnp_indel_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{"$filter_dbsnp_indel_index:$filter_dbsnp_indel_header_name"};
					}					
				}
			
				
				my ($chr,$coord) = split(':',$var_key);
				
				my %family_allele_data = ();
				
				
				my $aff_count = 0;
				my $unaff_count = 0;
	
				my $recessive = 1; #turn on by default until we see something disproving it
				my $mother_var = 0; #flag if mother variant
				my $father_var = 0; #flag if father variant
				
				#Lookup for pileup
				my $pileup_key;
				
				if ($var_type eq 'DEL') {
					my $indel_coord = $coord - 1;
					$pileup_key = "$chr:$indel_coord";
				} else {
					$pileup_key = "$chr:$coord";
				}
	
				#Try to figure out phasing; set everything to REF by default if no data
				for my $rel_key (keys %relations) {
					my ($sample) = split(":",$rel_key);
					if (exists $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}) {
						$family_allele_data{$rel_key}{zyg} = 	$self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg};
					} else {
						#Handling for the case where no pileup data available (GATK indel realign)
						$family_allele_data{$rel_key}{zyg} = 'ref';
					} 
				}
				
				#Get parents info if available 
				my %parent_data = ();
				
				for my $family_key (keys %family_allele_data) {
					my ($sample,$relation) = split(':',$family_key);
					if ($relation eq 'mother') {
						
						my ($mother_first_allele,$mother_second_allele) = $self->_get_alleles($family_allele_data{$family_key}{zyg},$self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup});
						$parent_data{'mother'}{'zyg'} = $family_allele_data{$family_key}{zyg};
						$parent_data{'mother'}{'first_allele'} = $mother_first_allele;
						$parent_data{'mother'}{'second_allele'} = $mother_second_allele;
						$parent_data{'mother'}{'pileup'} = $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup};
					} elsif ($relation eq 'father') {
						my ($father_first_allele,$father_second_allele) = $self->_get_alleles($family_allele_data{$family_key}{zyg},$self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup});
						$parent_data{'father'}{'zyg'} = $family_allele_data{$family_key}{zyg};
						$parent_data{'father'}{'first_allele'} = $father_first_allele;
						$parent_data{'father'}{'second_allele'} = $father_second_allele;
						$parent_data{'father'}{'pileup'} = $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup};
					}
				}
	
				my $disease_inheritance;
				#Check for inheritance patterns
				$disease_inheritance = $self->_disease_inheritance(\%family_allele_data, $pileup_key, $sg_name);
	
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{disease_inheritance} = $disease_inheritance;
	
	
				#Get mendelian inheritance pattern; need to account for which parent(s) are sequenced
				my $mendel_inheritance = 'mendelian_rules_followed';
				if ($self->_has_single_parent($sg_name)) {
					#Check if it's mendelian
					my $mendel_str = $self->_mendel_inheritance(\%family_allele_data,$sg_name);
					if ($mendel_str ne 'yes') {
						$mendel_inheritance = $mendel_str;
					} 
					
				} else {
					$mendel_inheritance = 'No Parental info';
				}
				
				#This is the inheritance based on what pileup shows
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{mendel_inheritance} = $mendel_inheritance;
				
				
				#Parent alleles to affected
				my $affected_alleles_from_one_parent = 'yes'; #default is yes; set to zero if misproven or ? if alleles unknown
				my %affected_allele_data = ();
				
				
				my $del_com_het = 0;
				#Get the mother and father alleles
				for my $family_key (keys %family_allele_data) {
					my ($sample,$relation) = split(':',$family_key);
					next if $relation eq 'mother' || $relation eq 'father';
					my $sz = $family_allele_data{$family_key}{zyg};
					
					my $pileup_data = $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup};
					
					my ($child_first_allele,$child_second_allele) = $self->_get_alleles($sz,$pileup_data);
					
					my $mother_allele = '?';
					my $father_allele = '?';	
				
					
					my $het_key;
					my $rare_het_key; #collection of RARE_ALLELE,NOVEL,NO_FREQ variants
					
					if ($var_type eq 'snv') {
						$het_key = "$chr:$coord:$var_type:gmaf=$gmaf_snv_str";
						$rare_het_key = "$chr:$coord:$var_type:gmaf=$gmaf_snv_str:filter_dbsnp=$filter_dbsnp_snv_str";
						
					} else {
						$het_key = "$chr:$coord:$var_type:gmaf=$gmaf_indel_str";
						$rare_het_key = "$chr:$coord:$var_type:gmaf=$gmaf_indel_str:filter_dbsnp=$filter_dbsnp_indel_str";
					}
				
					#$self->{source_info}{$sg_name}{var_data}{snv}{$snv_unique_key}{data}{$sample}{"$index:$snv_header_lookup{$i}"}
					
					
					if ($self->_has_two_parents($sg_name)) {
						my $mz = $parent_data{'mother'}{'zyg'};
						my $fz = $parent_data{'father'}{'zyg'};
						
					
						
						#Data for both parents; easy case
						if (($mz eq 'ref' && $fz eq 'ref' && $sz eq 'ref') || ($mz eq 'het' && $fz eq 'ref' && $sz eq 'ref') || ($mz eq 'ref' && $fz eq 'het' && $sz eq 'ref') || ($mz eq 'het' && $fz eq 'het' && $sz eq 'ref')) {
							#ref/ref cases
							$mother_allele = $parent_data{'mother'}{'first_allele'};
							$father_allele = $parent_data{'father'}{'first_allele'};
						} elsif (($mz eq 'hom' && $fz eq 'het' && $sz eq 'hom') || ($mz eq 'het' && $fz eq 'hom' && $sz eq 'hom') || ($mz eq 'hom' && $fz eq 'hom' && $sz eq 'hom') || ($mz eq 'het' && $fz eq 'het' && $sz eq 'hom')) {
							#hom/hom cases
							$mother_allele = $parent_data{'mother'}{'second_allele'};
							$father_allele = $parent_data{'father'}{'second_allele'};
						} elsif (($mz eq 'ref' && $fz eq 'het' && $sz eq 'het') || ($mz eq 'ref' && $fz eq 'hom' && $sz eq 'het') || ($mz eq 'het' && $fz eq 'hom' && $sz eq 'het')) {
							#father variant cases
							$mother_allele = $parent_data{'mother'}{'first_allele'};
							$father_allele = $parent_data{'father'}{'second_allele'};
							if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
								$affected_allele_data{father}++;
							} else {
								$affected_alleles_from_one_parent = 'no';
							}
						} elsif (($mz eq 'het' && $fz eq 'ref' && $sz eq 'het') || ($mz eq 'hom' && $fz eq 'ref' && $sz eq 'het') || ($mz eq 'hom' && $fz eq 'het' && $sz eq 'het')) {
							#mother variant cases
							$mother_allele = $parent_data{'mother'}{'second_allele'};
							$father_allele = $parent_data{'father'}{'first_allele'};
							if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
								$affected_allele_data{mother}++;
							} else {
								$affected_alleles_from_one_parent = 'no';
							}
						} 
						
						
						if ($sz eq 'het' && $self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/) {
							#Flag if it's definite compound het 
							if (($mz eq 'het' && $fz eq 'ref') || ($mz eq 'hom' && $fz eq 'ref') || ($mz eq 'hom' && $fz eq 'het')) {
								#mother variant allele
								$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{mother}{$het_key}++;
								if (defined $filter_dbsnp_snv_str){
									if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{mother}{$rare_het_key}++;
									}
								}
								if (defined $filter_dbsnp_indel_str){
									if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{mother}{$rare_het_key}++;
									}
								}
							} 
							
							if (($fz eq 'het' && $mz eq 'ref') || ($fz eq 'hom' && $mz eq 'ref') || ($fz eq 'hom' && $mz eq 'het')) {
								#father variant allele
								$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{father}{$het_key}++;
								if (defined $filter_dbsnp_snv_str){
									if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{father}{$rare_het_key}++;
									}
								}
								
								if (defined $filter_dbsnp_indel_str){
									if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{father}{$rare_het_key}++;
									}
								}

							}
			
							if ($fz eq 'het' && $mz eq 'het') {
								#ambiguous case
								$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{either}{$het_key}++;
								if (defined $filter_dbsnp_snv_str ){
									if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{either}{$rare_het_key}++;
									}
								}
								
								if (defined $filter_dbsnp_indel_str){
									if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{either}{$rare_het_key}++;
									}
								}

							}
							
							
						} 
						
						

						if ($father_allele eq '?' || $mother_allele eq '?') {
							$affected_alleles_from_one_parent = 'no';
						}
						
						#If unaffected and not reference or unknown
						if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 0) {
							if ($father_allele ne 'ref' || $mother_allele ne 'ref') {
								$affected_alleles_from_one_parent = 'no';
							}
						}
						
						#If affected and is reference or unknown
						if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
							if ($father_allele eq 'ref' && $mother_allele eq 'ref') {
								$affected_alleles_from_one_parent = 'no';
							}
						}
						
					} elsif ($self->_has_mother($sg_name)) {
						my $mz = $parent_data{'mother'}{'zyg'};
						
						#Data for mother available; can still deduce some cases
						if (($mz eq 'ref' && $sz eq 'ref') || ($mz eq 'het' && $sz eq 'ref') || ($mz eq 'ref' && $sz eq 'het')) {
							$mother_allele = $parent_data{'mother'}{'first_allele'}; #inherits ref allele
						} elsif (($mz eq 'het' && $sz eq 'hom') || ($mz eq 'hom' && $sz eq 'het') || ($mz eq 'hom' && $sz eq 'hom')) {
							$mother_allele = $parent_data{'mother'}{'second_allele'}; #inherits variant allele
						}
						
						
						if ($sz eq 'het') {
							#Flag if it's definite compound het 
							if ($mz eq 'hom') {
								if($self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/){
									#mother variant allele
									$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{mother}{$het_key}++;
									if (defined $filter_dbsnp_snv_str){
										if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{mother}{$rare_het_key}++;
										}
									}
									
									if (defined $filter_dbsnp_indel_str){
										if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{mother}{$rare_het_key}++;
										}
									}

								}
								if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
									$affected_allele_data{mother}++;
								} else {
									$affected_alleles_from_one_parent = 'no';
								}
							} 
							
							#Here the mutant allele must have come from the father
							if ($mz eq 'ref') {
								if($self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/){

									$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{father}{$het_key}++;
									if (defined $filter_dbsnp_snv_str) {
										if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{father}{$rare_het_key}++;
										}
									}
									
									if (defined $filter_dbsnp_indel_str){
										if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{father}{$rare_het_key}++;
										}
									}

								}
								if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
									$affected_allele_data{father}++;
								} else {
									$affected_alleles_from_one_parent = 'no';
								}
							}
							
							if ($mz eq 'het') {
								if($self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/){
									#ambiguous case
									$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{either}{$het_key}++;
									if (defined $filter_dbsnp_snv_str) {
										if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{either}{$rare_het_key}++;
										}
									}
									
									if (defined $filter_dbsnp_indel_str){
										if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{either}{$rare_het_key}++;
										}
									}

								}
							}
							
							
						}
						
						if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 0 && $sz ne 'ref' && $sz ne '?') {
							$affected_alleles_from_one_parent = 'no';
						}
						
						
					} elsif ($self->_has_father($sg_name)) {
						my $fz = $parent_data{'father'}{'zyg'};
						
						#Data for father available; can still deduce some cases
						if (($fz eq 'ref' && $sz eq 'ref') || ($fz eq 'het' && $sz eq 'ref') || ($fz eq 'ref' && $sz eq 'het')) {
							$father_allele = $parent_data{'father'}{'first_allele'}; #inherits ref allele
						} elsif (($fz eq 'het' && $sz eq 'hom') || ($fz eq 'hom' && $sz eq 'het') || ($fz eq 'hom' && $sz eq 'hom')) {
							$father_allele = $parent_data{'father'}{'second_allele'}; #inherits variant allele
						}
						
						
						if ($sz eq 'het') {
							#Flag if it's definite compound het 
							if ($fz eq 'hom') {
								if($self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/){
									#father variant allele
									$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{father}{$het_key}++;
									if (defined $filter_dbsnp_snv_str){
										if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
										$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{father}{$rare_het_key}++;
										}
									}
									
									if (defined $filter_dbsnp_indel_str){
										if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{father}{$rare_het_key}++;
										}
									}

								}
								if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
									$affected_allele_data{father}++;
								} else {
									$affected_alleles_from_one_parent = 'no';
								}
							} 
							
							#Here the mutant allele must have come from the mother
							if ($fz eq 'ref') {
								if($self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/){
									$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{mother}{$het_key}++;
									if (defined $filter_dbsnp_snv_str) {
										if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{mother}{$rare_het_key}++;
										}
									}
									
									if (defined $filter_dbsnp_indel_str){
										if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{mother}{$rare_het_key}++;
										}
									}
	
								}
								if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
									$affected_allele_data{mother}++;
								} else {
									$affected_alleles_from_one_parent = 'no';
								}
							}
							
							if ($fz eq 'het') {
								if($self->{source_info}{$sg_name}{gene_lookup}{$var_key} !~ /No gene overlap/){
									#ambiguous case
									$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{either}{$het_key}++;
									if (defined $filter_dbsnp_snv_str){		
										if ($filter_dbsnp_snv_str eq 'NOVEL' || $filter_dbsnp_snv_str eq 'RARE_ALLELE' || $filter_dbsnp_snv_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{either}{$rare_het_key}++;
										}
									}
									
									if (defined $filter_dbsnp_indel_str){
										if ($filter_dbsnp_indel_str eq 'NOVEL' || $filter_dbsnp_indel_str eq 'RARE_ALLELE' || $filter_dbsnp_indel_str eq 'NO_FREQ'){
											$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{either}{$rare_het_key}++;
										}
									}

								}
							}
							
							
						}
						
						if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 0 && $sz ne 'ref' && $sz ne '?') {
							$affected_alleles_from_one_parent = 'no';
						}
						
					} 
					
					#delete het entries if any pre-conditions are broken
					if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 0) {
						#unaffected can't be homozygous
						if ($sz eq 'hom') {
							$del_com_het = 1;
						}
					} elsif ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
						#all affected must be hets
						if ($sz ne 'het') {
							$del_com_het = 1;
						}
					}
					
					
					if ($self->_has_father($sg_name)) {
						$self->{source_info}{$sg_name}{allele_data}{$var_type}{$var_key}{$sample}{father} = $father_allele;
						if ($father_allele eq '?') {
							$affected_alleles_from_one_parent = '?';
						}
					} 
					
					if ($self->_has_mother($sg_name)) {
						$self->{source_info}{$sg_name}{allele_data}{$var_type}{$var_key}{$sample}{mother} = $mother_allele;
						if ($mother_allele eq '?') {
							$affected_alleles_from_one_parent = '?';
						}
					} 
								
				
					
				}	
				
				#If we broke comhet conditions for variant we need to remove existing entries
				if (exists $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}) {
					if ($del_com_het) {
					
						for my $sample (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}}) {
							for my $het_parent (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}}) {
								#print "deleting com het for $sample $het_parent $chr $coord $var_type \n";						
								my $gmaf = $var_type eq 'snv'?$gmaf_snv_str:$gmaf_indel_str;
								if (exists $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{$het_parent}{"$chr:$coord:$var_type:gmaf=$gmaf"}) {
									if (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{$het_parent}} > 1) {
										#If other entries exist only remove coordinate
										delete $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{$het_parent}{"$chr:$coord:$var_type:gmaf=$gmaf"};																								
									} else {
										#If last entry remove the het_parent field as well
										delete $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}{$het_parent};
									}
								}
								
								if (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets}{$sample}} == 0) {
									#If nothing left for gene remove entire entry
									delete $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{hets};
								}							
							}
							
						}
					} 
				}
							
				if (exists $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}) {
					if ($del_com_het) {
					
						for my $sample (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}}) {
							for my $het_parent (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}}) {
								#print "deleting com het for $sample $het_parent $chr $coord $var_type \n";						
								my $gmaf = $var_type eq 'snv'?$gmaf_snv_str:$gmaf_indel_str;
								my $filter_dbsnp = $var_type eq 'snv'?$filter_dbsnp_snv_str:$filter_dbsnp_indel_str;

								if (exists $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{$het_parent}{"$chr:$coord:$var_type:gmaf=$gmaf:filter_dbsnp=$filter_dbsnp"}) {
									if (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{$het_parent}} > 1) {
										#If other entries exist only remove coordinate
										delete $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{$het_parent}{"$chr:$coord:$var_type:gmaf=$gmaf:filter_dbsnp=$filter_dbsnp"};																								
									} else {
										#If last entry remove the het_parent field as well
										delete $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}{$het_parent};
									}
								}														
								if (keys %{$self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets}{$sample}} == 0) {
									#If nothing left for gene remove entire entry
									delete $self->{source_info}{$sg_name}{gene_counts}{$self->{source_info}{$sg_name}{gene_lookup}{$var_key}}{rare_hets};
								}							
							}							
						}
					} 
				}
								
				if ($affected_alleles_from_one_parent eq 'yes' && keys %affected_allele_data == 1) {
					if (exists $affected_allele_data{mother}) {
						$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} = "Mother gives affected allele";
					} elsif (exists $affected_allele_data{father}) {
						$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} = "Father gives affected allele";
					} else {
						modules::Exception->throw("ERROR: Incorrect key for affected single allele");
					}
				} elsif ($affected_alleles_from_one_parent eq '?') {
					$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} = '?';
				} else {
					$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} = 'No';
				}
			}
		}
	}
	
}

#Get phasing blocks
sub get_phasing_blocks {
	my ($self) = @_;
	
	
	my %variant_data = %{$self->{source_info}};

	#flags for detecting blocks
	my $mother_block_start = 0;
	my $father_block_start = 0;
	my $mother_block_end = 0;
	my $father_block_end = 0;
	my $current_var_count = 0;
	my %block_coords = ();

	#params for comparison
	my $min_variant_num = 2;
	my $min_block = 10000;
	

	if ($self->{filter_output}) {
		if ($self->{filters}{phase_block_size}) {
			#overwrite if necessary
			$min_block = $self->{filters}{phase_block_size};			
		}
		if ($self->{filters}{phase_var_num}) {
			#overwrite if necessary
			$min_variant_num = $self->{filters}{phase_var_num};			
		}
	} 

	for my $sg_name (keys %variant_data) {
		
		next unless $self->_has_single_parent($sg_name);
		
		for my $var_type (keys %{$self->{source_info}{$sg_name}{var_data}}) {
			
		
			my $prev_chr = 'N';
			for my $var_key (sort {my ($a_chr,$a_coord) = split(':',$a); my ($b_chr,$b_coord) = split(':',$b); $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %{$variant_data{$sg_name}{var_data}{$var_type}}) {
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{phase_info} = 'No'; #Set the default
				my ($chr,$coord) = split(':',$var_key);
				if ($prev_chr ne $chr) {
					if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
						if ($mother_block_start != 0) {
							my $block_count = keys %block_coords;
							my $size = $mother_block_end - $mother_block_start;
							for my $block_key (keys %block_coords) {
								$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Mother; ' . $chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
							}
						} elsif ($father_block_start != 0) {
							my $block_count = keys %block_coords;
							my $size = $father_block_end - $father_block_start;
							for my $block_key (keys %block_coords) {
								$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Father; ' . $chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
							}							
						} else {
							modules::Exception->throw("ERROR: Can't find block");
						}
					}
					$prev_chr = $chr;
					$current_var_count = 0;
					$father_block_start = 0;
					$father_block_end = 0;
					$mother_block_start = 0;
					$mother_block_end = 0;
					%block_coords = ();
					
				}
				
				if ($self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} =~ /Mother/) {
					#Mother allele
					if ($father_block_start != 0 ) {
						#Check if we're coming out of father block
						if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
							my $block_count = keys %block_coords;
							my $size = $father_block_end - $father_block_start;
							for my $block_key (keys %block_coords) {
								$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Father; ' . $chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
							}	
						}
						
						$current_var_count = 0;
						$father_block_start = 0;
						$father_block_end = 0;
						%block_coords = ();
					}
					
					if ($mother_block_start == 0) {
						$mother_block_start = $coord;
					}
					$mother_block_end = $coord;
					$block_coords{$var_key}++;
					$current_var_count++;
					
				} elsif ($self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} =~ /Father/) {
					#Father_allele found
					if ($mother_block_start != 0) {
						#Check if we're coming out of mother block
						if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
							my $block_count = keys %block_coords;
							my $size = $mother_block_end - $mother_block_start;
							for my $block_key (keys %block_coords) {
								$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Mother; ' . $chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
							}
						}
						
						$current_var_count = 0;
						$mother_block_start = 0;
						$mother_block_end = 0;
						%block_coords = ();
					}
		
					if ($father_block_start == 0) {
						$father_block_start = $coord;
					}			
					$father_block_end = $coord;
					$block_coords{$var_key}++;
					$current_var_count++;
					
				} elsif ($self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected} eq 'No') {
					#Trigger to check for end of a block
					if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
						if ($mother_block_start != 0) {
							my $block_count = keys %block_coords;
							my $size = $mother_block_end - $mother_block_start;
							for my $block_key (keys %block_coords) {
								$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Mother; ' . $chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
							}
						} elsif ($father_block_start != 0) {
							my $block_count = keys %block_coords;
							my $size = $father_block_end - $father_block_start;
							for my $block_key (keys %block_coords) {
								$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Father; ' . $chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
							}	
						} else {
							modules::Exception->throw("ERROR: Can't find block");
						}
					}					
					%block_coords = ();
					$current_var_count = 0;
					$father_block_start = 0;
					$father_block_end = 0;
					$mother_block_start = 0;
					$mother_block_end = 0;
				}
			}
		
			#Check final case if it's at end of last chromosome	
			if ($self->_check_block($mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block)) {
				if ($mother_block_start != 0) {
					my $block_count = keys %block_coords;
					my $size = $mother_block_end - $mother_block_start;
					for my $block_key (keys %block_coords) {
						$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Mother; ' . $prev_chr .':'.$mother_block_start.'-'.$mother_block_end.'; '.$size .'bp';
					}
				} elsif ($father_block_start != 0) {
					my $block_count = keys %block_coords;
					my $size = $father_block_end - $father_block_start;
					for my $block_key (keys %block_coords) {
						$self->{source_info}{$sg_name}{var_data}{$var_type}{$block_key}{phase_info} = $block_count . ' variants; Father; ' . $prev_chr .':'.$father_block_start.'-'.$father_block_end.'; '.$size .'bp';
					}	
				} else {
					modules::Exception->throw("ERROR: Can't find block");
				}			
			}
						
		}
		
	}
}




#Get the compound het info and update disease_inheritance field
sub get_compound_het {
	my ($self) = @_;
	
	my %variant_data = %{$self->{source_info}};
	
	for my $sg_name (keys %variant_data) {
		
		next unless $self->_has_single_parent($sg_name);
		
		for my $var_type (keys %{$variant_data{$sg_name}{var_data}}) {
			print "Parsing $var_type..\n";
			my $prev_chr = 'N';
			for my $var_key (sort {my ($a_chr,$a_coord) = split(':',$a); my ($b_chr,$b_coord) = split(':',$b); $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %{$variant_data{$sg_name}{var_data}{$var_type}}) {
				my @line_data;
				my ($chr,$coord) = split(':',$var_key);
				if ($prev_chr ne $chr) {
					print "Chr $chr\n" if $var_type eq 'snv'; #Report progress as this step takes a while
					$prev_chr = $chr;
				}
				#my $aff_allele_count = 0;
				#my $aff_gene_count = 0;
				my %gene_counts = %{$self->{source_info}{$sg_name}{gene_counts}};
				my %gene_lookup = %{$self->{source_info}{$sg_name}{gene_lookup}};
				
#				if (exists $gene_counts{$gene_lookup{$var_key}}{uniq_aff}) {
#					$aff_allele_count = $gene_counts{$gene_lookup{$var_key}}{uniq_aff}{$var_key} if (exists $gene_counts{$gene_lookup{$var_key}}{uniq_aff} && $gene_counts{$gene_lookup{$var_key}}{uniq_aff}{$var_key});
#					$aff_gene_count = keys %{$gene_counts{$gene_lookup{$var_key}}{uniq_aff}};
#				}
				
				#my $unaff_allele_count = 0;
	
				
#				if (exists $gene_counts{$gene_lookup{$var_key}}{uniq_unaff}) {
#					$unaff_allele_count = $gene_counts{$gene_lookup{$var_key}}{uniq_unaff}{$var_key} if (exists $gene_counts{$gene_lookup{$var_key}}{uniq_unaff} && $gene_counts{$gene_lookup{$var_key}}{uniq_unaff}{$var_key});
#				}
#				
				my $def_compound_het = 0;
				my $pos_compound_het = 0;
				my $def_compound_rare_het = 0;
				my $pos_compound_rare_het = 0;
				
				my %compound_hets = ();
				if (exists $gene_counts{$gene_lookup{$var_key}}{hets}) {

					#Here we have definite potential compound het conditions for this sample; check that this variant is in the list
					for my $sample (keys %{$gene_counts{$gene_lookup{$var_key}}{hets}}) {
					
						if (exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{mother} && exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{father}) {
							for my $parent (qw(mother father)) {
								for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{$parent}})  {
									my ($local_chr,$local_coord) = split(':',$het_str);
									if ($local_chr eq $chr && $local_coord == $coord) {
										$def_compound_het = 1;
									}
								}
							}
							
							
							next unless $def_compound_het;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_het_str = '(M='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{mother}}).',F='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{father}}).')';
							push @{$compound_hets{def}{$compound_het_str}}, $short_sample;
						} 
					}
					
					
					#Now check for possible compound hets
					for my $sample (keys %{$gene_counts{$gene_lookup{$var_key}}{hets}}) {
						if (exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{mother} && exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}) {
							#Here we have mother and ambiguous cases
							my $found_pos = 0;
							for my $parent (qw(mother either)) {
								for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{$parent}})  {
									my ($local_chr,$local_coord) = split(':',$het_str);
									if ($local_chr eq $chr && $local_coord == $coord) {
										$pos_compound_het = 1;
										$found_pos = 1;	
									}
								}
							}
							next unless $found_pos;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_het_str = '(M='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{mother}}).',AMB='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}}).')';
							push @{$compound_hets{pos}{$compound_het_str}}, $short_sample;
						} elsif (exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{father} && exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}) {
							#Here we have father and ambiguous cases
							my $found_pos = 0;
							for my $parent (qw(father either)) {
								for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{$parent}})  {
									my ($local_chr,$local_coord) = split(':',$het_str);
									if ($local_chr eq $chr && $local_coord == $coord) {
										$found_pos = 1;
										$pos_compound_het = 1;
									}
								}
							}
							next unless $found_pos;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_het_str = '(F='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{father}}).',AMB='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}}).')';
							push @{$compound_hets{pos}{$compound_het_str}}, $short_sample;
						} elsif (exists $gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either} && keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}} > 1) { 
							#Here we have at least 2 ambiguous cases
							my $found_pos = 0;
							for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}})  {
								my ($local_chr,$local_coord) = split(':',$het_str);
								if ($local_chr eq $chr && $local_coord == $coord) {
									$pos_compound_het = 1;
									$found_pos = 1;
								}
							}
							next unless $found_pos;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_het_str = '(AMB='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{hets}{$sample}{either}}).')';
							push @{$compound_hets{pos}{$compound_het_str}}, $short_sample;	
						
						}
						
					}
				}
				
				if (exists $gene_counts{$gene_lookup{$var_key}}{rare_hets}) {

					#Here we have definite potential compound het conditions for this sample; check that this variant is in the list
					for my $sample (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}}) {
						if (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{mother}} && keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{father}}) {
							for my $parent (qw(mother father)) {
								for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{$parent}})  {
									my ($local_chr,$local_coord) = split(':',$het_str);
									if ($local_chr eq $chr && $local_coord == $coord) {
										$def_compound_rare_het = 1;
									}
								}
							}
							
							
							next unless $def_compound_rare_het;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_rare_het_str = '(M='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{mother}}).',F='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{father}}).')';
							push @{$compound_hets{rare_def}{$compound_rare_het_str}}, $short_sample;
						} 
					}
					
					
					#Now check for possible compound hets
					for my $sample (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}}) {
						if (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{mother}} && keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}}) {
							#Here we have mother and ambiguous cases
							my $found_pos = 0;
							for my $parent (qw(mother either)) {
								for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{$parent}})  {
									my ($local_chr,$local_coord) = split(':',$het_str);
									if ($local_chr eq $chr && $local_coord == $coord) {
										$pos_compound_rare_het = 1;
										$found_pos = 1;	
									}
								}
							}
							next unless $found_pos;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_rare_het_str = '(M='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{mother}}).',AMB='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}}).')';
							push @{$compound_hets{rare_pos}{$compound_rare_het_str}}, $short_sample;
						} elsif (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{father}} && keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}}) {
							#Here we have father and ambiguous cases
							my $found_pos = 0;
							for my $parent (qw(father either)) {
								for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{$parent}})  {
									my ($local_chr,$local_coord) = split(':',$het_str);
									if ($local_chr eq $chr && $local_coord == $coord) {
										$found_pos = 1;
										$pos_compound_rare_het = 1;
									}
								}
							}
							next unless $found_pos;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_rare_het_str = '(F='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{father}}).',AMB='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}}).')';
							push @{$compound_hets{rare_pos}{$compound_rare_het_str}}, $short_sample;
						} elsif (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}} && keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}} > 1) { 
							#Here we have at least 2 ambiguous cases
							my $found_pos = 0;
							for my $het_str (keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}})  {
								my ($local_chr,$local_coord) = split(':',$het_str);
								if ($local_chr eq $chr && $local_coord == $coord) {
									$pos_compound_rare_het = 1;
									$found_pos = 1;
								}
							}
							next unless $found_pos;
							my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
							my $compound_rare_het_str = '(AMB='.join(",",keys %{$gene_counts{$gene_lookup{$var_key}}{rare_hets}{$sample}{either}}).')';
							push @{$compound_hets{rare_pos}{$compound_rare_het_str}}, $short_sample;	
						
						}
						
					}
				}	
	
				my $disease_inheritance =  $variant_data{$sg_name}{var_data}{$var_type}{$var_key}{disease_inheritance};
				my $def_com_het;	
				my $def_com_rare_het;
				my $pos_com_het;
				my $pos_com_rare_het;
			
				if (!$def_compound_het) {
					$def_com_het = 'No';
				} else {
					#Check that unaffected don't have exact same comhets
					my $diff_unaffected = 1;
					
					if (keys %{$compound_hets{def}} == 1) {
						for my $str (keys %{$compound_hets{def}}) {
							my $sample_string = join(",",@{$compound_hets{def}{$str}});
							if ($sample_string =~ /unaffected/) {
								$diff_unaffected = 0;
							}
						}
					}
					
					if ($diff_unaffected) {
					
						for my $str (keys %{$compound_hets{def}}) {
							$def_com_het .= join(",",@{$compound_hets{def}{$str}}) . $str . ',';
						}
						$def_com_het =~ s/,$//;
						if ($disease_inheritance eq 'none') {
							$disease_inheritance = "def_compound_het";
						} else {
							$disease_inheritance .= ",def_compound_het";
						}
					} else {
						$def_com_het = 'No';
					}
				}
				
				if (!$pos_compound_het) {
					$pos_com_het = 'No';
				} else {
					
					#print Dumper $compound_hets{pos};
					my $diff_unaffected = 1;
					
					if (keys %{$compound_hets{pos}} == 1) {
						for my $str (keys %{$compound_hets{pos}}) {
							my $sample_string = join(",",@{$compound_hets{pos}{$str}});
							if ($sample_string =~ /unaffected/) {
								$diff_unaffected = 0;
							}
						}
					}
					
					if ($diff_unaffected) {
						for my $str (keys %{$compound_hets{pos}}) {
							$pos_com_het .= join(",",@{$compound_hets{pos}{$str}}) . $str . ',';
						}
						$pos_com_het =~ s/,$//;
						if ($disease_inheritance eq 'none') {
							$disease_inheritance = "pos_compound_het";
						} else {
							$disease_inheritance .= ",pos_compound_het";
						}
						
					} else {
						$pos_com_het = 'No';
					}
				}	

				###### get $def_com_rare_het and $pos_com_rare_het
		
				if (!$def_compound_rare_het) {
					$def_com_rare_het = 'No';
				} else {
					#Check that unaffected don't have exact same comhets
					my $diff_unaffected = 1;
					
					if (keys %{$compound_hets{rare_def}} == 1) {
						for my $str (keys %{$compound_hets{rare_def}}) {
							my $sample_string = join(",",@{$compound_hets{rare_def}{$str}});
							if ($sample_string =~ /unaffected/) {
								$diff_unaffected = 0;
							}
						}
					}
					
					if ($diff_unaffected) {
					
						for my $str (keys %{$compound_hets{rare_def}}) {
							$def_com_rare_het .= join(",",@{$compound_hets{rare_def}{$str}}) . $str . ',';
						}
						$def_com_rare_het =~ s/,$//;
						if ($disease_inheritance eq 'none') {
							$disease_inheritance = "def_compound_rare_het";
						} else {
							$disease_inheritance .= ",def_compound_rare_het";
						}
					} else {
						$def_com_rare_het = 'No';
					}
				}
				
				if (!$pos_compound_rare_het) {
					$pos_com_rare_het = 'No';
				} else {
					
					#print Dumper $compound_hets{pos};
					my $diff_unaffected = 1;
					
					if (keys %{$compound_hets{rare_pos}} == 1) {
						for my $str (keys %{$compound_hets{rare_pos}}) {
							my $sample_string = join(",",@{$compound_hets{rare_pos}{$str}});
							if ($sample_string =~ /unaffected/) {
								$diff_unaffected = 0;
							}
						}
					}
					
					if ($diff_unaffected) {
						for my $str (keys %{$compound_hets{rare_pos}}) {
							$pos_com_rare_het .= join(",",@{$compound_hets{rare_pos}{$str}}) . $str . ',';
						}
						$pos_com_rare_het =~ s/,$//;
						if ($disease_inheritance eq 'none') {
							$disease_inheritance = "pos_compound_rare_het";
						} else {
							$disease_inheritance .= ",pos_compound_rare_het";
						}
						
					} else {
						$pos_com_rare_het = 'No';
					}
				}				
				#Get the mother and father alleles first (per sample)
				my %allele_data = ();
				my $mother_allele_str;
				my $father_allele_str;
				
				
				for my $sample (sort keys %{$variant_data{$sg_name}{allele_data}{$var_type}{$var_key}}) {
	
					my ($short_sample) = $sample =~ /_([a-z]+\d+)$/;
					
					if ($self->_has_father($sg_name)) {
						my $father_allele = $self->{source_info}{$sg_name}{allele_data}{$var_type}{$var_key}{$sample}{father};
						$father_allele_str .= $short_sample.' ('.$father_allele.'),';
					}
					if ($self->_has_mother($sg_name)) {
						my $mother_allele = $self->{source_info}{$sg_name}{allele_data}{$var_type}{$var_key}{$sample}{mother};
						$mother_allele_str .= $short_sample.' ('.$mother_allele.'),';
					}
					
				}
				
				
				
				if ($self->_has_mother($sg_name) && $mother_allele_str) {
					$mother_allele_str =~ s/,$//;
					$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{mother_allele_str} = $mother_allele_str;
				} elsif ($self->_has_mother($sg_name)) {
					print Dumper $variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data};
					print Dumper $self->{source_info}{$sg_name}{allele_data}{$var_type}{$var_key};
					modules::Exception->throw("ERROR: No allele info for $var_key and mother_allele_str $mother_allele_str\n");
				}
				
				if ($self->_has_father($sg_name) && $father_allele_str) {
					$father_allele_str =~ s/,$//;
					$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{father_allele_str} = $father_allele_str;		
				} elsif ($self->_has_father($sg_name)) {
					print Dumper $variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data};
					print Dumper $self->{source_info}{$sg_name}{allele_data}{$var_type}{$var_key};
					modules::Exception->throw("ERROR: No allele info for $var_key and father_allele_str $father_allele_str\\n");
				}
				
				#push @line_data,  $likely_aff_allele_count . ' ('.$aff_allele_count.')', $likely_unaff_allele_count . ' ('.$unaff_allele_count.')', $mendel_inheritance, $disease_inheritance, $gene, $gene_total, $gene_uniq,  $mother_allele_pos_compound_rare_hetr_allele_str,$def_com_het, $pos_com_het;
	
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{disease_inheritance} = $disease_inheritance;
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{def_com_het} = $def_com_het;
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{pos_com_het} = $pos_com_het;
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{def_com_rare_het} = $def_com_rare_het;
				$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{pos_com_rare_het} = $pos_com_rare_het;
			}
		}
	}
}






#Generate the full summary lines
sub generate_line_data {
	my ($self) = @_;

	my %mutant_snv_data = ();
	my %mutant_indel_data = ();
		
	my %variant_data = %{$self->{source_info}};
	
	for my $sg_name (keys %variant_data) {
		#These headers are needed for snvs and indels
		my @new_headers = (@{$self->{source_info}{$sg_name}{extra_headers}}, @{$self->{source_info}{$sg_name}{sample_headers}});
		
		
		my @all_snv_headers = @new_headers;
		my @all_indel_headers = @new_headers;	
		
		my %snv_headers = %{$self->{source_info}{$sg_name}{snv_headers}};
		my %indel_headers = %{$self->{source_info}{$sg_name}{indel_headers}};
		
		#Add the variant type specific headers
		for my $snv_header (sort {$snv_headers{$a} <=> $snv_headers{$b}} keys %snv_headers) {
			push @all_snv_headers, $snv_header;
		}
		
		for my $indel_header (sort {$indel_headers{$a} <=> $indel_headers{$b}} keys %indel_headers) {
			push @all_indel_headers, $indel_header;
		}
		
		my %unsorted_lines = ();
		$self->{source_info}{$sg_name}{snv_headers} =  join("\t",@all_snv_headers);
		$self->{source_info}{$sg_name}{indel_headers} = join("\t",@all_indel_headers);
		
		if ($self->{score}){
			$self->{source_info}{$sg_name}{snv_headers} .= "\t"."variant_score_total"."\t"."gene_filter";
			$self->{source_info}{$sg_name}{indel_headers} .= "\t"."variant_score_total"."\t"."gene_filter";
		}
		
		for my $var_type (keys %{$variant_data{$sg_name}{var_data}}) {
			print "Parsing $var_type..\n";
			
			my $var_type_xml;
			
			if( $var_type eq 'DEL' || $var_type eq 'INS'){
				$var_type_xml = 'indel';
			}else{
				$var_type_xml = 'snv';
			}
			
			my $proband_samplename;
	    	my @phenotype_keywords ;
	    	my @gene_filter ;
	        my @exon_overlap ;
	        my %dbsnp_freq_score ;
	        my $polyphen_score_cutoff;
	        my $sift_score_cutoff;
	        my $cadd_score_cutoff;
	        my $filter_exac_freq_cutoff;
	        my $final_status;
	                    	
			if ($self->{family} && $self->{score}){
				
				for my $sample ( keys %{$variant_data{$sg_name}{samples}} ) {
					if($variant_data{$sg_name}{samples}{$sample}{relation} eq 'proband'){
						$proband_samplename = $sample;
					}
				}
				
				#Adding scores for CPI samples
				my $config = modules::Pipeline->get_report_conf(-report_xml=>$self->{report_xml});					
	    
	            # get score calculated for variants with these conditions: 
	            # 1) present in proband (affected) AND in certain proportion of affected group
				# 2) match exon_types (@exon_overlap) in report config
				# 3) match $final_status in report config
				# 4) meet filter_exac_freq_cutoff (or NA) as in report config
						
				if ($config->exists($var_type_xml,'common','scoring_filters','keywords')) {
	    			@phenotype_keywords = split(",",$config->read($var_type_xml,'common','scoring_filters','keywords'));
	    		}
				if ($config->exists($var_type_xml,'common','scoring_filters','gene_filter')) {
	    			@gene_filter = split(",",$config->read($var_type_xml,'common','scoring_filters','gene_filter'));
	    		}
	                  
	            if ($config->exists($var_type_xml,'common','scoring_filters','exon_overlap')) {
	                @exon_overlap = split(",",$config->read($var_type_xml,'common','scoring_filters','exon_overlap'));
	            }
	                    
	            if ($config->exists($var_type_xml,'common','scoring_filters','dbfreq_score')) {
	                my @dbfreq_score = split(",",$config->read($var_type_xml,'common','scoring_filters','dbfreq_score'));
	                for my $dbfreq_bin (@dbfreq_score) {
	                 	my @fields = split(":",$dbfreq_bin);
	                    $dbsnp_freq_score{$fields[0]} = $fields[1];
	                }
	            }
	                    
	            if ($config->exists($var_type_xml,'common','scoring_filters','polyphen_score')) {
	                $polyphen_score_cutoff = $config->read($var_type_xml,'common','scoring_filters','polyphen_score');
	            }
	                   
	            if ($config->exists($var_type_xml,'common','scoring_filters','sift_score')) {
	                $sift_score_cutoff = $config->read($var_type_xml,'common','scoring_filters','sift_score');
	            }
	
	            if ($config->exists($var_type_xml,'common','scoring_filters','cadd_score')) {
	                $cadd_score_cutoff = $config->read($var_type_xml,'common','scoring_filters','cadd_score');
	            }
	                  
	            if ($config->exists($var_type_xml,'common','scoring_filters','filter_exac')) {
	                $filter_exac_freq_cutoff = $config->read($var_type_xml,'common','scoring_filters','filter_exac');
	            }
	                    
	            if ($config->exists($var_type_xml,'common','scoring_filters','final_status')) {
	                $final_status = $config->read($var_type_xml,'common','scoring_filters','final_status');
	            }
			}				
			
			my $prev_chr = 'N';
			for my $var_key (sort {my ($a_chr,$a_coord) = split(':',$a); my ($b_chr,$b_coord) = split(':',$b); $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %{$variant_data{$sg_name}{var_data}{$var_type}}) {
				my @line_data;
				my ($chr,$coord,$end,$var_base) = split(':',$var_key);
				if ($prev_chr ne $chr) {
					print "Chr $chr\n" if $var_type eq 'snv'; #Report progress as this step takes a while
					$prev_chr = $chr;
				}
				
				
				my %gene_counts = %{$self->{source_info}{$sg_name}{gene_counts}};
				my %gene_lookup = %{$self->{source_info}{$sg_name}{gene_lookup}};
				my %gene_lengths = %{$self->{source_info}{$sg_name}{gene_length}}; 
				
				my $unaff_allele_count = 0;
	
				#Now get a count that reflects the zygosities you see in the pileups
				my $likely_aff_allele_count = 0;
				my $likely_unaff_allele_count = 0;
				my %var_samples = ();
				#Lookup for pileup
				my $pileup_key;
				if ($var_type eq 'DEL') {
					my $indel_coord = $coord - 1;
					$pileup_key = "$chr:$indel_coord";
				} else {
					$pileup_key = "$chr:$coord";
				}
				
				if ($self->{pileups}) {
					for my $sample (sort keys %{$self->{source_info}{$sg_name}{samples}}) {
						
						if ($self->{pileups}) {
							#Use pileups for determining aff/unaff if available
							if (exists $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample} && exists $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg} && ($self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg} eq 'het' || $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg} eq 'hom')) {
								$var_samples{$sample}++; #build up sample list when we have pileups
								if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
									$likely_aff_allele_count++;
								} else {
									$likely_unaff_allele_count++;
								}
							} 
						} 
					}
				} else {
				
					#Otherwise we need to use the report info to determine if variant is present
#					if (exists $self->{source_info}{$sg_name}{var_counts}{$var_type}{$var_key}{aff}) {
#						$likely_aff_allele_count = $self->{source_info}{$sg_name}{var_counts}{$var_type}{$var_key}{aff};
#					} 
#					if (exists $self->{source_info}{$sg_name}{var_counts}{$var_type}{$var_key}{unaff}) {
#						$likely_unaff_allele_count = $self->{source_info}{$sg_name}{var_counts}{$var_type}{$var_key}{unaff};
#					}
					
					
					#Do it this way to avoid double counting 'COMBINED' cases 
					for my $sample (keys %{$self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{data}}) {
						if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
							$likely_aff_allele_count++;
						} else {
							$likely_unaff_allele_count++;
						}	
					}
					
				}
				
				my $sample_str;
				if ($self->{pileups}) {
					if (keys  %var_samples) {
						$sample_str =  join(",",sort keys %var_samples)
					} else {
						$sample_str = 	join(',',sort keys %{$variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}}) . ' (NOT IN PILEUP)';					
					}
				} else {
					$sample_str = join(',',sort keys %{$variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}});
					
				}
				
				if (! defined $likely_aff_allele_count && ! defined $likely_unaff_allele_count) {
					modules::Exception->throw("ERROR: Can't retrieve aff/unaff count for $var_key variant\n");
				}
				
				
				if (exists $gene_counts{$gene_lookup{$var_key}}{uniq_unaff}) {
					$unaff_allele_count = $gene_counts{$gene_lookup{$var_key}}{uniq_unaff}{$var_key} if (exists $gene_counts{$gene_lookup{$var_key}}{uniq_unaff} && $gene_counts{$gene_lookup{$var_key}}{uniq_unaff}{$var_key});
				}
				
				my $gene = $gene_lookup{$var_key};
				my $gene_uniq = keys  %{$gene_counts{$gene}{uniq_total}};
				my $gene_total = $gene_counts{$gene}{non_uniq_total};
				my $gene_var_percent;
				if ($gene_lengths{$gene} =~ /\d/) {
					$gene_var_percent =  sprintf("%.2f",($gene_uniq / ($gene_lengths{$gene} * 3)) * 100) . "%";	
				} else {
					$gene_var_percent = "N/A";
				}
				
				my $aff_total = defined $self->{source_info}{$sg_name}{aff_count}?$self->{source_info}{$sg_name}{aff_count}:0;
				my $unaff_total = defined $self->{source_info}{$sg_name}{unaff_count}?$self->{source_info}{$sg_name}{unaff_count}:0;
				
				if ($self->{filter_output}) {
					if ($self->{filters}{min_num_aff}) {
						next unless $likely_aff_allele_count >= $self->{filters}{min_num_aff};
					}
					if ($self->{filters}{inheritance}) {
						next unless $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{disease_inheritance} eq $self->{filters}{inheritance};
					}
					if ($self->{filters}{denovo}) {
						next unless $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{mendel_inheritance} =~ /de novo/;
					}
					if ($self->{filters}{comhet}) {
						next unless $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{def_com_het} ne 'No';
					}
					
				}
				
				
				push @line_data, $sample_str;
				
				if ($self->{source_info}{$sg_name}{affected}) {
					push @line_data,  $likely_aff_allele_count . ' out of '.$aff_total;
				} 
				if ($self->{source_info}{$sg_name}{unaffected}) {
					push @line_data,  $likely_unaff_allele_count . ' out of '.$unaff_total;
				}
				
				my $mendelian = 1;
				#Get the precalculated values				
				if ($self->{family}) {
					my $mendel_inheritance = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{mendel_inheritance};
					
					$mendelian = $mendel_inheritance eq 'mendelian_rules_followed'?1:0; #Overwite if necessary
					my $def_com_het = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{def_com_het};
					my $pos_com_het = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{pos_com_het};
					my $def_com_rare_het = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{def_com_rare_het};
					my $pos_com_rare_het = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{pos_com_rare_het};
										
					my $disease_inheritance =  $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{disease_inheritance};
					
					
					
					if ($self->_has_single_parent($sg_name)) {
						push @line_data, $mendel_inheritance;
					}
					#if ($self->{family} && $self->{unaffected}) {
					push @line_data, $disease_inheritance;
					#}
					push @line_data, $gene, $gene_total, $gene_uniq, $gene_var_percent;
					
					if ($self->_has_two_parents($sg_name)) {
						my $mother_allele_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{mother_allele_str};
						my $father_allele_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{father_allele_str};
						my $parent_allele_common_to_affected = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected};
						my $phase_block = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{phase_info};
						push @line_data, $mother_allele_str, $father_allele_str, $parent_allele_common_to_affected, $phase_block, $def_com_het, $pos_com_het, $def_com_rare_het, $pos_com_rare_het;
					} elsif ($self->_has_mother($sg_name)) {
						my $mother_allele_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{mother_allele_str};
						my $parent_allele_common_to_affected = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected};
						my $phase_block = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{phase_info};
						push @line_data, $mother_allele_str, $parent_allele_common_to_affected, $phase_block, $def_com_het, $pos_com_het, $def_com_rare_het, $pos_com_rare_het;
					} elsif ($self->_has_father($sg_name)) {
						my $father_allele_str = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{father_allele_str};
						my $parent_allele_common_to_affected = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{parent_allele_affected};
						my $phase_block = $self->{source_info}{$sg_name}{var_data}{$var_type}{$var_key}{phase_info};
						push @line_data, $father_allele_str, $parent_allele_common_to_affected, $phase_block, $def_com_het, $pos_com_het, $def_com_rare_het, $pos_com_rare_het;
					}  
					
				} else {
					push @line_data, $gene, $gene_total, $gene_uniq, $gene_var_percent;
				}
				
				
				
				
				
				#Now get the pileup data
				if ($self->{pileups}) {
					for my $sample (sort keys %{$self->{source_info}{$sg_name}{samples}}) {
						push @line_data, $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup} . ' ('.$self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg}.')';	
					}
				}
				 
				if ($var_type eq 'snv') {
					push @line_data, $chr, $coord;
				} else {
					my ($chr,$start,$end) = split(":",$var_key);
					push @line_data, $chr, $start, $end;
				}
				
				my $pass = 0;
				my $dbsnp_freq;
				my %values = ();
				for my $sample (sort keys %{$variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}}) {
					for my $col_str (sort {my ($a_index) = split(':',$a); my ($b_index) = split(':',$b); $a_index <=> $b_index} keys %{$variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}}) {
						my ($index,$col_name) = split(':',$col_str);
						if ($col_name eq 'dbsnp_var_allele_freq') {
							if ($variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{$col_str} =~ /\d/) {
								$dbsnp_freq = $variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{$col_str};
							} else {
								$dbsnp_freq = 1;
							}
						}
						#next if $col_name eq 'alleles';
						$pass = 1 if $variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{$col_str} =~ /OVERALL_PASS/;
						$values{$index}{$variant_data{$sg_name}{var_data}{$var_type}{$var_key}{data}{$sample}{$col_str}}++;
						#push @line_data, join(';', $variant_data{$sg_name}{$var_type}{$var_key}{$sample}{data}{$col_str});
					}
				}
				
				#Add no_data for entries missing data like clr_score
				my $index_counter = 0;		
		
				for my $index (sort {$a<=>$b} keys %values) {
					if ($index != $index_counter) {
						while ($index_counter < $index) {
							push @line_data, "No Data" unless $index <= 3; #account for chr and coord str 
							$index_counter++;
						}
					} 
					
					my $value_str = join('; ', keys %{$values{$index}});
					
					push @line_data, $value_str;
					$index_counter++;
				}
				
				my $line = join("\t",@line_data);
	
				if ($var_type eq 'snv') {

						my $pass_exon_overlap = 0;
						my $pass_status = 0;
						my $pass_filter_exac = 0;   
						my $freqScore = 0;
						my $damagingScore = 0;
						my $goScore = 0;
						my $phenotypeScore = 0; 
						my $dbsnp_freq;
						my $sift_score;
						my $polyphen_score;
						my $cadd_phred;
						my $final_score = 0;
						my $match_gene_filter = 0;
												
						for (my $col_count = 0 ; $col_count < @all_snv_headers ; $col_count++ ) {
							my $tmp_col = $all_snv_headers[$col_count];
							$tmp_col =~ s/^%/percent/;
							$tmp_col =~ s/\(/_/g;
							$tmp_col =~ s/\)/_/g;
							$tmp_col =~ s/ //g;
							$tmp_col =~ s/\//_of_/;
							$mutant_snv_data{"$chr:$coord"}{$end}{$var_base}{$tmp_col} = $line_data[$col_count];
					                         		
							if($self->{score} && $tmp_col eq 'snv_exon_type'){
								for my $exon_type (@exon_overlap){
									if($line_data[$col_count] eq $exon_type){
										$pass_exon_overlap = 1;
									}
								}
							}
							if($self->{score} && $tmp_col eq 'gene'){
								for my $gene (@gene_filter){
									if($line_data[$col_count] eq $gene){
										$match_gene_filter = 1;
									}
								}
							}
							if($self->{score} && $tmp_col eq 'final_status'){
								if($line_data[$col_count] eq $final_status){
									$pass_status = 1;
								}
							}
							if($self->{score} && $tmp_col eq 'filter_exac_snv'){
								if($line_data[$col_count] eq 'N/A'){
									$pass_filter_exac = 1;
								}
								else{
									my ($a,$freq) = split(":",$line_data[$col_count]);
									if ($freq < $filter_exac_freq_cutoff){
										$pass_filter_exac = 1;
									}
								}
							}
							if($self->{score} && $tmp_col eq 'dbsnp_var_allele_freq'){
								$dbsnp_freq = $line_data[$col_count];
							}						          
							if($self->{score} && $tmp_col eq 'sift_score'){
								$sift_score = $line_data[$col_count];
							}
							if($self->{score} && $tmp_col eq 'polyphen_score'){
								$polyphen_score = $line_data[$col_count];
							}
							if($self->{score} && $tmp_col eq 'cadd_phred'){
								$cadd_phred = $line_data[$col_count];
							}
							if($self->{score} && $tmp_col eq 'go'){
								for my $keyword (@phenotype_keywords){
									if($line_data[$col_count] =~ $keyword || $line_data[$col_count] eq 'NO_GO'){
										$goScore = 1;
									}
								}
							}					
							if($self->{score} && $tmp_col eq 'phenotype'){
								for my $keyword (@phenotype_keywords){
									if($line_data[$col_count] =~ $keyword){
										$phenotypeScore = 2;
									}elsif($line_data[$col_count] eq 'NO_PHENOTYPE'){
										$phenotypeScore = 1;
									}
								}
							}							
						} 
						
						if($self->{score}){
							
							if($sample_str =~ $proband_samplename && $pass_exon_overlap && $pass_status && $pass_filter_exac ){
							
								for my $cutoff (keys %dbsnp_freq_score){
									if($cutoff eq 'N/A' || $cutoff eq 'No_freq'){
										if ($dbsnp_freq eq $cutoff){
											#print "cutoff : $cutoff, dbsnp freq: $dbsnp_freq, score: $dbsnp_freq_score{$cutoff}\n";
											$freqScore += $dbsnp_freq_score{$cutoff};
										}
									} elsif($dbsnp_freq ne 'N/A' && $dbsnp_freq ne 'No_freq'){
										if($dbsnp_freq < $cutoff){
											#print "cutoff : $cutoff, dbsnp freq: $dbsnp_freq, score: $dbsnp_freq_score{$cutoff}\n";
											$freqScore += $dbsnp_freq_score{$cutoff};
										}
									}
								}
								
								if($sift_score eq 'N/A' || $polyphen_score eq 'N/A'){
									$damagingScore += 1;
								}elsif ($sift_score < $sift_score_cutoff || $polyphen_score > $polyphen_score_cutoff){
									$damagingScore += 1;
								}
								if($cadd_phred eq 'N/A'){
									$damagingScore += 1;
								}elsif($cadd_phred > $cadd_score_cutoff){
									$damagingScore += 1;
								}
								
								$final_score = $freqScore + $damagingScore + $goScore + $phenotypeScore;
								#print "$chr:$coord freqS: $freqScore, damagingS: $damagingScore, goS $goScore, phenoS: $phenotypeScore, finalS = $final_score\n";
							}
							
							$mutant_snv_data{"$chr:$coord"}{$end}{$var_base}{'variant_score_total'} = $final_score;
							$mutant_snv_data{"$chr:$coord"}{$end}{$var_base}{'gene_filter'} = $match_gene_filter;
							
							push @line_data, ($final_score,$match_gene_filter);
							$line = join("\t",@line_data);
						}
				
				} else {

						my $pass_exon_overlap = 0;
						my $pass_status = 0;
						my $pass_filter_exac = 0;   
						my $freqScore = 0;
						my $damagingScore = 0;
						my $goScore = 0;
						my $phenotypeScore = 0; 
						my $dbsnp_freq;
#						my $sift_score;
#						my $polyphen_score;
#						my $cadd_phred;
						my $final_score = 0;
						my $match_gene_filter = 0;
											
						for (my $col_count = 0 ; $col_count < @all_indel_headers ; $col_count++ ) {
							my $tmp_col = $all_indel_headers[$col_count];
							$tmp_col =~ s/^%/percent/;
							$tmp_col =~ s/\(/_/g;
							$tmp_col =~ s/\)/_/g;
							$tmp_col =~ s/ //g;
							$tmp_col =~ s/\//_of_/;
							$mutant_indel_data{"$chr:$coord"}{$end}{$var_base}{$tmp_col} = $line_data[$col_count];

							if($self->{score} && $tmp_col eq 'exon_overlap'){
								for my $exon_type (@exon_overlap){
									if($line_data[$col_count] eq $exon_type){
										$pass_exon_overlap = 1;
									}
								}
							}
							if($self->{score} && $tmp_col eq 'gene'){
								for my $gene (@gene_filter){
									if($line_data[$col_count] eq $gene){
										$match_gene_filter = 1;
									}
								}
							}
							if($self->{score} && $tmp_col eq 'final_status'){
								if($line_data[$col_count] eq $final_status){
									$pass_status = 1;
								}
							}
							if($self->{score} && $tmp_col eq 'filter_exac_indel'){
								if($line_data[$col_count] eq 'N/A'){
									$pass_filter_exac = 1;
								}
								else{
									my ($a,$freq) = split(":",$line_data[$col_count]);
									if ($freq < $filter_exac_freq_cutoff){
										$pass_filter_exac = 1;
									}
								}
							}
							if($self->{score} && $tmp_col eq 'dbsnp_var_allele_freq'){
								$dbsnp_freq = $line_data[$col_count];
							}						          
#							if($self->{score} && $tmp_col eq 'sift_score'){
#								$sift_score = $line_data[$col_count];
#							}
#							if($self->{score} && $tmp_col eq 'polyphen_score'){
#								$polyphen_score = $line_data[$col_count];
#							}
#							if($self->{score} && $tmp_col eq 'cadd_phred'){
#								$cadd_phred = $line_data[$col_count];
#							}
							if($self->{score} && $tmp_col eq 'go'){
								for my $keyword (@phenotype_keywords){
									if($line_data[$col_count] =~ $keyword || $line_data[$col_count] eq 'NO_GO'){
										$goScore = 1;
									}
								}
							}					
							if($self->{score} && $tmp_col eq 'phenotype'){
								for my $keyword (@phenotype_keywords){
									if($line_data[$col_count] =~ $keyword){
										$phenotypeScore = 2;
									}elsif($line_data[$col_count] eq 'NO_PHENOTYPE'){
										$phenotypeScore = 1;
									}
								}
							}	
						}
						if($self->{score}){
							
							if($sample_str =~ $proband_samplename && $pass_exon_overlap && $pass_status && $pass_filter_exac ){
							
								for my $cutoff (keys %dbsnp_freq_score){
									if($cutoff eq 'N/A' || $cutoff eq 'No_freq'){
										if ($dbsnp_freq eq $cutoff){
											#print "cutoff : $cutoff, dbsnp freq: $dbsnp_freq, score: $dbsnp_freq_score{$cutoff}\n";
											$freqScore += $dbsnp_freq_score{$cutoff};
										}
									} elsif($dbsnp_freq ne 'N/A' && $dbsnp_freq ne 'No_freq'){
										if($dbsnp_freq < $cutoff){
											#print "cutoff : $cutoff, dbsnp freq: $dbsnp_freq, score: $dbsnp_freq_score{$cutoff}\n";
											$freqScore += $dbsnp_freq_score{$cutoff};
										}
									}
								}
								
#								if($sift_score eq 'N/A' || $polyphen_score eq 'N/A'){
#									$damagingScore += 1;
#								}elsif ($sift_score < $sift_score_cutoff || $polyphen_score > $polyphen_score_cutoff){
#									$damagingScore += 1;
#								}
#								if($cadd_phred eq 'N/A'){
#									$damagingScore += 1;
#								}elsif($cadd_phred > $cadd_score_cutoff){
#									$damagingScore += 1;
#								}
								
								$final_score = $freqScore + $damagingScore + $goScore + $phenotypeScore;
								#print "$chr:$coord freqS: $freqScore, damagingS: $damagingScore, goS $goScore, phenoS: $phenotypeScore, finalS = $final_score\n";
							}
							
							$mutant_indel_data{"$chr:$coord"}{$end}{$var_base}{'variant_score_total'} = $final_score;
							$mutant_indel_data{"$chr:$coord"}{$end}{$var_base}{'gene_filter'} = $match_gene_filter;
							
							push @line_data, ($final_score,$match_gene_filter);
							$line = join("\t",@line_data);
						}						
				}
				
						
				#for sorting
				my $new_chr = $chr;
				
				if ($chr eq 'X') {
					$new_chr = 23;
				} elsif ($chr eq 'Y') {
					$new_chr = 24;
				}
				
				
				
				if ($var_type eq 'snv') {
					$self->{source_info}{$sg_name}{snv_lines}{$pass}{$likely_aff_allele_count}{$likely_unaff_allele_count}{$dbsnp_freq}{$new_chr}{$coord} = $line;
				} else {
					$self->{source_info}{$sg_name}{indel_lines}{$pass}{$likely_aff_allele_count}{$likely_unaff_allele_count}{$dbsnp_freq}{$new_chr}{$coord} = $line;
				}
			}
		}
		#create the xml object
						
		my ($snv_xml,$outdir) = fileparse($self->{source_info}{$sg_name}{snv_xml});
		my ($indel_xml,undef) = fileparse($self->{source_info}{$sg_name}{indel_xml});
		my $var_xml = modules::VariantXML->new(-outdir=>$outdir);
		$var_xml->create_var_xml(-file_name=>$snv_xml,-data=>\%mutant_snv_data,-chr=>'all');
		$var_xml->create_var_xml(-file_name=>$indel_xml,-data=>\%mutant_indel_data,-chr=>'all');
		
	}
}

#Special handling for generating sv reports
sub generate_sv_report {
	my ($self) = @_;
	
	
	
	for my $sg_name (keys %{$self->{source_info}}) {
		my $group_count = 1;
		my $sv_out = $self->{source_info}{$sg_name}{sv_file};
		my $aff_total = $self->{source_info}{$sg_name}{aff_count};
		my $unaff_total = defined $self->{source_info}{$sg_name}{unaff_count}?$self->{source_info}{$sg_name}{unaff_count}:0;
		my $proband = $self->_get_proband_name($sg_name);
		
		my @all_sv_headers = qw(samples priority number_affected_variant/total number_unaffected_variant/total total_variants_in_gene);
		my %sv_headers = %{$self->{source_info}{$sg_name}{sv_headers}};
		for my $sv_header (sort {$sv_headers{$a} <=> $sv_headers{$b}} keys %sv_headers) {
			push @all_sv_headers, $sv_header;
		}
		
		
		open(SV,">$sv_out") || modules::Exception->throw("Can't open file $sv_out\n");
		print SV join("\t",
						'Proband/Sample?',
						'Unique_group_count',
						@all_sv_headers) . "\n\n";
		
		my %ordered_lines = ();
		my @priority = qw(HIGH MEDIUM LOW);
		
		#First report high priority
		for my $sv_type (sort keys %{$self->{source_info}{$sg_name}{var_data}{sv}}) {
			for my $sv_key (sort {$a<=>$b} keys %{$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}}) {
				my $samples = join(",",$proband,
											keys %{$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{samples}});
				
				my $affected_count = () = $samples =~ /_affected/g;
				my $unaffected_count = () = $samples =~ /unaffected/g;
				
				my $gene_name = $self->{source_info}{$sg_name}{gene_lookup}{$sv_key};
				for my $priority ( @priority ) {
					if (exists $self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{$priority}) {
						for my $line_count ( sort {$a<=>$b}  keys %{$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{$priority}} ) {
							my @proband_line = ();
							my $gene_count = keys %{$self->{source_info}{$sg_name}{gene_counts}{$gene_name}{uniq_aff}};
							
							push @proband_line,'proband',$group_count,$samples,$priority,$affected_count.'/'.$aff_total,$unaffected_count.'/'.$unaff_total,$gene_count;
							for my $key_combos (sort {my ($a_index) = split(':',$a); my ($b_index) = split(':',$b); $a_index <=> $b_index} keys %{$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{$priority}{$line_count}}) {
								push @proband_line,$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{$priority}{$line_count}{$key_combos};
							}
							
							$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$sv_key}{proband}{$line_count} = \@proband_line;
						}
						
						for my $sample (keys %{$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{samples}}) {
							
							#print Dumper $self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key};
							
							for my $sample_line (@{$self->{source_info}{$sg_name}{var_data}{sv}{$sv_type}{$sv_key}{samples}{$sample}}) {
								my @sample_line = join("\t",
														'sample',
														$group_count,
														$sample,
														$sample_line
														);
								push @{$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$sv_key}{sample}}, @sample_line; 
							}
						}
						$group_count++;
					}	
				}
			}
		}
		
		my @priorities = qw(HIGH MEDIUM LOW);
		
		#print Dumper \%ordered_lines;
		
		
		for my $affected_count (sort {$b<=>$a} keys %ordered_lines) {
			for my $unaffected_count (sort {$a<=>$b} keys %{$ordered_lines{$affected_count}}) {
				for my $priority (@priorities) {
					
					for my $sv_type (keys %{$ordered_lines{$affected_count}{$unaffected_count}{$priority}}) {
						for my $event_number (sort {$a<=>$b} keys %{$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}}) {
							for my $line_count (sort {$a<=> $b} keys %{$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$event_number}{proband}}) {
								if ($line_count > 0) {
									my @entries = @{$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$event_number}{proband}{$line_count}};
									$entries[3] = '';
									$entries[4] = '';
									$entries[5] = '';
									$entries[6] = '';
#									print SV join("\t",
#													$entries[0],
#													$entries[1],
#													'',
#													'',
#													'',
#													'',
#													'',
#													@entries[6..-1]
#													) ."\n";
									print SV join("\t",@entries) ."\n";
								} else {
									print SV join("\t",
												@{$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$event_number}{proband}{$line_count}}
												) ."\n";
								}
								
							}
							print SV join("\n",
										@{$ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$event_number}{sample}}
										) ."\n" if exists $ordered_lines{$affected_count}{$unaffected_count}{$priority}{$sv_type}{$event_number}{sample};
						}
					}
					print SV "\n";
				}
				print SV "\n";
			}
			print SV "\n";
		}
		
		
		
		
	}
}

	
#write to files	
sub write_to_files {
	my ($self) = @_;
		
	for my $sg_name (keys %{$self->{source_info}}) {
		my %unsorted_snv_lines = %{$self->{source_info}{$sg_name}{snv_lines}} if exists $self->{source_info}{$sg_name}{snv_lines};
		my %unsorted_indel_lines = %{$self->{source_info}{$sg_name}{indel_lines}} if exists $self->{source_info}{$sg_name}{indel_lines}; #May not exist for small targeted projects
		my $snv_out = $self->{source_info}{$sg_name}{snv_file};
		my $indel_out = $self->{source_info}{$sg_name}{indel_file};
		
#		if (-e $snv_out) {
#			#If it's a project we append to file
#			open(SNV,">>$snv_out") || modules::Exception->throw("Can't open file to write $snv_out\n");
#		} else {
			open(SNV,">$snv_out") || modules::Exception->throw("Can't open file to write $snv_out\n");
			print SNV $self->{source_info}{$sg_name}{snv_headers} . "\n\n";
#		}
		
		#Sort by non-mendelian first; then highest likely affected variant and lowest unaffected non-variant 
		for my $pass (sort {$b<=>$a} keys %unsorted_snv_lines) {
			if ($pass == 1) {
				print SNV "\nNOVEL_OR_RARE SNVS\n\n";
			} else {
				print SNV "\n\nLOW_PRIORITY SNVS\n\n";
			}
			for my $aff_count (sort {$b<=>$a} keys %{$unsorted_snv_lines{$pass}}) {
				for my $unaff_count (sort {$a<=>$b} keys %{$unsorted_snv_lines{$pass}{$aff_count}}) {
					for my $dbsnp_freq (sort {$a<=>$b} keys %{$unsorted_snv_lines{$pass}{$aff_count}{$unaff_count}}) {
						for my $chr (sort {$a <=> $b} keys %{$unsorted_snv_lines{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}}) {
							for my $coord (sort {$a<=>$b} keys %{$unsorted_snv_lines{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}{$chr}}) {
								#print "Mendel $mendelian Pass $pass Aff $aff_count Unaff $unaff_count Chr $chr Coord $coord\n";
								print SNV $unsorted_snv_lines{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}{$chr}{$coord}. "\n";
							}
						}
					}
				}
			}
		}
		
		close SNV;
	
#		if (-e $indel_out) {
#			open(INDEL,">>$indel_out") || modules::Exception->throw("Can't open file to write $indel_out\n");
#		} else {
			open(INDEL,">$indel_out") || modules::Exception->throw("Can't open file to write $indel_out\n");
			print INDEL $self->{source_info}{$sg_name}{indel_headers} . "\n\n";
#		}
		
		for my $pass (sort {$b<=>$a} keys %unsorted_indel_lines) {
			if ($pass == 1) {
				print INDEL "\nNOVEL_OR_RARE INDELS\n\n";
			} else {
				print INDEL "\n\nLOW_PRIORITY INDELS\n\n";
			}
			for my $aff_count (sort {$b<=>$a} keys %{$unsorted_indel_lines{$pass}}) {
				for my $unaff_count (sort {$a<=>$b} keys %{$unsorted_indel_lines{$pass}{$aff_count}}) {
					for my $dbsnp_freq (sort {$a<=>$b} keys %{$unsorted_indel_lines{$pass}{$aff_count}{$unaff_count}}) {
						for my $chr (sort {$a <=> $b} keys %{$unsorted_indel_lines{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}}) {
							for my $coord (sort {$a<=>$b} keys %{$unsorted_indel_lines{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}{$chr}}) {
								#print "Mendel $mendelian Pass $pass Aff $aff_count Unaff $unaff_count Chr $chr Coord $coord\n";
								print INDEL $unsorted_indel_lines{$pass}{$aff_count}{$unaff_count}{$dbsnp_freq}{$chr}{$coord}. "\n";
							}
						}
					}
				}
			}
		}
		
		close INDEL;
	 }
}	
	

#Create db entries	
sub write_to_db {
	my ($self) = @_;
	
	for my $sg_name (keys %{$self->{source_info}}) {
		my $out;
		
		if ($self->{sv}) {
			$out = $self->{source_info}{$sg_name}{sv_file};
		} else {
			$out = $self->{source_info}{$sg_name}{snv_file};
		}
	
		
			
		my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime time;
		my $timestr = sprintf( "%02d-%02d-%d_%02d:%02d",
								   $mday, 
								   $mon + 1, 
								   $year + 1900,
								   $hour, 
								   $min );
		
		my $sqltime = sprintf( "%d-%02d-%02d",$year+1900,$mon+1,$mday);

		if ($self->{report_type} eq 'project') {
			
			my $project_name = $self->{project_name};
			my ($project_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
			my $sample_total = $self->{project_info}{sample_number};
			my %project_summary = (
								total_samples=>$sample_total,
								file_name=>$out,
								project_id=>$project_obj->id,
								summary_date => $sqltime
								);
			my $summary_id = modules::Adaptors::Project_Summary->insert(\%project_summary);
			print STDERR "Create summary db entry $summary_id\n";
			return; #Just create one project entry		
		} else {
			my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$sg_name);
			#Insert the db object to prevent from running again
			my $sample_total = $self->{source_info}{$sg_name}{sample_number}; #= $snv_out =~ /(\d+)_sample/;
			my %group_summary = (
								total_samples=>$sample_total,
								file_name=>$out,
								source_group_id=>$source_group_obj->id,
								summary_date => $sqltime
								);
			my $summary_id = modules::Adaptors::Group_Summary->insert(\%group_summary);
			print STDERR "Create summary db entry $summary_id\n";
		}
		
		
		
	}
}	
	
#Write files to mdss for source_group cases
sub mdss {
	my ($self) = @_;
	
	my $sys_call = modules::SystemCall->new();
	my $clus_conf = modules::Pipeline::get_cluster_conf();	
	
	for my $sg_name (keys %{$self->{source_info}}) {
		
		my $source_type = modules::Pipeline::get_source_type(-source_group_name=>$sg_name);
		my $source_name = modules::Pipeline::get_source_name($sg_name);
		my $mdss_results_dir = $clus_conf->read('common','mdss','mdss_results');
		my $mdss_mkdir_cmd = "mdss mkdir $mdss_results_dir/$source_type/$source_name/source_group_summaries";
		$sys_call->run($mdss_mkdir_cmd);
		my $sg_dir = $self->{source_info}{$sg_name}{dir_base};
		my $mdss_command = "mdss put $sg_dir/* $mdss_results_dir/$source_type/$source_name/source_group_summaries";
		$sys_call->run($mdss_command);  
		
	}
}	
	
sub _get_proband_name {
	my ($self,$sg_name) = @_;
	my $proband;
	
	my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$sg_name);
			
	my $source_type = modules::Pipeline::get_source_type(-source_group_name=>$sg_name);
	my @samples = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);	
	
	if ($source_type =~ /human_related/) {
		#Get proband here
		for my $sample_obj ( @samples ) {
			my $sample_name = $sample_obj->sample_name;
			
			if ($self->{source_info}{$sg_name}{samples}{$sample_name}{relation} eq 'proband'){
				return $sample_name;
			}
			
		}
				
	} else {
		#Arbitrarily pick first sample
		for my $sample_obj ( @samples ) {
			my $sample_name = $sample_obj->sample_name;
			return $sample_name;
		}	
	}
	modules::Exception->throw("ERROR: Can't find proband for sg_name $sg_name");	
	
}


sub _check_block {
	my ($self,$mother_block_start,$mother_block_end,$father_block_start,$father_block_end,$min_variant_num,$current_var_count,$min_block) = @_;
	
	if ($current_var_count >= $min_variant_num) {
		#print "if ($mother_block_start != 0 && ($mother_block_end-$mother_block_start) > $min_block\n";
		if ($mother_block_start != 0 && ($mother_block_end-$mother_block_start) > $min_block) {
			return 1;
		}
		if ($father_block_start != 0 && ($father_block_end-$father_block_start) > $min_block) {
			return 1;
		}
	} 
		
	return 0;
	
	
	
}	


#Get alleles from pileup string and zygosity
sub _get_alleles {
	my ($self,$zyg,$pileup_str) = @_;
	if ($zyg eq 'No data') {
		return ('?','?');
	} elsif ($zyg eq 'ref') {
		return ('ref','ref');
	} elsif ($zyg eq 'hom') {
		my ($var_allele) = split(':',$pileup_str);
		return ($var_allele,$var_allele);
	} else {
		#het case
		my ($block1,$block2) = split('_',$pileup_str);
		
		if (!defined $block2) {
			modules::Exception->warning("ERROR: Het zyg and pileup str $pileup_str");
			return('?','?');
		}
		
		my ($first_allele) = split(':',$block1);
		my ($second_allele) = split(':',$block2);
		#Return ref as first allele always
		if ($block2 =~ /ref/i){
			return ($second_allele,$first_allele)
		} else {
			return ($first_allele,$second_allele);
		}
	}
}

#Checks whether Mendelian inheritance is upheld with one or more parent sequenced
sub _mendel_inheritance {
	my ($self,$fam_data,$sg_name) = @_;
	
	
	if ($self->_has_two_parents($sg_name)) {
		#Sequence data from both parents
		my $mother_zyg = my $father_zyg;
		for my $family_key (keys %{$fam_data}) {
			my ($sample,$relation) = split(':',$family_key);
			if ($relation eq 'mother') {
				$mother_zyg = $fam_data->{$family_key}{zyg};
				if ($mother_zyg eq 'No data') {
					return 'missing allele info'; #Inconclusive at best so return 0
				}	
			} elsif ($relation eq 'father') {
				$father_zyg = $fam_data->{$family_key}{zyg};
				if ($father_zyg eq 'No data') {
					return 'missing allele info'; #Inconclusive at best so return 0
				}		
			} else {
				if ($fam_data->{$family_key}{zyg} eq 'No data') {
					return 'missing allele info'; #Inconclusive at best so return 0
				}
			}
		}
		
		
		for my $family_key (keys %{$fam_data}) {
			my ($sample,$relation) = split(':',$family_key);
			next if $relation ne 'proband';
			#next if $relation eq 'mother' || $relation eq 'father';
			my $sample_zyg = $fam_data->{$family_key}{zyg};
			
			if ($sample_zyg eq 'ref') {
				if ($mother_zyg eq 'hom') {
					return "hom_mother ref_child";
				} elsif ($father_zyg eq 'hom') {
					return "hom_father ref_child";
				}
			} elsif ($sample_zyg eq 'het') {
				if ($mother_zyg eq 'ref' && $father_zyg eq 'ref') { 
					return "de novo (ref_parents het_child)";
				} 
				if ($mother_zyg eq 'hom' && $father_zyg eq 'hom') {
					return "hom_parents het_child";
				}
				
			} elsif ($sample_zyg eq 'hom') {
				if ($mother_zyg eq 'ref') {
					return "ref_mother hom_child";
				}
				if ($father_zyg eq 'ref') {
					return "ref_father hom_child";
				}
			} elsif ($sample_zyg eq 'No data') {
				return 'missing allele info'; #Inconclusive at best so return 0
			} else {
				modules::Exception->throw("ERROR: No zyg for $family_key");
			}	
				
		}
	} elsif ($self->_has_single_parent($sg_name)) {
		#Sequence from only one parent
		my $parent_zyg;
		for my $family_key (keys %{$fam_data}) {
			my ($sample,$relation) = split(':',$family_key);
			if ($relation eq 'mother') {
				$parent_zyg = $fam_data->{$family_key}{zyg};
				if ($parent_zyg eq 'No data') {
					return 'missing allele info'; #Inconclusive at best so return 0
				}	
			} elsif ($relation eq 'father') {
				$parent_zyg = $fam_data->{$family_key}{zyg};
				if ($parent_zyg eq 'No data') {
					return 'missing allele info'; #Inconclusive at best so return 0
				}
			}
		}
		
		for my $family_key (keys %{$fam_data}) {
			my ($sample,$relation) = split(':',$family_key);
			next if $relation ne 'proband';
			#next if $relation eq 'mother' || $relation eq 'father';
			my $sample_zyg = $fam_data->{$family_key}{zyg};
			
			if ($sample_zyg eq 'ref') {
				if ($parent_zyg eq 'hom') {
					return 'hom_parent ref_child';
				}
			} elsif ($sample_zyg eq 'hom') {
				if ($parent_zyg eq 'ref') {
					return 'ref_parent hom_child';
				}
			} 	
				
		}
		
	}  else {
		modules::Exception->throw("ERROR: No parental data available");
	}
	return 'yes';
}						

#Get the disease inheritance patterns for dominant or recessive
sub _disease_inheritance {
	my ($self,$fam_data,$pileup_key,$sg_name) = @_;
	
	my ($chr) = split(":",$pileup_key);
	
	my $x_or_auto = $chr eq 'X'?'X':'auto';
	
	my $mother_var = 0;
	my $father_var = 0;
	my $mother_aff = 0;
	my $father_aff = 0;
	
	my $aff_count = 0;
	my $unaff_count = 0;
	my $recessive = 1;
	my $dominant = 1;
	my $disease_inheritance = 'none';
	

	#Use the sample info from the db as just b/c variant isn't called a snv doesn't mean it isn't variant (e.g. just below 40 cutoff score)
	for my $sample (keys %{$self->{source_info}{$sg_name}{samples}}) {

		if (!exists $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg}) {
			modules::Exception->warning("ERROR: No data for $sample for $pileup_key\n");
			$self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{pileup} = 'No data';
			$self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg} = 'No data';
		}

		#Skip the entries that are ref or no data
		if ($self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg} eq 'ref' || $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg} eq 'No data') {
			next;
		}
		
		
		my $relation = $self->{source_info}{$sg_name}{samples}{$sample}{relation};
		my $zyg = $self->{source_info}{$sg_name}{pileup_str}{$pileup_key}{$sample}{zyg};
		
		#Need to find mother and father to show it's auto-recessive
		if ($relation eq 'mother') {
			$mother_var = 1;
			if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
				$mother_aff = 1;
			}
		}
		if ($relation eq 'father') {
			$father_var = 1;
			if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
				$father_aff = 1;
			}
		}
		
		
		if ($self->{source_info}{$sg_name}{samples}{$sample}{affected} == 1) {
			#Affected samples
			$aff_count++;
			if ($zyg ne 'hom') {
				#Affected must be homs for auto-recessive
				$recessive = 0;
			}
			if ($zyg ne 'het') {
				#Affected must be het for auto-dominant
				$dominant = 0;
			}
		} else {
			#unaffected samples
			$unaff_count++;
			if ($zyg eq 'hom') {
				#Unaffected can't be homs for auto-recessive
				$recessive = 0;
			}
			#Unaffected can't be variant for auto-dominant
			$dominant = 0;
		}
	}
		
	if ($self->_has_two_parents($sg_name)) {
		#Sequence data from both parents
		if ($self->{source_info}{$sg_name}{aff_count} == $aff_count && $recessive) {
			#If all affected are hom mutant and unaffected either WT or het
			if ($mother_var && $father_var) {
				#Here it's recessive
				$disease_inheritance = $x_or_auto . '-recessive';
			} 
		
		} elsif ($self->{source_info}{$sg_name}{aff_count} == $aff_count && $dominant) {
			#If all affecteds are mutant (het or hom) and all unaffected WT; also check one parent mutant and mutant parent(s) affected
			if ($mother_var && $father_var) {
				$disease_inheritance = $x_or_auto . '-dominant' if ($mother_aff && $father_aff);
			} elsif ($mother_var) {
				$disease_inheritance = $x_or_auto . '-dominant' if $mother_aff;
			} elsif ($father_var) {
				$disease_inheritance = $x_or_auto . '-dominant' if $father_aff;
			} else {
				$disease_inheritance = $x_or_auto . '-dominant';
			}
		} 
		
	} elsif ($self->_has_single_parent($sg_name)) {
		if ($self->{source_info}{$sg_name}{aff_count} == $aff_count && $recessive) {
			if ($mother_var || $father_var) {
				#Here it's recessive
				$disease_inheritance = $x_or_auto . '-recessive';
			} 
		} elsif ($self->{source_info}{$sg_name}{aff_count} == $aff_count && $dominant) {
			#If all affecteds are mutant (het or hom) and all unaffected WT; also check one parent mutant and mutant parent(s) affected
			if ($mother_var) {
				$disease_inheritance = $x_or_auto . '-dominant' if $mother_aff;
			} elsif ($father_var) {
				$disease_inheritance = $x_or_auto . '-dominant' if $father_aff;
			} else {
				$disease_inheritance = $x_or_auto . '-dominant';
			}
		} 
	} else {
		#here we can't say anything about parents but can still see patterns
		if ($self->{source_info}{$sg_name}{aff_count} == $aff_count && $recessive) {
			$disease_inheritance = $x_or_auto . '-recessive';	
		} elsif ($self->{source_info}{$sg_name}{aff_count} == $aff_count && $dominant) {
			$disease_inheritance = $x_or_auto . '-dominant';
		} 
	}
	return $disease_inheritance;

	
}

sub _has_mother {
	my ($self,$sg_name) = @_;
	if ($self->{source_info}{$sg_name}{mother}) {
		return 1;
	} else {
		return 0;
	}	
}

sub _has_father {
	my ($self,$sg_name) = @_;
	if ($self->{source_info}{$sg_name}{father}) {
		return 1;
	} else {
		return 0;
	}
}

sub _has_two_parents {
	my ($self,$sg_name) = @_;
	if ($self->{source_info}{$sg_name}{father} && $self->{source_info}{$sg_name}{mother}) {
		return 1;
	} else {
		return 0;
	}
}

sub _has_single_parent {
	my ($self,$sg_name) = @_;
	if ($self->{source_info}{$sg_name}{father} || $self->{source_info}{$sg_name}{mother}) {
		return 1;
	} else {
		return 0;
	}
}

return 1;
