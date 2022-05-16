package modules::PipelineSample;

#module contains pipeline functions specific to samples

use strict;
use modules::QualityEncoding;
use modules::Exception;
use Data::Dumper;
use modules::Cluster;
use modules::ConfigXML;
use modules::Adaptors::Syscall;
use modules::Adaptors::Source;
use modules::Adaptors::Sample;

sub new
{
	my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

	#xml objects that contain the config data we need to access
    my @required_args = (
			         	'-source_name',
			         	'-cluster_obj',
			         	'-source_type'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    #Set the global variables
    $self->{source_name} = $args{'-source_name'};
    $self->{source_type} = $args{'-source_type'};
    $self->{cluster_obj} = $args{-cluster_obj};
    
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
	$self->_pipe_xml($pipe_config);
	
	my $cluster_xml = "$svndir/conf/cluster.xml";

	if ( !-e $cluster_xml ) {
		modules::Exception->throw("File $cluster_xml doesn't exist");	
	}
	
	my $cluster_config = modules::ConfigXML->new($cluster_xml);
	$self->_cluster_xml($cluster_config);    
	
    return $self;
}

#Set steps_xml
sub _pipe_xml
{
    my ($self, $pipe_xml) = @_;

    if (defined $pipe_xml) {
		$self->{'pipe_xml'} = $pipe_xml;
    } elsif (! defined $self->{'pipe_xml'}) {
		modules::Exception->throw("pipe_xml not set");
    }

    return $self->{'pipe_xml'};
}

#Set steps_xml
sub _cluster_xml
{
    my ($self, $cluster_xml) = @_;

    if (defined $cluster_xml) {
		$self->{'cluster_xml'} = $cluster_xml;
    } elsif (! defined $self->{'cluster_xml'}) {
		modules::Exception->throw("cluster_xml not set");
    }

    return $self->{'cluster_xml'};
}

sub _svndir
{
	my ($self) = @_;
	return $self->_cluster_xml->read('common','svndir');
}


sub _threadnum
{
	my ($self) = @_;
	return $self->_cluster_xml->read('common','qsub_vars','thread_num');
}

sub _samtools
{
	my ($self) = @_;
	return $self->_pipe_xml->read($self->{source_type},'binaries','samtools','binary');
}

sub _bwa
{
	my ($self) = @_;
	return $self->_pipe_xml->read($self->{source_type},'binaries','bwa','binary');
}

sub _qsub_dir
{
	my ($self) = @_;
	return $self->_cluster_xml->read($self->{source_type},'base_directories','base_qsub_directory');
}
sub _sym_bam_path # value set in create_sample_xml
{
	my ($self) = @_;

	modules::Exception->throw("Undefined _sym_bam_path value") 
	    unless defined $self->{_sym_bam_path};

	return $self->{_sym_bam_path};
}


#Generic parser for converting steps.xml template file into xml used by pipeline
sub _parse_xml 
{
	my ($self,$file,$source_type,$args) = @_;
	my $xml_value;
	if ($file eq 'cluster') {
		if ($self->_cluster_xml->exists(@{$args})) {
			$xml_value = $self->_cluster_xml->read(@{$args});		
		} else {
			$xml_value = $self->_cluster_xml->read($source_type,@{$args});
		}
	} elsif ($file eq 'pipe') {
		if ($self->_pipe_xml->exists(@{$args})) {
			$xml_value = $self->_pipe_xml->read(@{$args});		
		} else {
			$xml_value = $self->_pipe_xml->read($source_type,@{$args});
		}
	} else {
		modules::Exception->throw("ERROR: File argument must be pipe or cluster\n");
	}
	return $xml_value;
}

#Checks if a step exists in the config file when -resume flag used
sub _step_exists {
	my ($self,$check_step,$steps) = @_;
	
	for my $step ( @{$steps} ) {
	    if ($check_step eq $step) {
	    	return 1;
	    }
	}
	return 0;
}


sub create_sample_xml 
{
	my $self = shift;
	
	if (!defined $self->{source_name}) {
		modules::Exception->throw("ERROR: Must call create_sample_xml via pipeline object");
	}
	
    my %args = @_;
	
    my @required_args = (
			             '-sample_name',
			             '-steps_xml_template',
			             '-xml_out',
			             '-encoding',
			             '-bams'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $sample_name = $args{-sample_name};
    my $steps_xml_template = $args{-steps_xml_template};
    my $xml_out = $args{-xml_out};
    my $encoding = $args{-encoding};
    my $source_name = $self->{source_name};
    my $source_group_name = modules::Pipeline::get_source_group_name($sample_name);
    my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_name);
    my $bam_files = join(",",@{$args{-bams}});
	open(TEMPLATE, $steps_xml_template) || modules::Exception->throw("ERROR: Can't open file $steps_xml_template");
	open(XMLOUT, ">$xml_out")  || modules::Exception->throw("ERROR: Can't open file for writing $xml_out");   
	
	while (<TEMPLATE>) {
		if (/COMMANDS/) {
			#Here we write out the sample specific variables
			print XMLOUT "<!ENTITY nameSample \"$sample_name\">\n";
			print XMLOUT "<!ENTITY nameSourceGroup \"$source_group_name\">\n";
			print XMLOUT "<!ENTITY nameSource \"$source_name\">\n";
			print XMLOUT "<!ENTITY encoding \"$encoding\">\n\n";
			print XMLOUT $_;
		} elsif (/(.*ENTITY.*)\"(.*)->\((.*)\)/) {
			my $start = $1;
			my $file = $2;
			my @xml_args = split(" ",$3);
			my $value = $self->_parse_xml($file,$self->{source_type},\@xml_args);
			print XMLOUT "$start \"$value\">\n";
		} elsif (/BAMFILES/) {
			#Need to get the input bam files from the lanes into the xml
			$_ =~ s/BAMFILES/$bam_files/g;
			print XMLOUT $_;
		} elsif (/SYM_BAM/) {
			my $bamlinks_dir = "$source_name/bam_links/";
			#Need to have standard names for bam files for snv_calling later
			$_ =~ s/SYM_BAM/$bamlinks_dir$sample_name.bam/g;
			$self->{_sym_bam_path} = "$bamlinks_dir$sample_name.bam";
			print XMLOUT $_;
		} elsif (/TARGETED/) {
#			if ($sequence_type eq 'targeted') {
#				$_ =~ s/TARGETED/-S/;
#			} else {
#				$_ =~ s/TARGETED//;
#			}
			$_ =~ s/TARGETED//; 
			print XMLOUT $_;
		} else {
			print XMLOUT $_;
		}
	}
	close XMLOUT;
}




#Create the snv_calling qsubs
sub create_parallel_cancer_snvcall_qsubs
{
	my $self = shift;
    my %args = @_;
	
    my @required_args = (
			             '-normal_sample',
			             '-tumour_samples',
			             '-source_group_name',
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }	
    my $normal_sample = $args{-normal_sample};
	my @tumour_samples = @{$args{-tumour_samples}};
	
	my $sg_name = $args{-source_group_name};
	my $source_name = $self->{source_name};

       #These are from the xml files
    my $svndir = $self->_svndir();
 	 	
    my $qsub_dir = $self->_qsub_dir().'/'.$source_name.'/snvcalls';
    
   	if ( !-d $qsub_dir ) {
   		system("mkdir -p $qsub_dir");	
   	}
   	
   	my $samtools_bin = $self->_samtools();
   	my $bcftools_bin = $self->_pipe_xml->read($self->{source_type},'binaries','bcftools','binary');
   	
   	my $rundir = $self->_cluster_xml->read($self->{source_type},'base_directories','base_run_directory');
   	my $fasta_ref = $self->_cluster_xml->read($self->{source_type},'svn','fasta_file');
   	my @chrs = split(" ",$self->_pipe_xml->read($self->{source_type},'annotation_version','chr'));
   	my $mpileup_args = $self->_pipe_xml->read($self->{source_type},'binaries','samtools','mpileup','args');
   	my $bcftools_tumour_args = $self->_pipe_xml->read($self->{source_type},'binaries','bcftools','view','paired_args');
   	my $bcftools_normal_args = $self->_pipe_xml->read($self->{source_type},'binaries','bcftools','view','unpaired_args');
 	my $record_snv_bin = "$svndir/scripts/record_snv_call.pl";  	
 
 	#DON'T DO THIS ANYMORE AS HAPPENS WHEN SNVS CALLED SIMULTANEOUSLY WITH TUMOUR SAMPLES
 	if ($args{-call_normal}) {
	   	#Create the normal sample files for normal snv calls
	   	for my $chr ( @chrs ) {
	   			my $vcf_file = "$rundir/$source_name/$sg_name"."_snvcalls/$normal_sample.$chr.vcf";
	   			my $snv_call_qsub_tmp = $normal_sample.".snvcall.$chr.qsub.tmp"; #tmp b/c we need to add the run id later
	   			my $snv_call_qsub = $qsub_dir."/".$normal_sample.".snvcall.$chr.qsub"; #pass this to file
	   			#Mpileup command to create vcf
	   			my $pileup_command = "$samtools_bin $mpileup_args $fasta_ref -r $chr $rundir/$source_name/bam_links/$normal_sample.bam | $bcftools_bin $bcftools_normal_args - | grep -v NNNNN  > $vcf_file";
	   			#Create a record of vcf being created; need to allow mpileup jobs to finish before proceeding to next pipeline steps
	   			my $record_command = "$record_snv_bin -vcf $vcf_file -chr $chr -step_name single_vcf -qsub_file $snv_call_qsub -runid RUNID"; #RUNID gets replaced by submit_snvs.pl
	   			my @commands = ($pileup_command,$record_command);
	    		$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$snv_call_qsub_tmp,-commands=>\@commands,-cpus=>1, -mem=>'4Gb',-walltime=>'48:00:00');
	   	}
 	}
 
 	#Now create all the tumour sample files (paired with normal) for each tumour snv call
   	for my $tumour_sample ( @tumour_samples ) {
   		for my $chr ( @chrs ) {
   			my $vcf_file = "$rundir/$source_name/".$sg_name."_snvcalls/$tumour_sample.$chr.vcf";
   			my $snv_call_qsub_tmp = $tumour_sample.".snvcall.$chr.qsub.tmp"; #tmp b/c we need to add the run id later
   			my $snv_call_qsub = $qsub_dir."/".$tumour_sample.".snvcall.$chr.qsub"; #pass this to file
   			#Mpileup command to create vcf
   			my $pileup_command = "$samtools_bin $mpileup_args $fasta_ref -r $chr $rundir/$source_name/bam_links/$normal_sample.bam $rundir/$source_name/bam_links/$tumour_sample.bam | $bcftools_bin $bcftools_tumour_args - | grep -v NNNNN  > $vcf_file";
   			#Create a record of vcf being created; need to allow mpileup jobs to finish before proceeding to next pipeline steps
   			my $record_command = "$record_snv_bin -vcf $vcf_file -chr $chr -step_name single_vcf -qsub_file $snv_call_qsub -runid RUNID"; #RUNID gets replaced by submit_snvs.pl
   			my @commands = ($pileup_command,$record_command);
    		$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$snv_call_qsub_tmp,-commands=>\@commands,-cpus=>1, -mem=>'4Gb',-walltime=>'48:00:00');
   		}
   	}
}

#Create the snv_calling qsubs
sub create_parallel_snvcall_qsubs
{
    my $self = shift;
    my %args = @_;
	
    my @required_args = ('-sample',
			 '-source_group_name',
			 '-steps_xml'
	);

    foreach my $required_arg (@required_args){

	if (! defined $args{$required_arg}){
	    modules::Exception->throw("Required argument [$required_arg] not set");
	}
    }	
    my $phred64 = defined $args{-phred64}?1:0;
    my $sample = $args{-sample};	
    my $sg_name = $args{-source_group_name};
    my $source_name = $self->{source_name} if defined $self->{source_name};

    # Retrieve config from configuration xml files
    my $svndir = $self->_svndir();
 	 	
    my $qsub_dir = $self->_qsub_dir().'/'.$source_name.'/snvcalls';
    
    
    if ( !-d $qsub_dir ) {
		system("mkdir -p $qsub_dir");	
    }
   	
    my $samtools_bin = $self->_samtools();
    my $bcftools_bin = $self->_pipe_xml->read($self->{source_type},'binaries','bcftools','binary');
   	
    my $rundir = $self->_cluster_xml->read($self->{source_type},'base_directories','base_run_directory');
    my $fasta_ref = $self->_cluster_xml->read($self->{source_type},'svn','fasta_file');
    my @chrs = split(" ",$self->_pipe_xml->read($self->{source_type},'annotation_version','chr'));
    my $mpileup_args = $self->_pipe_xml->read($self->{source_type},'binaries','samtools','mpileup','args');
    my $bcftools_args = $self->_pipe_xml->read($self->{source_type},'binaries','bcftools','view','unpaired_args');
    my $record_snv_bin = "$svndir/scripts/record_snv_call.pl";  	
 
    ### HERE - variant call commands SET HERE  ###

    if ($self->{source_type} !~ /gatk/i) { # Default workflow for calling SNVs - for now, SAMtools commands

		for my $chr ( @chrs ) {
		    my $vcf_file = "$rundir/$source_name/$sg_name"."_snvcalls/$sample.$chr.vcf";
		    my $snv_call_qsub_tmp = $sample.".snvcall.$chr.qsub.tmp"; #tmp b/c we need to add the run id later
		    my $snv_call_qsub = $qsub_dir."/".$sample.".snvcall.$chr.qsub"; #pass this to file
		    #Mpileup command to create vcf
		    my $pileup_command = "$samtools_bin $mpileup_args $fasta_ref -r $chr $rundir/$source_name/bam_links/$sample.bam | $bcftools_bin $bcftools_args - | grep -v NNNNN  > $vcf_file";
		    #Create a record of vcf being created; need to allow mpileup jobs to finish before proceeding to next pipeline steps
		    my $record_command = "$record_snv_bin -vcf $vcf_file -chr $chr -step_name single_vcf -qsub_file $snv_call_qsub -runid RUNID"; #RUNID gets replaced by submit_snvs.pl
		    my @commands = ($pileup_command,$record_command);
		    $self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$snv_call_qsub_tmp,-commands=>\@commands,-cpus=>1, -mem=>'4Gb',-walltime=>'48:00:00');
		}


    } else { # GATK workflow for SNV calling

		# Some GATK-specific config
		my $java_bin = $self->_pipe_xml->read($self->{source_type},'binaries','javabin','binary');
		$java_bin.=' -Xms3200m -Xmx3600m'; #restrict java heap to 3.6G to fit job within 4G mem
		my $gatkjar = $self->_pipe_xml->read($self->{source_type},'binaries','gatk','jar');
		my $num_threads = $self->_threadnum();
		my $genome_fasta = $fasta_ref;
		#my $filenamestub = $rundir.'/'.$source_name.'/'.$sg_name;
		
		#Get all the GATK arguments from pipe.xml
		my $gatkdatafilesdir = $self->_cluster_xml->read($self->{source_type},'svn','gatk_resources');
		my $intervals_dir =  $self->_cluster_xml->read($self->{source_type},'svn','chr_intervals_dir');
		my $realign_target_args =  $self->_pipe_xml->read($self->{source_type},'binaries','gatk','args_realign_target');
		my $realign_indel_args = $self->_pipe_xml->read($self->{source_type},'binaries','gatk','args_realign_indel');
		my $base_recal_args = $self->_pipe_xml->read($self->{source_type},'binaries','gatk','args_base_recal');
		my $print_reads_args = $self->_pipe_xml->read($self->{source_type},'binaries','gatk','args_print_reads');
		my $haplo_caller_args = $self->_pipe_xml->read($self->{source_type},'binaries','gatk','args_haplo_call');
		my $genotype_gvcfs_args = $self->_pipe_xml->read($self->{source_type},'binaries','gatk','args_genotype_gvcf');

		my $genome_indel_file = $gatkdatafilesdir .'/'.$self->_pipe_xml->read($self->{source_type},'binaries','gatk','genome_indel_file');
		my $mills_indel_file = $gatkdatafilesdir .'/'.$self->_pipe_xml->read($self->{source_type},'binaries','gatk','mills_indel_file');
		my $dbsnp_file = $gatkdatafilesdir .'/'.$self->_pipe_xml->read($self->{source_type},'binaries','gatk','dbsnp_file');
		my $full_bam_file = "$rundir/$source_name/bam_links/$sample.bam";
	
		#To minimise number of total jobs we group into 10 blocks each taking 24 hours to run
		my %chr_groups = (
						1=>[1,21],
						2=>[2,18],
						3=>[3,14],
						4=>[4,15],
						5=>[5,13,'Y'],
						6=>[6,9],
						7=>[7,10],
						8=>[8,11],
						9=>[12,'X',22],
						10=>[16,17,19,20]
						);
		
		 my $steps_xml = $args{-steps_xml};
    	if ( !-e $steps_xml ) {
			modules::Exception->throw("File $steps_xml doesn't exist");	
    	}
    	my $steps_config = modules::ConfigXML->new($steps_xml);
		my $module = $steps_config->read('steps','step','single_vcf','module');
	
		for my $chr_count (sort {$a<=>$b} keys %chr_groups) {
		    # some config for where output files will go
		    my $snv_call_qsub_tmp = $sample.".snvcall.$chr_count.qsub.tmp"; #tmp cos we need to add the run id later
		    my $snv_call_qsub = $qsub_dir."/".$sample.".snvcall.$chr_count.qsub"; #pass this to file
	
	
		    ### Alert! - All GATK commands configured here
			my @gatk_commands = ();
			
			for my $chr (@{$chr_groups{$chr_count}}) {
			    my $vcf_file = "$rundir/$source_name/$sg_name"."_snvcalls/$sample.$chr.vcf";
			    my $intervals_file = "$intervals_dir/".$chr.".intervals";
			    
			    my $realign_target_cmd = "$java_bin -jar $gatkjar $realign_target_args -R $genome_fasta -I $full_bam_file --intervals $intervals_file -known $genome_indel_file -o GATKOUT.realigner.$chr.intervals";
				$realign_target_cmd .= ' --fix_misencoded_quality_scores' if $phred64;
			   	push @gatk_commands, $realign_target_cmd;
			   	my $realign_indel_cmd = "$java_bin -jar $gatkjar $realign_indel_args -R $genome_fasta -I $full_bam_file -known $genome_indel_file --intervals $intervals_file -targetIntervals GATKOUT.realigner.$chr.intervals -o GATKOUT.sorted_dupmarked_readgrp_realigned.".$chr.".bam ";
			   	$realign_indel_cmd .= ' --fix_misencoded_quality_scores' if $phred64;
			   	push @gatk_commands, $realign_indel_cmd;
			   	my $base_recal_cmd = "$java_bin -jar $gatkjar $base_recal_args -R $genome_fasta -L $intervals_file -I GATKOUT.sorted_dupmarked_readgrp_realigned.".$chr.".bam -knownSites $dbsnp_file -knownSites $mills_indel_file -o GATKOUT.recal.".$chr.".table";
			   	#$base_recal_cmd .= ' --fix_misencoded_quality_scores' if $phred64;
			   	push @gatk_commands, $base_recal_cmd;
			   	my $print_reads_cmd = "$java_bin -jar $gatkjar $print_reads_args -R $genome_fasta -L $intervals_file -I GATKOUT.sorted_dupmarked_readgrp_realigned.".$chr.".bam -BQSR GATKOUT.recal.".$chr.".table -o GATKOUT.sorted_dupmarked_readgrp_realigned_recal.".$chr.".bam";
			   	#$print_reads_cmd .= ' --fix_misencoded_quality_scores' if $phred64;
			   	push @gatk_commands, $print_reads_cmd;
			   	my $haplo_caller_cmd = "$java_bin -jar $gatkjar $haplo_caller_args -R $genome_fasta -L $intervals_file -I GATKOUT.sorted_dupmarked_readgrp_realigned_recal.".$chr.".bam -o GATKOUT.haplotypecaller.".$chr.".g.vcf.gz";
			   	#$haplo_caller_cmd .= ' --fix_misencoded_quality_scores' if $phred64;
			   	push @gatk_commands, $haplo_caller_cmd;
			   	my $genotype_gvcfs_cmd = "$java_bin -jar $gatkjar $genotype_gvcfs_args -R $genome_fasta -L $intervals_file -V GATKOUT.haplotypecaller.".$chr.".g.vcf.gz -o $vcf_file";
			   	#$genotype_gvcfs_cmd .= ' --fix_misencoded_quality_scores' if $phred64;
			   	push @gatk_commands, $genotype_gvcfs_cmd;
			   
			    
				#Create a record of vcf being created; need to allow mpileup jobs to finish before proceeding to next pipeline steps
				my $record_command = "$record_snv_bin -vcf $vcf_file -chr $chr -step_name single_vcf -qsub_file $snv_call_qsub -runid RUNID"; #RUNID gets replaced by submit_snvs.pl
				push @gatk_commands, $record_command;
			}
	
	
		    $self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,
		    										-qsub_file=>$snv_call_qsub_tmp,
		    										-commands=>\@gatk_commands,
		    										-cpus=>1,
		    										-mem=>'4Gb',
		    										-walltime=>'24:00:00',
		    										-modules=>$module);
		}
    }
}


#Generate the bwa align lane to sort bam commands
sub align_lane_qsub 
{
	my $self = shift;
	
	if (!defined $self->{source_name}) {
		modules::Exception->throw("ERROR: Must call create_sample_xml via pipeline object");
	}
	
    my %args = @_;
	
    my @required_args = (
			             '-read1',
			             '-read2',
			             '-lane_name',
						 '-lane_id',
						 '-sample_name',
						 '-outdir'
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $read1 = $args{-read1};
    my $read2 = $args{-read2};
    my $lane_name = $args{-lane_name};
    my $lane_id = $args{-lane_id};
    my $outdir = $args{-outdir};
    my $sample_name = $args{-sample_name};
    my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_name);
    
    #These are from the xml files
    my $svndir = $self->_svndir();
 	my $threadnum = $self->_threadnum();
 	my $max_memory = $self->_cluster_xml->read('common','qsub_vars','max_memory');
 	
 	#Account for samtools poor memory estimation
 	my ($mem_number) = $max_memory =~ /(\d+)/;
 	my $max_sam = $mem_number * 1000000000 * 0.33; 
 	
    my $qsub_dir = defined $args{-qsub_dir}?$args{-qsub_dir}:$self->_qsub_dir().'/'.$self->{source_name}.'/lanes';
   	if ( !-d $qsub_dir ) {
   		system("mkdir -p $qsub_dir");	
   	}
   	my $index = $self->_cluster_xml->read($self->{source_type},'svn','bwa_index');
    

    my %commands = ();
	my $encoding;
    
    
	if ($args{-skip_quality}) {
		$encoding = 'phred33';	
	} else {
    	$encoding = modules::QualityEncoding::encoding(-readfile=>"$read1", -reads_to_check=>100000);
	}
    
    #If the files need decompression first (only bz2 need this so shouldn't happen anymore)
    if (defined $args{-commands}) {
    	%commands = %{$args{-commands}};
		for my $command (keys %commands) {
			push @{$commands{bwa}}, $command;
		}
    	$read1 =~ s/.bz2//; #Change the read files to be the uncompressed versions
    	$read2 =~ s/.bz2//;
    }
    
    #Create the file name; all match lane name
    my $read1_sai = $lane_name.'_r1.sai';
    my $read2_sai = $lane_name.'_r2.sai';
    my $sam = $lane_name.'.sam';
	(my $bam = $sam) =~ s/.sam$//; 
	my $aln_args;
	
	
	if ($encoding eq 'phred33') {
		$aln_args = "aln";
	} else {
		$aln_args = "aln -I";
	}

    my $bwa_qsub_file = $lane_name.'.bwa.qsub';
    my $sam_qsub_file =  $lane_name.'.sam.qsub';
	my $bam_qsub_file =  $lane_name.'.bam.qsub';
	my $bam_stats_qsub_file = $lane_name.'.bam_stats.qsub';
	my $align_lane_qsub_file = $lane_name.'.align_lane.qsub';
	#my $qsub_wrapper_file =  $lane_name.'.wrapper.qsub';
	my @jobs = ($bwa_qsub_file,$sam_qsub_file,$bam_qsub_file,$bam_stats_qsub_file,$align_lane_qsub_file);

	my $bwa_full_qsub = $qsub_dir . '/' . $bwa_qsub_file;
	my $sam_full_qsub = $qsub_dir . '/' . $sam_qsub_file;
	my $bam_full_qsub = $qsub_dir . '/' . $bam_qsub_file;
	my $bam_stats_full_qsub = $qsub_dir .'/'.$bam_stats_qsub_file;
	my $samtools_bin = $self->_samtools();
	my $bwa_bin = $self->_bwa();
	my $index_short = $self->_pipe_xml->read($self->{source_type},'annotation_version','ref_genome');

	my $bwamem = $bwa_bin.'_mem';
	my $rg_str = '"@RG\tID:'.$lane_name.'\tSM:'.$sample_name.'\tPL:ILLUMINA"'; #default RG fields

	# Implementing Aaron Statham's KCCG fastq lane naming conventions in lookup file
	if($lane_name=~/^MGRBp1_sample(\d+)_/) {
		my $sampleNum=$1;
		my @lookup=();
		open IN,'/g/data3/wq2/MGRB_Phase1_FASTQs/MGRB_Phase1_Translation_Table_NCI_2016-09-04.csv';
		while(<IN>) {
			chomp;
			push @lookup,[split/,/];
		}
		close IN;
		my($sm,$lb,$mc,$pu,$r1,$r2)=@{$lookup[$sampleNum]};
		$rg_str = '"@RG\tID:'.$lane_name.'\tSM:'.$sm.'\tLB:'.$lb.'\tmc:'.$mc.'\tPU:'.$pu.'\tCN:KCCG\tPL:ILLUMINA\tPM:HiSeqX"';
	}

	my $hthreadnum = ($threadnum/2)?$threadnum/2:1;
	push @{$commands{bwa}}, "$bwamem mem -M -R $rg_str -t $threadnum $index $read1 $read2 |  $samtools_bin view -u -S - | /g/data/u86/software/bin/novosort --md --kt -c $hthreadnum -x 6 -m 10G -t \$PBS_JOBFS -i -o $outdir/$bam.bam -";
	push @{$commands{sam}},	"#awk '{if (NF>=15 && substr(\$14,2,8) != substr(\$15,2,8)) {print \$0} else if (NF!=15) {print \$0}}' $outdir/$sam | $samtools_bin view -S -b - | $samtools_bin sort -m$max_sam - $outdir/$bam ";				
	

	push @{$commands{bam}}, "#$samtools_bin index $outdir/$bam.bam";
	push @{$commands{bam_stats}}, "$svndir/scripts/get_bam_stats.pl -bam $outdir/$bam.bam -output $outdir/$bam.stats -sample_name $sample_name";
	push @{$commands{align_lane}}, "$svndir/scripts/align_lanes_db.pl -lane_id $lane_id -phred_quality $encoding -qsub_bam_stats $bam_stats_full_qsub -qsub_bwa $bwa_full_qsub -qsub_sam $sam_full_qsub -qsub_bam $bam_full_qsub";
	
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_next_file=>$sam_qsub_file,-qsub_file=>$bwa_qsub_file,-commands=>$commands{bwa},-cpus=>$threadnum, -mem=>$max_memory, -walltime=>'48:00:00',-jobfs=>'200Gb');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_next_file=>$bam_qsub_file,-qsub_file=>$sam_qsub_file,-commands=>$commands{sam},-cpus=>3, -mem=>$max_memory, -walltime=>'96:00:00',-jobfs=>'100Gb');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_next_file=>$bam_stats_qsub_file,-qsub_file=>$bam_qsub_file,-commands=>$commands{bam},-cpus=>1, -mem=>'4Gb', -walltime=>'4:00:00');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_next_file=>$align_lane_qsub_file,-qsub_file=>$bam_stats_qsub_file,-commands=>$commands{bam_stats},-cpus=>1, -mem=>'1Gb');
	$self->{cluster_obj}->create_single_qsub(-qsub_dir=>$qsub_dir,-qsub_file=>$align_lane_qsub_file,-commands=>$commands{align_lane},-cpus=>1, -mem=>'1Gb');
	#$self->{cluster_obj}->create_qsub_wrapper(-qsub_dir=>$qsub_dir,-qsub_file=>$qsub_wrapper_file,-jobs=>\@jobs);
    
    return ($encoding,"$outdir/$bam.bam");
}

#Create all the individual qsub files and a wrapper for a range of steps
sub create_run_qsubs
{
    my $self = shift;
	
    if (!defined $self->{source_name}) {
	modules::Exception->throw("ERROR: Must call create_run_qsub via PipelineSample object");
    }

    my %args = @_;
    my @required_args = (-steps_xml,
			 -start_step,
			 -sample_name,
			 -pipe_block,
	);

    foreach my $required_arg (@required_args){

	if (! defined $args{$required_arg}){
	    modules::Exception->throw("Required argument [$required_arg] not set");
	}
    }
    
    my $svndir = $self->_svndir();
    my $pipeline_script = $svndir . '/scripts/run_pipeline_steps.pl';
    my $qsub_dir = defined $args{-qsub_dir}?$args{-qsub_dir}:$self->_qsub_dir().'/'.$self->{source_name}.'/runs';
    my $new = '';
    my $start_step = $args{-start_step};
    my $sample_name = $args{-sample_name};
    my $pipe_block = $args{-pipe_block};
    
	
    #Get the steps xml
    my $steps_xml = $args{-steps_xml};
    if ( !-e $steps_xml ) {
	modules::Exception->throw("File $steps_xml doesn't exist");	
    }
    my $steps_config = modules::ConfigXML->new($steps_xml);
    my $steps = $steps_config->read('steps_order', 'step');
    my @steps = ();
    
    #Get the steps list	
    if (ref($steps) eq 'ARRAY'){ # Cope with steps being a single step or an array of steps
    	@steps = @$steps;
    } else {
    	@steps = ($steps);
    }
	
    my %step_numbers = my %step_lookup = ();

    #Get the step order to check last_step is after resume_step
    my $step_count = 0;
    for my $step ( @steps ) {
    	$step_numbers{$step} = $step_count;
    	$step_lookup{$step_count} = $step;
    	$step_count++;
    }
	
    #Check the step arguments exist in the xml
    if (!$self->_step_exists($start_step,\@steps)) {
	modules::Exception->throw("ERROR: Step $start_step doesn't exist in xml file $steps_xml");
    }
	
    my $end_step;
    if (defined $args{-end_step}) {
	$end_step = $args{-end_step};
	if (!$self->_step_exists($end_step,\@steps)) {
	    modules::Exception->throw("ERROR: Step $end_step doesn't exist in xml file $steps_xml");
	}
    } else {
	$end_step = $start_step;
    }
    
    my $first_step_count = $step_numbers{$start_step};
    my $last_step_count = $step_numbers{$end_step};
	
    my @dependent_jobs = ();
	
    for ( my $count = $first_step_count ; $count <= $last_step_count ; $count++ ) {
	my $local_step = $step_lookup{$count};
	my $qsub_file = $sample_name.'.pipe'.$pipe_block.'.'.$local_step.'.qsub';
	
	my $last_block_step = $local_step eq $end_step?1:0;
	
	#Check if we're dealing with the first step
	my $first_step = $steps_config->exists('steps', 'step', $local_step, 'first_step')?1:0;
	my $pipe_args = "-xml $steps_xml";	
	if ($first_step) {
	    $pipe_args .= " -last_step $local_step";
	} else {
	    $pipe_args .= " -resume $local_step -latest_sample $sample_name -last_step $local_step";
	}
	
	if ( !-d $qsub_dir ) {
	    system("mkdir -p $qsub_dir");	
	}
	
	my @commands = ("$pipeline_script $pipe_args");
	
	#now add the block of code to confirm the step ran; prevents submitting new jobs after failed step; 
	#One case it will fail in is where the step is being run and the job dies before connecting to database
	my $check_command = "\noutput=\$($svndir/scripts/check_step_ran.pl -sample $sample_name -step $local_step -xml $steps_xml 2>&1)";
	push @commands, $check_command, 'regex="[EXCEPTION]"','if [[ $output =~ $regex ]]',"then\n\texit\nfi\n";
	
	
	#Get the resources required for the step in question from the xml
	my $mem = $steps_config->exists('steps', 'step', $local_step, 'mem')? $steps_config->read('steps', 'step', $local_step, 'mem'):"4GB";
	my $cpus = $steps_config->exists('steps', 'step', $local_step, 'cpus')? $steps_config->read('steps', 'step', $local_step, 'cpus'):"1";
	my $walltime = $steps_config->exists('steps', 'step', $local_step, 'walltime')? $steps_config->read('steps', 'step', $local_step, 'walltime'):"24:00:00";	
	my $queue = $steps_config->exists('steps', 'step', $local_step, 'copyq')?'copyq':'normal';
	if ($last_block_step) {
	    if ($steps_config->exists('steps','step',$local_step,'module')) {

		my $module = $steps_config->read('steps','step',$local_step,'module');

		$self->{cluster_obj}->create_single_qsub(-qsub_dir  => $qsub_dir,
							 -qsub_file => $qsub_file,
							 -commands  => \@commands,
							 -cpus      => $cpus,
							 -mem       => $mem,
							 -walltime  => $walltime,
							 -queue     => $queue,
							 -modules   => $module);
	    } else {
		$self->{cluster_obj}->create_single_qsub(-qsub_dir  => $qsub_dir,
							 -qsub_file => $qsub_file,
							 -commands  => \@commands,
							 -cpus      => $cpus,
							 -mem       => $mem,
							 -walltime  => $walltime,
							 -queue     => $queue);
	    }
	} else {
	    #here we need to add the next job submission
	    my $next_job = $step_lookup{$count+1};
	    my $next_qsub_job = $sample_name.'.pipe'.$pipe_block.'.'.$next_job.'.qsub';

	    if ($steps_config->exists('steps','step',$local_step,'module')) {

		my $module = $steps_config->read('steps','step',$local_step,'module');

		$self->{cluster_obj}->create_single_qsub(-qsub_dir       => $qsub_dir,
							 -qsub_next_file => $next_qsub_job,
							 -qsub_file      => $qsub_file,
							 -commands       => \@commands,
							 -cpus           => $cpus,
							 -mem            => $mem,
							 -walltime       => $walltime,
							 -queue          => $queue,
							 -modules        => $module);

	    } else {

		$self->{cluster_obj}->create_single_qsub(-qsub_dir       => $qsub_dir,
							 -qsub_next_file => $next_qsub_job,
							 -qsub_file      => $qsub_file,
							 -commands       => \@commands,
							 -cpus           => $cpus,
							 -mem            => $mem,
							 -walltime       => $walltime,
							 -queue          => $queue);

	    }
	}
	#push @dependent_jobs, $qsub_file;
	last if $local_step eq $end_step;
    }
    
    #Now create the wrapper if it's more than one command
    #my $qsub_wrapper_file = $sample_name.'.pipe'.$pipe_block.'.wrapper.qsub';
    #$self->{cluster_obj}->create_qsub_wrapper(-qsub_dir=>$qsub_dir,-qsub_file=>$qsub_wrapper_file,-jobs=>\@dependent_jobs) if $end_step ne $start_step;
}

#Validate line
sub validate_source_line {
	
	my $self = shift;
	
	if (!defined $self->{source_name}) {
		modules::Exception->throw("ERROR: Must call create_run_qsub via pipeline object");
	}
	
    my %args = @_;
    my @required_args = (
			            -line
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
	
	#format of lines samples
	#human_single,Melanoma,A11,R1,RAM,genome,A11_LCL,affected
	#human_related,2,3,'Rare disease',RARE_cohort4,R1,BRF,exome,2723930,2773179:2822150,RARE_cohort4_proband,RARE_cohort4_unaffected_father:RARE_cohort4_unaffected_mother,proband,father:mother,male,both
	
	my $line = $args{-line};
	my @fields = split(",",$line);
	my $error = 0;
	
	#Number of entries expected 
	my %entry_number = (
					'human_related' => 14,
					'human_single' => 8
					);
	
	if (@fields != $entry_number{$fields[0]}) {
		$error = "Line entries and expected numbers don't match";
	}
	
	(my $project_name = $fields[1]) =~ s/'//g;
	my ($proj_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
	if (!$proj_obj) {
		$error = "Can't find project $project_name in db";
	}
	
	if ($fields[5] !~ /exome/ && $fields[5] !~ /genome/ && $fields[5] !~ /targeted/) {
		$error = "$fields[5] must be exon, genome, or targeted";
	}
	
	if ($fields[0] eq 'human_single') {
		(my $affected = $fields[7]) =~ s/ //g;
		if ($affected ne 'affected' && $affected ne 'unaffected') {
			$error = "$affected must be affected or unaffected";
		}
		
	} elsif ($fields[0] eq 'human_related') {
		(my $affected = $fields[13]) =~ s/ //g;
		if ($affected ne 'affected' && $affected ne 'both') {
			$error = "$affected must be affected or both";
		}
		my @affected_names = split(':',$fields[6]);
		my @affected_readdirs = split(':',$fields[8]);
		my @affected_relations = split(':',$fields[10]);
		my @unaffected_names = split(':',$fields[7]);
		my @unaffected_readdirs = split(':',$fields[9]);
		my @unaffected_relations = split(':',$fields[11]);
		
		
		
		if (@affected_readdirs != @affected_relations) {
			$error = "Number of affected entries don't match"
		}
		if (@unaffected_readdirs != @unaffected_relations) {
			$error = "Number of unaffected entries don't match"
		}
		if (@affected_readdirs != @affected_names) {
			$error = "Number of affected entries don't match"
		}
		if (@unaffected_readdirs != @unaffected_names) {
			$error = "Number of unaffected entries don't match"
		}
		
	} else {
		$error = "Must be human_related or human_single";
	}
	return $error;
}


return 1;
