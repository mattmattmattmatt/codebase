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
use modules::SystemCall;
use Pod::Usage;
use Cwd;

use vars qw(%OPT);


GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"project_name=s",	# eg) 'Rare disease','APOSLE' etc
	   		"family_tree=s",		# eg) "Request169-FamilyTree.txt"
	   		"importdir=s",		# initial import directory
	   		"run",				# execute commands
	   		"family_tree_dir=s",
			"get_last_source_name",
	   		"new_source_name_related=s",
			"new_source_name_single=s",
			);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

prepare_samples.pl -project_name APOSLE -family_tree Request169-FamilyTree.txt -importdir /home/vicky/tmp_importdir -run [options]
Required flags: family_tree, project_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

prepare_samples.pl -> prepare samples to update '../conf/sample_info.csv' and create read data directories

=head1 DESCRIPTION

June 06, 2014

a script that reads family tree info from APF database and updates "sample_info.csv" for a new human sample and creates read data directories

=head1 AUTHOR

Vicky Cho

=head1 EXAMPLE

.prepare_samples.pl -source_type human_related -project_name APOSLE -family_tree Request169-FamilyTree.txt

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

my $pipe_config = modules::Pipeline::get_pipe_conf();
my @possible_relations = split(",", $pipe_config->read('common','related_types'));

my $clus_conf = modules::Pipeline::get_cluster_conf();

my $sample_csv = modules::Pipeline::get_sample_conf();
#my $sample_csv = "/home/vicky/workspace/v2.1/conf/human_sample_info_test.csv";

#print "Sample config file to update: $sample_csv\n\n";

my $project_name = $OPT{project_name} ? $OPT{project_name}:_prompt_option(-name => 'project name', -default => 'APOSLE');

my ($proj_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
if (!$proj_obj) {
	modules::Exception->throw("ERROR: Can't find project $project_name in db \n");
}

my $latest_external_name_single = get_latest_external_source_name('human_single_gatk',$project_name);
my $latest_external_name_related = get_latest_external_source_name('human_related_gatk',$project_name);

print "The latest human_single source name for project $project_name is $latest_external_name_single \n";
print "The latest human_related source name for project $project_name is $latest_external_name_related \n\n";

if($OPT{get_last_source_name}){
	exit;
}

my $family_tree;
my $family_tree_dir = $OPT{family_tree_dir} ? $OPT{family_tree_dir}:'/g/data/u86/snv_pipeline_runs/v2.3_human_genome/apf_family_tree/';
if($OPT{family_tree}){ 
	$family_tree = $OPT{family_tree};
}else {
	print "\nListing family tree files: ls $family_tree_dir* ...\n";
	
	for my $file (glob $family_tree_dir.'*') {
		my @filepath = split('/',$file);
		print "\t".$filepath[scalar(@filepath)-1]."\n";
	}
	
	print "\n";	
	
	$family_tree = _prompt_option(-name => 'Family tree file');
	$family_tree = $family_tree_dir.$family_tree;

#print $family_tree;
}

my $import_reads = _prompt_option(-name => 'Import read files?', -default => 'Yes');
$import_reads = $import_reads =~ /Yes|yes|YES/?1:0;

my $read1_pattern = _prompt_option(-name => 'read1 pattern', -regex => '(R1|_1)', -default => 'R1');
my $sequencing_centre = _prompt_option(-name => 'sequencing centre', -regex => '(BRF|AGRF|RAM|BGI|Macrogen|Novogene)', -default => 'BRF');
my $sequence_type = _prompt_option(-name => 'sequence type', -regex => '(exome|genome|targeted)', -default => 'exome');

my %samples = ();
my @header = ();

open(FAMILY_TREE,$family_tree) || modules::Exception->throw("Can't open file $family_tree");

#$/="\n";	## fix chomp acting wierd..
$/="\r\n";	## fix chomp acting wierd..

while(<FAMILY_TREE>){
	if($. == 1 ){ 
		chomp;
		@header = split(/\|/); ### RequestId|SampleId|IndividualName|Proband|Gender|Affected|FatherName|MotherName|FamilyId

	}else{
		chomp;
		s/,/_/g; ## some samples have multiple reqids or HUM ids separated by ",", replace separater to "_"
		
		my %info = ();		
		@info{@header}= split(/\|/);
		#print Dumper @info{@header};
		$samples{$info{IndividualName}} = {%info};		###### 'samples' hash of hashes, key = IndividualName
		
		if(!defined $info{FamilyId}){
			$samples{$info{IndividualName}}{FamilyId}="";
			#modules::Exception->throw("FamilyId has not been assigned");
		}
	}
}


my $nS = scalar keys %samples;
print "\nTotal number of samples in this request = $nS \n";

my @affected_sample_names;
my %families;				# group families using FamilyId
my %single_samples;			# request may contain single samples 
my %proband_sex;			# maps familyID to proband sex

#print Dumper \%samples;

for my $i ( keys %samples) {

	##-- there is a family within this request when proband exists and when FamilyId doesn't exists!
	if($samples{$i}{Proband} eq 'Yes' && $samples{$i}{FamilyId} ne ""){
	
		my $proband_familyId = $samples{$i}{FamilyId};
		$families{$proband_familyId}{$samples{$i}{IndividualName}}='proband';		
		if(!defined($samples{$i}{FatherName}) | $samples{$i}{FatherName} ne ''){
			$families{$proband_familyId}{$samples{$i}{FatherName}}='father';
		}		
		if(!defined($samples{$i}{MotherName}) | $samples{$i}{MotherName} ne ''){
			$families{$proband_familyId}{$samples{$i}{MotherName}}='mother';
		}			
		$proband_sex{$proband_familyId} = lc($samples{$i}{Gender});		
		
		print "\n##### Processing family $proband_familyId : proband sample $samples{$i}{IndividualName} \n";
					
		##-- find other members of the family using FamilyId
		##-- Round1: Searching sister/brother/sibling/daugher/son/child of the proband 
		for my $j (keys %samples){
			if($samples{$j}{FamilyId} eq $proband_familyId){
				print "\n         ##### Found a member of family $proband_familyId - processing sample $samples{$j}{IndividualName} \n";
								
				if(!exists $families{$proband_familyId}{$samples{$j}{IndividualName}}){					
					#sister/brother/siblings if $samples{$j} share either father or mother; 
					#daughter/son/child if $samples{$j} has proband as one of the parents;
					if( defined($families{$proband_familyId}{$samples{$j}{FatherName}}) & $samples{$j}{FatherName} ne ''){
						if($families{$proband_familyId}{$samples{$j}{FatherName}} eq 'father'){
							if($samples{$j}{Gender} eq 'Female'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='sister';
								next;
							}elsif($samples{$j}{Gender} eq 'Male'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='brother';
								next;
							}elsif( !defined($samples{$j}{Gender}) | $samples{$j}{Gender} eq ''){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='sibling';
#print Dumper \%families;	
								next;
							}								
						}elsif($families{$proband_familyId}{$samples{$j}{FatherName}} eq 'proband'){
							if($samples{$j}{Gender} eq 'Female'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='daughter';
								next;
							}elsif($samples{$j}{Gender} eq 'Male'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='son';
								next;
							}elsif( !defined($samples{$j}{Gender}) | $samples{$j}{Gender} eq ''){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='child';
								next;
							}							
						}else{
							$families{$proband_familyId}{$samples{$j}{IndividualName}}='other';
						}		
					}
					
					if(defined($families{$proband_familyId}{$samples{$j}{MotherName}}) & $samples{$j}{MotherName} ne ''){
						if($families{$proband_familyId}{$samples{$j}{MotherName}} eq 'mother'){
							if($samples{$j}{Gender} eq 'Female'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='sister';
								next;
							}elsif($samples{$j}{Gender} eq 'Male'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='brother';
								next;
							}elsif( !defined($samples{$j}{Gender}) | $samples{$j}{Gender} eq ''){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='sibling';
								next;
							}							
						}
						elsif($families{$proband_familyId}{$samples{$j}{MotherName}} eq 'proband'){
							if($samples{$j}{Gender} eq 'Female'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='daughter';
								next;
							}elsif($samples{$j}{Gender} eq 'Male'){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='son';
								next;
							}elsif( !defined($samples{$j}{Gender}) | $samples{$j}{Gender} eq ''){
								$families{$proband_familyId}{$samples{$j}{IndividualName}}='child';
								next;
							}								
						}
						else{
							$families{$proband_familyId}{$samples{$j}{IndividualName}}='other';
						}
#print Dumper \%families;
					}
					else{					
						$families{$proband_familyId}{$samples{$j}{IndividualName}}='other';
#print Dumper \%families;
					}
					print "         ##### sample $samples{$j}{IndividualName} has defined relationship as $families{$proband_familyId}{$samples{$j}{IndividualName}} \n";
#print Dumper \%families;
				}
			}			
		}
		
		my @family_members = keys %{$families{$proband_familyId}};
				
		##-- Looping through all members with the same FamilyId, work out relationship for 'other' 
		##-- Round2: Searching spouse/grandfather/grandmother/neice/nephew/nibling of the proband 
		foreach my $family_member (@family_members){
			if( $families{$proband_familyId}{$family_member} eq 'other'){
#				print "\n \t\t\t unassigned family member name: $family_member \n";
				
				for my $i (@family_members){
					#--'other' is a spouse of the proband ; there's a child having this member as father/mother
					if(defined($samples{$i}{FatherName}) & defined($families{$proband_familyId}{$samples{$i}{MotherName}})){
						if($samples{$i}{FatherName} eq $family_member & $families{$proband_familyId}{$samples{$i}{MotherName}} eq 'proband'){
							$families{$proband_familyId}{$family_member}='spouse';
						}
					}elsif(defined($samples{$i}{MotherName}) & defined($families{$proband_familyId}{$samples{$i}{FatherName}})){
						if($samples{$i}{MotherName} eq $family_member & $families{$proband_familyId}{$samples{$i}{FatherName}} eq 'proband'){
							$families{$proband_familyId}{$family_member}='spouse';
						}
					}
								
					#--'other' is a grandfather/grandmother ; this member is father/mother of proband's father/mother
					
					if(($families{$proband_familyId}{$samples{$i}{IndividualName}} eq 'father') | ($families{$proband_familyId}{$samples{$i}{IndividualName}} eq 'mother')){
						if($samples{$i}{FatherName} eq $family_member){		  	
					  		$families{$proband_familyId}{$family_member}='grandfather';
						}
						if($samples{$i}{MotherName} eq $family_member){
							$families{$proband_familyId}{$family_member}='grandmother';
						}
					}					
					
					#--'other' is a niece/nephew ; this member has proband's sibling as father/mother
					if(defined($families{$proband_familyId}{$samples{$i}{FatherName}}) & $samples{$i}{FatherName} ne ''){
						if($families{$proband_familyId}{$samples{$i}{FatherName}} eq 'sibling' | $families{$proband_familyId}{$samples{$i}{FatherName}} eq 'brother'){
							
							if($samples{$family_member}{Gender} eq 'Female'){
								$families{$proband_familyId}{$family_member} = 'niece';
							}elsif($samples{$family_member}{Gender} eq 'Male'){
								$families{$proband_familyId}{$family_member} = 'nephew';
							}elsif(!defined($samples{$family_member}{Gender}) | $samples{$family_member}{Gender} eq ''){
								$families{$proband_familyId}{$family_member} = 'nibling';
							}
						}
					}elsif(defined($families{$proband_familyId}{$samples{$i}{MotherName}}) & $samples{$i}{MotherName} ne ''){
						if($families{$proband_familyId}{$samples{$i}{MotherName}} eq 'sibling' | $families{$proband_familyId}{$samples{$i}{MotherName}} eq 'sister'){
							if($samples{$family_member}{Gender} eq 'Female'){
								$families{$proband_familyId}{$family_member} = 'niece';
							}elsif($samples{$family_member}{Gender} eq 'Male'){
								$families{$proband_familyId}{$family_member} = 'nephew';
							}elsif(!defined($samples{$family_member}{Gender}) | $samples{$family_member}{Gender} eq ''){
								$families{$proband_familyId}{$family_member} = 'nibling';
							}	
						}						
					}
				}
			}
		}
		
		
		##-- Looping through all members with the same FamilyId, work out relationship for 'other'
		##-- Round3 - needs to be done after 'grandfather/grandmother' level is defined
		foreach my $family_member (@family_members){		
			if( $families{$proband_familyId}{$family_member} eq 'other'){
#				print "\n \t\t\t unassigned family member name: $family_member \n";
				
				#--'other' is an aunt/uncle/ ; aunt/uncle has grandfather/grandmother as mother/father
				if(defined($samples{$family_member}{FatherName}) & $samples{$family_member}{FatherName} ne ''){
					
					if($families{$proband_familyId}{$samples{$family_member}{FatherName}} eq 'grandfather'){
						if($samples{$family_member}{Gender} eq 'Male'){
							$families{$proband_familyId}{$family_member} = 'uncle';
						}elsif($samples{$family_member}{Gender} eq 'Female'){
							$families{$proband_familyId}{$family_member} = 'aunt';
						}
					}
				}elsif(defined($samples{$family_member}{MotherName}) & $samples{$family_member}{MotherName} ne ''){
					if($families{$proband_familyId}{$samples{$family_member}{MotherName}} eq 'grandmother'){
						if($samples{$family_member}{Gender} eq 'Male'){
							$families{$proband_familyId}{$family_member} = 'uncle';
						}elsif($samples{$family_member}{Gender} eq 'Female'){
							$families{$proband_familyId}{$family_member} = 'aunt';
						}		
					}		
				}
				#--'other' is a grandchild 
				if(defined($samples{$family_member}{FatherName}) & $samples{$family_member}{FatherName} ne ''){
					if($families{$proband_familyId}{$samples{$family_member}{FatherName}} eq 'son' | $families{$proband_familyId}{$samples{$family_member}{FatherName}} eq 'child'){
						if($samples{$family_member}{Gender} eq 'Male'){
							$families{$proband_familyId}{$family_member} = 'grandson';	
						}elsif($samples{$family_member}{Gender} eq 'Female'){
							$families{$proband_familyId}{$family_member} = 'granddaughter';	
						}elsif(!defined($samples{$family_member}{Gender}) | $samples{$family_member}{Gender} eq ''){
							$families{$proband_familyId}{$family_member} = 'grandchild';	
						}
					}
				}elsif(defined($samples{$family_member}{MotherName}) & $samples{$family_member}{MotherName} ne ''){
					if($families{$proband_familyId}{$samples{$family_member}{MotherName}} eq 'daughter' | $families{$proband_familyId}{$samples{$family_member}{MotherName}} eq 'child'){
						if($samples{$family_member}{Gender} eq 'Male'){
							$families{$proband_familyId}{$family_member} = 'grandson';	
						}elsif($samples{$family_member}{Gender} eq 'Female'){
							$families{$proband_familyId}{$family_member} = 'granddaughter';	
						}elsif(!defined($samples{$family_member}{Gender}) | $samples{$family_member}{Gender} eq ''){
							$families{$proband_familyId}{$family_member} = 'grandchild';	
						}	
					}		
				}
			}
		}
		
		
		##-- Looping through all members with the same FamilyId, if there's more than 1 same relation type within the famiy, append a number eg daughter1, daughter2
		
		my %relation_types_counts;
		
		foreach my $family_member (@family_members){
			my $relation = $families{$proband_familyId}{$family_member};
			if(!exists $relation_types_counts{$relation}){
				$relation_types_counts{$relation}=1;
			}else{
				$relation_types_counts{$relation}++; 
			}
			
		}
		
		foreach my $relation_types (keys %relation_types_counts){
			if($relation_types_counts{$relation_types}>1){
				
				my @sample_names_to_add_suffix;
				
				foreach my $family_member (@family_members){
					if( $families{$proband_familyId}{$family_member} eq $relation_types){
						push @sample_names_to_add_suffix, $family_member;
					}	
				}
#				print "######## \n";
#				print Dumper \@sample_names_to_add_suffix;
				
				for(my $count =1; $count <= $relation_types_counts{$relation_types}; $count++ ){
					my $relation_with_suffix = $families{$proband_familyId}{$sample_names_to_add_suffix[$count-1]}.$count;
					$families{$proband_familyId}{$sample_names_to_add_suffix[$count-1]}=$relation_with_suffix;
				}
			}
		}
		
#		print "relation types included in this family are : \n";
#		print Dumper \%relation_types_counts;
		print Dumper \$families{$proband_familyId};
	}

	if($samples{$i}{FamilyId} eq ""){	#if($samples{$i}{FamilyId}=~'single'){
		$single_samples{$samples{$i}{SampleId}} = $samples{$i}{IndividualName};
	}
	
	if($samples{$i}{Affected} eq 'Yes'){
		push(@affected_sample_names, $samples{$i}{SampleId});
	}
	
}

my $source_type;

if (keys %families){	# if there's any families within this request, proceed

	#-- assign new source name; use new_source_name_related as a counter starting from the latest number+1
	my $separator = "\_cohort";
	my $new_source_name_related;
	
	if($latest_external_name_related eq ""){
		if ($OPT{new_source_name_related}){
			$new_source_name_related = $OPT{new_source_name_related};
		} else {
			modules::Exception->throw("ERROR: Provide new_source_name_related to assign source name");
		}			
	} else {
		if($OPT{new_source_name_related}){	#force renaming source name when provided
			$new_source_name_related = $OPT{new_source_name_related};
		} else{
			my @source_name = split($separator,$latest_external_name_related);
			$new_source_name_related = join $separator,$source_name[0],++$source_name[1];
		}	
	}
	
	#-- loop through %families  : for each member of family- 1) mkdir in readdir, 2) mv read files from importdir, then for $family 3) append sample info into $sample_csv file

	for my $family (keys %families){
		my @affected_external_sample_names;
		my @unaffected_external_sample_names;
		my @affected_readdir;
		my @unaffected_readdir;
		my @affected_relations;
		my @unaffected_relations;
		my $no_affected;
		my $proband_sex = (defined $proband_sex{$family} & $proband_sex{$family} ne '')?$proband_sex{$family}:'unknown';

		$source_type ='human_related_gatk';		
#		my $new_source_name_related = $family;
	
#		print "new source name for sample $family is $new_source_name_related \n";
		print "Running family id $family \n";
#		print "proband sex is $proband_sex{$family} \n";
		my $readdir = $clus_conf->read($source_type,'base_directories','base_read_directory');
	#	my $readdir = '/home/vicky/tmp_readdir';

		my @family_members = keys %{$families{$family}};
		
#		print "Members in family $family are :\n";
		
		my $family_dir = $readdir."/".$new_source_name_related;
		if (-d $family_dir){
			print "\nWARNING: Directory $family_dir already exists \n";
		}
							
		if(defined $OPT{run} & !-d $family_dir & $import_reads){
			mkdir($family_dir);	
		}
		
		##-- Looping through all members
		foreach my $family_member (@family_members){	
			if($family_member =~ /fake|FAKE|Fake/){
				next;
			}else{
			
				my $HUMid = $samples{$family_member}{SampleId};
				my $sample_readdir;
				
				if(grep{$_ eq $HUMid} @affected_sample_names){	
					
					$sample_readdir = $new_source_name_related."_affected_".$families{$family}{$family_member};
					
					push @affected_external_sample_names, $family_member;
					push @affected_readdir, $sample_readdir;
					
					my $relation = $families{$family}{$family_member};
					$relation =~ s/\d//g;
					if(!grep{$_ eq $relation} @possible_relations){
						$relation = 'other';
					}
					$relation = @possible_relations ? $relation:'other';
					push @affected_relations, $relation;
					
				}else{
					
					$sample_readdir = $new_source_name_related."_unaffected_".$families{$family}{$family_member};
					
					push @unaffected_external_sample_names, $family_member;
					push @unaffected_readdir, $sample_readdir;
					
					my $relation = $families{$family}{$family_member};
					$relation =~ s/\d//g;
					push @unaffected_relations, $relation;
				}
						
				my $mkdir_outdir ="$readdir/$new_source_name_related/$sample_readdir";
				
				if(-d $mkdir_outdir){
					print "\nWARNING: Directory $mkdir_outdir already exists \n";
				}
				
				if(defined $OPT{run} & !-d $mkdir_outdir & $import_reads){	
					system("mkdir $mkdir_outdir");
				}

				#my $readfile_dirname = "Sample_$HUMid";	#default
				my $readfilename = $family_member;	#default

				#my $mv_command = "mv $OPT{importdir}/$readfile_dirname/* $mkdir_outdir";
				my $mv_command = "mv $OPT{importdir}/$family_member/$readfilename*fq.gz $mkdir_outdir";	

				my @files = glob($OPT{importdir}."/".$family_member."/".$readfilename."*fq.gz");
	
				if($import_reads){
					if(!@files){
						modules::Exception->warning("WARNING: $readfilename*fq.gz does not exist in $OPT{importdir}$family_member\n");
					}else{
					
						print "\n**Check commands for making new directory and moving files for $family_member (will run when import_reads=yes)\nmkdir $mkdir_outdir \n$mv_command\n";				
				
						if(defined $OPT{run} & $import_reads){	
							system $mv_command;
							print "\n**Read files are imported!\n";
						}
					}
				}	
				
			}
		}
	
		my $affected_external_sample_names = @affected_external_sample_names ? join (":", @affected_external_sample_names):'';
		my $unaffected_external_sample_names = @unaffected_external_sample_names ? join (":", @unaffected_external_sample_names):'';
		my $affected_readdir = @affected_readdir ? join (":", @affected_readdir):'';
		my $unaffected_readdir = @unaffected_readdir ? join (":", @unaffected_readdir):'';
		
		my $affected_relations = @affected_relations ? join (":",@affected_relations):'';
		my $unaffected_relations = @unaffected_relations ? join (":",@unaffected_relations):'';
		
		my @all_samples_external_names = @affected_external_sample_names;
		push(@all_samples_external_names,@unaffected_external_sample_names);
		
		my @apf_req_ids; 
		for my $name (@all_samples_external_names){
			push(@apf_req_ids,$samples{$name}{RequestId});
		}
		my $apf_req_ids_full = @apf_req_ids ? "apf_req_ids:".join("^",@apf_req_ids):'';
		
		$no_affected = @unaffected_external_sample_names ? "both":"affected";
		
		my $sampleinfoline = 	"$source_type,'$project_name',$new_source_name_related,$read1_pattern,$sequencing_centre,$sequence_type,'$affected_external_sample_names','$unaffected_external_sample_names','$affected_readdir','$unaffected_readdir','$affected_relations','$unaffected_relations',$proband_sex,$no_affected,$apf_req_ids_full";		

		print "\n\n***sample_info line for cohort $family\n\n";
		print "$sampleinfoline";
		print "\n\n****************************************************\n\n";		
	
		if(defined $OPT{run}){
			open(SAMPLEFILE,">>$sample_csv") || die("Cannot Open File $sample_csv");
			print SAMPLEFILE "$sampleinfoline\n";
			close(SAMPLEFILE);
			print "\n***Sample file $sample_csv is updated!\n";
		}
		
		# increment 1 new_source_name_related for next sample;
		my @source_name = split($separator,$new_source_name_related);
                $new_source_name_related = join $separator,$source_name[0],++$source_name[1];
	}
}



if(keys %single_samples){	# if there's any single sample exists, proceed
	#-- assign new source name; use new_source_name_related as a counter starting from the latest number+1
	my $separator = "\_single";
	my $new_source_name_single;
	
        if($latest_external_name_single eq ""){
                if ($OPT{new_source_name_single}){
                        $new_source_name_single = $OPT{new_source_name_single};
                } else {
                        modules::Exception->throw("ERROR: Provide new_source_name_single to assign source name");
                }
        } else {
		if ($OPT{new_source_name_single}){	#force renaming source name when provided
			$new_source_name_single = $OPT{new_source_name_single};
		} else{
			my @source_name = split($separator,$latest_external_name_single);
			$new_source_name_single = join $separator, $source_name[0],++$source_name[1];
		}
	}


	#--loop through @single_samples - 1)mkdir in readdir, 2)mv read files from importdir, 3)append sample info into $sample_csv file
	#Modified by Aaron and Sean at 7th of March 2017, so the meta information will be sorted properly.
	for my $single_sample (sort keys %single_samples){	# $single_sample = HUMId
		$source_type='human_single_gatk';
		my $apf_req_id = "apf_req_ids:".$samples{$single_samples{$single_sample}}{RequestId};
	
#		my $new_source_name_single = $samples{$single_samples{$single_sample}}{FamilyId};

		print "new source name for sample $single_sample is $new_source_name_single \n";
	
		my $readdir = $clus_conf->read($source_type,'base_directories','base_read_directory');
	#	my $readdir = '/home/vicky/tmp_readdir';
	
		my $sample_readdir = $readdir."/".$new_source_name_single;
		
		if(-d $sample_readdir){
			print "\nWARNING: Directory $sample_readdir already exists\n";
		}
	
		if(defined $OPT{run} & $import_reads){	mkdir($readdir."/".$new_source_name_single);}
		
		my $affected_status = "unaffected";
		if(grep{$_ eq $single_sample} @affected_sample_names){
			$affected_status = "affected";
		}

		#my $readfile_dirname = "Sample_$single_sample";	#default
		my $readfilename = $single_samples{$single_sample};	#default

#		my $mv_command = "mv $OPT{importdir}/$readfile_dirname/* $readdir/$new_source_name_single/";
		my $mv_command = "mv $OPT{importdir}/$readfilename/*fq.gz $readdir/$new_source_name_single/";

		my @files = glob($OPT{importdir}."/".$readfilename."/*fq.gz");		
		
		if($import_reads){
			if(!@files){
				modules::Exception->warning("WARNING: *fq.gz does not exist in $OPT{importdir}$readfilename/\n");
			}else{
		
				print "\n**Check commands for moving read files for $single_samples{$single_sample} (will run when import_reads=yes)\n$mv_command \n";

				if(defined $OPT{run} & $import_reads){	
					system $mv_command;
					print "\n**Read files are imported!\n";
				}
			}
		}
		
		my $sampleinfoline = "$source_type,'$project_name',$new_source_name_single,$read1_pattern,$sequencing_centre,$sequence_type,$single_samples{$single_sample},$affected_status,$apf_req_id";

		print "\n\n**sample_info line for sample $single_sample\n\n";
		print "$sampleinfoline";
		print "\n\n****************************************************\n\n";	
			
		if(defined $OPT{run}){
			open(SAMPLEFILE,">>$sample_csv") || die("Cannot Open File");
			print SAMPLEFILE "$sampleinfoline\n";
			close(SAMPLEFILE);
			print "\n***Sample file $sample_csv is updated!\n\n";
		}
                # increment 1 new_source_name_single for next sample;
                my @source_name = split($separator,$new_source_name_single);
                $new_source_name_single = join $separator,$source_name[0],++$source_name[1];
	}
}


#print Dumper \%families;
#print "\n\n Single samples in this requests are:\n\n";
#print Dumper \%single_samples;
#print Dumper \%samples;





#Get latest external_source_name for a given project and source_name
sub get_latest_external_source_name {
	my $source_type = $_[0];
	my $project_name = $_[1];
	my $latest_source_name = "";
	
	my ($proj_obj) = modules::Adaptors::Project->search(project_name=>$project_name);
	if (!$proj_obj) {
		modules::Exception->throw("ERROR: Can't find project $project_name in db \n");
	}
	my $proj_id = $proj_obj->id;
	my ($source_obj) = modules::Adaptors::Source->search_latest_source_name($proj_id,$source_type);

	if(!$source_obj){
		print "There hasn't been any source created for this project for type $source_type \n";
		$latest_source_name = "";
	}else{

		$latest_source_name = $source_obj->external_source_name;
	}
	return $latest_source_name;
}

# Repeat _get_option prompt until valid regex
sub _prompt_option {

	my %args = @_;
    my @required_args = (
			             -name
						 );

    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $got_valid_input = 0;
    my $user_input;
    
    while(!$got_valid_input){
		$user_input = _get_option(	-name => $args{-name},
									-default => $args{-default},
									-format => $args{-format},
									-required => $args{-required},
									-possibilities => $args{-possibilities},
									-yes => $args{-yes});
		
		if (defined $args{-regex} && (!defined $user_input || (defined $user_input && $user_input !~ /^$args{-regex}$/i)) ) {
			print "Please check your input and try again ...\n";
		} else {
			$got_valid_input = 1;
		}
	}	
	return $user_input;
}

# get option from STDIN - Matt
sub _get_option {    
    my %args = @_;   
    my $variable = $args{-default};
    my $name = $args{-name};
    my @possibilities = ();
    if (defined $args{-possibilities}) {@possibilities = @{$args{-possibilities}};}
    my $format_str = defined $args{-format}?$args{-format}:'';
    my $not_empty = defined $args{-required}?1:0;   
    my $user_input;
    #If passed on command line
    if (defined $OPT{$name}) {
		$variable = $OPT{$name};
	} else {
		#Otherwise enter the value for STDIN
		if ($variable) {
			print "$name $format_str [default=$variable]? ";
		} else {	
			print "$name $format_str? ";
		}
		$user_input = <STDIN>;
		chomp ($user_input);

		if ($user_input) {
			$variable = $user_input;
		} 
		if (!$variable && $not_empty) {
			modules::Exception->throw("ERROR: $name cannot be empty");
		}
		if($variable =~ /^y$/i || $variable =~ /^yes$/i){$variable='yes';}
		if($variable =~ /^n$/i || $variable =~ /^no$/i){$variable='no';}
	}
	if (@possibilities) {
		my $acceptable_value = 0;
		for my $acceptable ( @possibilities ) {
		  	if ($acceptable eq $variable) {
		  		$acceptable_value = 1;
		  	}
		}
		if (!$acceptable_value) {
			my $acceptable_str = join(" OR ",@possibilities);
			modules::Exception->throw("ERROR for input: $name must be $acceptable_str");
		}
	}
		
	return $variable;	
}


