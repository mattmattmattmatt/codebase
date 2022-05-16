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
use modules::Adaptors::Lane;
use modules::Adaptors::Human_Related_Sample;
use modules::Adaptors::Human_Single_Sample;
use modules::SystemCall;
use Pod::Usage;
use Cwd;

use vars qw(%OPT);


GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"single_source=s",
	   		"only_qsubs",
	   		"submit",
	   		"add",
	   		"update_db",
	   		"no_mdss",
	   		"test"
	   		);
	   		
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

process_samples.pl -submit submit_jobs -add add_new_samples_to_db -single_source source_name -update_db update_database -sample_config xml_steps_file(default=../conf/human_sample_info.csv) -only_qsubs only_update_qsubs -test add_test_flag_for_trunk [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

process_sample.pl -> Manage samples from config file

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./process_samples.pl -single_sample RA_single1 -add -submit

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

my $sys_call = modules::SystemCall->new();

my $qsubs = defined $OPT{only_qsubs}?1:0;
my $add = defined $OPT{add}?1:0;
my $submit = defined $OPT{submit}?1:0;
my $debug = defined $OPT{single_source}?$OPT{single_source}:0;
my $db = defined $OPT{update_db}?1:0;

#just update qsubs
if ($qsubs && ($add || $submit)) {
	modules::Exception->throw("ERROR: Can't run qsub update with -add or -submit");
} 

#update db flag
if ($db && ($add || $submit)) {
	modules::Exception->throw("ERROR: Can't run qsub update with -add or -submit");
} 

if ($submit && !$add) {
	modules::Exception->throw("ERROR: Can't run -submit without -add");
}

my $add_single_sample = $svndir . "/scripts/add_single_human_sample.pl";
my $add_multi_sample = $svndir . "/scripts/add_cohort_sample.pl";


if (!-x $add_single_sample) {
	modules::Exception->throw("ERROR: Can't locate executable $add_single_sample");
} 
if (!-x $add_multi_sample) {
	modules::Exception->throw("ERROR: Can't locate executable $add_multi_sample");
}


my $sample_csv = modules::Pipeline::get_sample_conf();

my %ENTRY_NUMBER = ('human_related'     => 14,
		    		'human_related_gatk' => 14,	
                    'human_single'      => 8,
                    'human_single_gatk' => 8,
                    'human_related_gatk_apf' => 15,
					'human_single_gatk_apf' => 9
					);

open(SAMPLES,$sample_csv) || modules::Exception->throw("Can't open file $sample_csv");

my @add_commands = ();
my @commands = ();
#add_cohort_sample.pl -affected_readdirs GERMAN_cohort1_proband,GERMAN_cohort1_affected_sibling -unaffected_readdirs GERMAN_cohort1_unaffected_father,GERMAN_cohort1_unaffected_mother,GERMAN_cohort1_unaffected_brother -affected_sample_external_names EFD17102007,JAD14062010 -unaffected_sample_external_names Father-KD26041972,mother-NP03051972,ABD12092004 -affected_relations proband,sibling -unaffected_relations father,mother,brother -sequencing_centre BRF -skip_quality -only_qsubs
#add_single_human_sample.pl -project_name 'German auto-immune' -read1_pattern R1 -sequencing_centre BRF -external_name GERMAN_single8 -external_sample_name NC41199 -skip_quality -only_qsubs

SAMPLE:
while (<SAMPLES>) {
	next if /^#/;
	next unless /,/;
	chomp;
	if (my $error = &Line_Problem($_)) {
		modules::Exception->throw("ERROR $error with line $_");
	}
	my @fields = split(',');

	my $command;
	(my $source_name = $fields[2]) =~ s/'//g;
	
	if ($debug) {
		next unless $source_name eq $debug;
	}
	
	my ($source_obj) = modules::Adaptors::Source->search(external_source_name=>$source_name);
	
	
	
	if ($db && !$source_obj) {
		next SAMPLE; #Skip as not created yet
	} 

	if ($fields[0] eq 'human_single' || $fields[0] eq 'human_single_gatk') {

		$command = join(" ",
							$add_single_sample,
							'-project_name',
							$fields[1],
							'-external_name',
							$fields[2],
							'-read1_pattern',
							$fields[3],
							'-sequencing_centre',
							$fields[4],
							'-sequence_type',
							$fields[5],
							'-external_sample_name',
							$fields[6]
							);
		if ($fields[7] eq 'unaffected') 
		{
			$command .= " -unaffected";
		}
		
   	if ($fields[0] eq 'human_single_gatk') 
   	{
    	$command .= ' -gatk';
		}

		if($fields[1] =~ /APOSLE/)
		{
			my $steps_xml;
			if($fields[0] eq 'human_single_gatk')
			{
				$steps_xml = "$svndir/conf/human_single_gatk_APOSLE.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /German auto-immune/)
		{
			my $steps_xml;
			if($fields[0] eq 'human_single_gatk') 
			{
				$steps_xml = "$svndir/conf/human_single_gatk_GERMAN.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /CACPI/) #marcin
		{
			my $steps_xml;
			if($fields[0] eq 'human_single_gatk')
			{
				$steps_xml = "$svndir/conf/human_single_gatk_CACPI.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /CPIEEU/) #marcin
		{
			my $steps_xml;
			if($fields[0] eq 'human_single_gatk')
			{
				$steps_xml = "$svndir/conf/human_single_gatk_CPIEEU.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /CPINJ/) #marcin
		{
			my $steps_xml;
			if($fields[0] eq 'human_single_gatk')
			{
				$steps_xml = "$svndir/conf/human_single_gatk_CPINJ.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}






	
		if(exists $fields[8])
		{
			if($fields[8] =~ /apf_req_id/)
			{
				my $apf_req_id = (split(/:/, $fields[8]))[1];
				$command .= " -apf_req_id $apf_req_id";
			}
		}

		if ($db) {

			(my $format_sample_name = $fields[6]) =~ s/'//g;
			my @sample_objs = modules::Adaptors::Sample->search(external_sample_name=>$format_sample_name);
			my $sample_obj;
			for my $tmp_obj (@sample_objs) {
				if ($tmp_obj->sample_name =~ /$fields[2]/) {
					$sample_obj = $tmp_obj;
				}
			}
			if (!$sample_obj) {
				modules::Exception->throw("ERROR: no sample obj for external name $fields[6]");
			}
			my $sample_name = $sample_obj->sample_name;
			(my $format_proj_name = $fields[1]) =~ s/'//g;
			my ($proj_obj) = modules::Adaptors::Project->search(id=>$source_obj->project_id);
			if ($proj_obj->project_name !~ $format_proj_name) {
				#update project
				my ($new_proj_obj) = modules::Adaptors::Project->search(project_name=>$format_proj_name);
				if (!$new_proj_obj) {
					modules::Exception->throw("ERROR: Can't retrieve proj_obj for $fields[1]");
				}
				
				&update_db('project_id',$new_proj_obj->id,$source_obj,$sample_name);
			}
			my @lanes = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);
			for my $lane_obj (@lanes) {
				if ($lane_obj->sequencing_centre ne $fields[4]) {
					#update sequencing centre
					&update_db('sequencing_centre',$fields[4],$lane_obj,$sample_name) if $lane_obj->lane_name =~ /$sample_name/;
				}
			}
			
			if ($sample_obj->sequence_type ne $fields[5]) {
				#update sequence type
				&update_db('sequence_type',$fields[5],$sample_obj,$sample_name);
			}
			
			#Stuff specific to human_single_samples	
			my ($human_single_obj) = modules::Adaptors::Human_Single_Sample->search(sample_id=>$sample_obj->id);
			if (!$human_single_obj) {
				print Dumper \@fields;
				print Dumper $sample_obj;
				my $id = $sample_obj->id;
				modules::Exception->throw("ERROR: Can't retrieve human_single_sample for $source_name with id $id and sample_name $sample_name");
			}
		
			my $affected_db = $sample_obj->sample_type =~ /_affected/?1:0;
			if ($fields[7] eq 'affected' && !$affected_db) {
				#update to affected
				&update_db('affected',1,$human_single_obj,$sample_name);
				&update_db('sample_type','single_affected',$sample_obj,$sample_name);
			} elsif ($fields[7] eq 'unaffected' && $affected_db) {
				#update to unaffected
				&update_db('affected',0,$human_single_obj,$sample_name);
				&update_db('sample_type','single_unaffected',$sample_obj,$sample_name);
			}
			next;	
		}
		
		
	} else { # everything else -- incl.  human_related ...

		(my $affected_external_names = $fields[6]) =~ s/:/,/g;
		(my $affected_readdirs =$fields[8]) =~ s/:/,/g;
		(my $affected_relations =$fields[10]) =~ s/:/,/g;
		
		my $unaffected_external_names;
		my $unaffected_readdirs;
		my $unaffected_relations;
		my $unaffected = 0;
		if ($fields[13] =~ /both/) {
			$unaffected = 1;
			($unaffected_external_names = $fields[7]) =~ s/:/,/g;
			($unaffected_readdirs =$fields[9]) =~ s/:/,/g;
			($unaffected_relations =$fields[11]) =~ s/:/,/g;
			
			if ($db) {
				$fields[7] =~ s/'//g;	
				my @unaffected_external_names = split(":",$fields[7]);
				$fields[11] =~ s/'//g;
				my @unaffected_relations = split(":",$fields[11]);
				my $count = 0;	
				
				for my $unaff_name (@unaffected_external_names) {
					$unaff_name =~ s/'//g;
					my @sample_objs = modules::Adaptors::Sample->search(external_sample_name=>$unaff_name);
					my $unaff_sample_obj;
					
					for my $tmp_obj (@sample_objs) {
						my $tmp_name = $tmp_obj->sample_name;
						if ($tmp_obj->sample_name =~ /$fields[2]/) {
							$unaff_sample_obj = $tmp_obj;
						}
					}

					if (!$unaff_sample_obj) {
						modules::Exception->throw("ERROR: no sample obj for external name $fields[6]");
					}
					my $unaff_sample_name = $unaff_sample_obj->sample_name;
					(my $format_proj_name = $fields[1]) =~ s/'//g;
					my ($proj_obj) = modules::Adaptors::Project->search(id=>$source_obj->project_id);
					if ($proj_obj->project_name !~ $format_proj_name) {
						#update project
						my ($new_proj_obj) = modules::Adaptors::Project->search(project_name=>$format_proj_name);
						if (!$new_proj_obj) {
							modules::Exception->throw("ERROR: Can't retrieve proj_obj for $fields[1]");
						}				
						
						&update_db('project_id',$new_proj_obj->id,$source_obj,$unaff_sample_name);
						
					}
					my @lanes = modules::Adaptors::Lane->search(sample_id=>$unaff_sample_obj->id);
					for my $lane_obj (@lanes) {
						if ($lane_obj->sequencing_centre ne $fields[4]) {
							#update sequencing centre
							&update_db('sequencing_centre',$fields[4],$lane_obj,$unaff_sample_name) if $lane_obj->lane_name =~ /$unaff_sample_name/;
							
						}
					}
					
					if ($unaff_sample_obj->sequence_type ne $fields[5]) {
						#update sequence type
						&update_db('sequence_type',$fields[5],$unaff_sample_obj,$unaff_sample_name);
					}
					
					#Stuff specific to human_related_samples	
					my ($human_related_obj) = modules::Adaptors::Human_Related_Sample->search(sample_id=>$unaff_sample_obj->id);
					if (!$human_related_obj) {
						modules::Exception->throw("ERROR: Can't retrieve human_related_sample for $source_name");
					}
					my $affected_db = $unaff_sample_obj->sample_type =~ /_affected/?1:0;
					if ($affected_db) {
						#Supposed to be unaffected
						&update_db('affected',0,$human_related_obj,$unaff_sample_name);
						&update_db('sample_type','related_unaffected',$unaff_sample_obj,$unaff_sample_name);
					}
					
					if ($human_related_obj->relation ne $unaffected_relations[$count]) {
						&update_db('relation',$unaffected_relations[$count],$human_related_obj,$unaff_sample_name);
					}
					
					$count++;
				}		
			}
		} 
		
		if ($db) {
			
		 	$fields[6] =~ s/'//g;	
			my @affected_external_names = split(":",$fields[6]);
			$fields[10] =~ s/'//g;
			my @affected_relations = split(":",$fields[10]);	
			my $count = 0;
			
			for my $aff_name (@affected_external_names) {
				$aff_name =~ s/'//g;
				my @sample_objs = modules::Adaptors::Sample->search(external_sample_name=>$aff_name);
				my $aff_sample_obj;
				
				for my $tmp_obj (@sample_objs) {
					if ($tmp_obj->sample_name =~ /$fields[2]/) {
						$aff_sample_obj = $tmp_obj;
					}
				}
			
				if (!$aff_sample_obj) {
					modules::Exception->throw("ERROR: no sample obj for external name $fields[6]");
				}
				my $aff_sample_name = $aff_sample_obj->sample_name;
				(my $format_proj_name = $fields[1]) =~ s/'//g;
				my ($proj_obj) = modules::Adaptors::Project->search(id=>$source_obj->project_id);
				if ($proj_obj->project_name !~ $format_proj_name) {
					#update project
					my ($new_proj_obj) = modules::Adaptors::Project->search(project_name=>$format_proj_name);
					if (!$new_proj_obj) {
						modules::Exception->throw("ERROR: Can't retrieve proj_obj for $fields[1]");
					}
					&update_db('project_id',$new_proj_obj->id,$source_obj,$aff_sample_name);
				}
				my @lanes = modules::Adaptors::Lane->search(sample_id=>$aff_sample_obj->id);
				for my $lane_obj (@lanes) {
					if ($lane_obj->sequencing_centre ne $fields[4]) {
						#update sequencing centre
						&update_db('sequencing_centre',$fields[4],$lane_obj,$aff_sample_name) if $lane_obj->lane_name =~ /$aff_sample_name/;
					}
				}
				
				if ($aff_sample_obj->sequence_type ne $fields[5]) {
					#update sequence type
					&update_db('sequence_type',$fields[5],$aff_sample_obj,$aff_sample_name);
				}
				
				#Stuff specific to human_related_samples	
				my ($human_related_obj) = modules::Adaptors::Human_Related_Sample->search(sample_id=>$aff_sample_obj->id);
				if (!$human_related_obj) {
					modules::Exception->throw("ERROR: Can't retrieve human_related_sample for $source_name");
				}
				
				if ($human_related_obj->relation eq 'proband') {
					if ($fields[12] =~ /\w/ && $fields[12] ne $human_related_obj->sex) {
						#update proband sex
						&update_db('sex',$fields[12],$human_related_obj,$aff_sample_name);
					}
				}
			
				my $affected_db = $aff_sample_obj->sample_type =~ /_affected/?1:0;
				if (!$affected_db) {
					#update to affected
					&update_db('affected',1,$human_related_obj,$aff_sample_name);
					&update_db('sample_type','related_affected',$aff_sample_obj,$aff_sample_name);
				} 
				
				if ($human_related_obj->relation ne $affected_relations[$count]) {
					print "$human_related_obj->relation ne $affected_relations[$count])\n";
					&update_db('relation',$affected_relations[$count],$human_related_obj,$aff_sample_name);
				}
				
				$count++;
				
			}
			
			next;	
		}
		
		$command = join(" ",
							$add_multi_sample,
							'-project_name',
							$fields[1],
							'-external_name',
							$fields[2],
							'-read1_pattern',
							$fields[3],
							'-sequencing_centre',
							$fields[4],
							'-sequence_type',
							$fields[5],
							'-affected_sample_external_names',
							$affected_external_names,
							'-affected_readdirs',
							$affected_readdirs,
							'-affected_relations',
							$affected_relations
							);
		
						
		if (!$unaffected) {
			$command .= " -no_unaffected";
		} else {
			$command .= ' -unaffected_sample_external_names '.$unaffected_external_names. ' -unaffected_readdirs '. $unaffected_readdirs .' -unaffected_relations '. $unaffected_relations;
		}					
		if ($fields[12] =~ /\w/) {
			$command .= " -proband_sex $fields[12]";
		}
		if ($fields[0] eq 'human_related_gatk') {
                     $command .= ' -gatk';
		}
		
		if($fields[1] =~ /APOSLE/) 
		{
			my $steps_xml;
			if($fields[0] eq 'human_related_gatk')
			{
				$steps_xml = "$svndir/conf/human_related_gatk_APOSLE.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /German auto-immune/)
		{
			my $steps_xml;
			if($fields[0] eq 'human_related_gatk'){
				$steps_xml = "$svndir/conf/human_related_gatk_GERMAN.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		
		elsif($fields[1] =~ /CACPI/) #marcin
		{
			my $steps_xml;
			if($fields[0] eq 'human_related_gatk')
			{
				$steps_xml = "$svndir/conf/human_related_gatk_CACPI.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /CPIEEU/) #marcin
		{
			print "adding steps from human_related_gatk CPIEEU\n";
			my $steps_xml;
			if($fields[0] eq 'human_related_gatk')
			{
				$steps_xml = "$svndir/conf/human_related_gatk_CPIEEU.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}
		elsif($fields[1] =~ /CPINJ/) #marcin
		{
			print "adding steps from human_related_gatk CPINJ\n";
			my $steps_xml;
			if($fields[0] eq 'human_related_gatk')
			{
				$steps_xml = "$svndir/conf/human_related_gatk_CPINJ.steps.xml";
			}
			$command .= " -steps_xml $steps_xml";
		}	
		
		if(exists $fields[14])
		{
			if ($fields[14] =~ /apf_req_id/)
			{
				my $apf_req_id = (split(/:/, $fields[14]))[1];
				$command .= " -apf_req_id $apf_req_id";
			}
		}
		
		if($fields[1] =~ /APOSLE/)
		{
			my $report_xml = "$svndir/conf/report_APOSLE.xml";
			$command .= " -score -report_xml $report_xml";
		}
		elsif($fields[1] =~ /CACPI/) #marcin
		{
			my $report_xml = "$svndir/conf/report_CACPI.xml";
			$command .= " -score -report_xml $report_xml";
		}
		elsif($fields[1] =~ /CPIEEU/) #marcin
		{
			print "adding reports from report CPIEEU\n";
			my $report_xml = "$svndir/conf/report_CPIEEU.xml";
			$command .= " -score -report_xml $report_xml";
		}
		elsif($fields[1] =~ /CPINJ/) #marcin
		{
			print "adding reports from report CPINJ\n";
			my $report_xml = "$svndir/conf/report_CPINJ.xml";
			$command .= " -score -report_xml $report_xml";
		}

	}
	
	if ($OPT{no_mdss}) {
		$command .= ' -no_mdss';
	}
	
	if ($OPT{test}) {
		$command .= " -test";
	}
	
	
	if ($qsubs) {
		$command .= ' -skip_quality -only_qsubs';
		push @commands, $command;
	} elsif (!$source_obj) {
		if ($OPT{submit}) {
			$command .= " -submit";
		}
		push @add_commands, $command;
	} elsif ($debug) {
		push @commands, $command;
	}
}

#Debug just print any commands and quit
if ($debug) {
	for my $command ( @commands,@add_commands ) {
	    print "$command\n";
	}
	exit;
}

if (@add_commands) {
	print "\n\nAdd the following samples:\n\n";
	for my $add_command (@add_commands) {
		print "$add_command\n";
		if ($OPT{add}) {
			$sys_call->run($add_command);
		}
	}
}


if ($qsubs) {
	print "\n\nRun the following qsub only commands:\n\n";
	for my $command ( @commands ) {
	    print "$command\n";
	    $sys_call->run($command);
	}
}


#Validate line
sub Line_Problem {
	my ($line) = @_;
	my @fields = split(",",$line);
	if($line =~ 'apf_req_id'){
		$fields[0] .= '_apf';
	}
	my $error = 0;
	if (@fields != $ENTRY_NUMBER{$fields[0]}) {
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
	
	if ($fields[0] eq 'human_single' || $fields[0] eq 'human_single_gatk' || $fields[0] eq 'human_single_gatk_apf') {
		(my $affected = $fields[7]) =~ s/ //g;
		if ($affected ne 'affected' && $affected ne 'unaffected') {
			$error = "$affected must be affected or unaffected";
		}

	} elsif ($fields[0] eq 'human_related' || $fields[0] eq 'human_related_gatk' || $fields[0] eq 'human_related_gatk_apf') {
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

#Check whether a db value needs changing
sub update_db {
	my ($column,$new_value,$update_obj,$sample_name) = @_;
	print "Update $sample_name column $column to $new_value\n";
	$update_obj->$column($new_value);
	$update_obj->update();
}
