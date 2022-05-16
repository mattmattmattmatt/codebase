#! /usr/bin/perl -w

use strict;
use MIME::Base64;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;
use File::Copy;
use File::Basename;
use modules::Adaptors::Release_File;
use modules::Adaptors::Run;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;


GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"local_path=s",
	   	"remote_scp=s",
	   	"bam",
	   	"run"
);

	   
pod2usage(1) if ($OPT{help});


=pod

=head1 SYNOPSIS

copy_results_to_local.pl -tmp_dir [dir] -local_path [path] -remote_scp [usr@host:path] -bam -run

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

copy_results_to_local.pl -> copy production run results from raijin to local dir

=head1 DESCRIPTION

This script queries the DB for samples information, finds and downloads result or bam files from raijin.

-bam option downloads bam files only

-run option to execute commands 

This script requires the sshpass utility to be installed for automated scp

This script requires an rc file with the pword on the first line encoded as base64 (encode it here: http://www.base64decode.org/)

=head1 AUTHOR

Vicky Cho

=head1 EXAMPLE

copy_results_to_local.pl -local_path /home/vicky/tmp_results/
-remote_scp vxc221@raijin.nci.org.au -run

copy_results_to_local.pl -local_path /home/vicky/tmp_results/
-remote_scp vxc221@raijin.nci.org.au -bam -run

=cut

my $cluster_config = modules::Pipeline::get_cluster_conf();

my $sys_call = modules::SystemCall->new();

my $run = defined $OPT{run}?1:0;

if (!defined $OPT{local_path} || !-d $OPT{local_path}) {
	modules::Exception->throw("ERROR: Cannot read local_path param $OPT{local_path}");
} elsif(substr($OPT{local_path},-1,1) ne '/' ){
	modules::Exception->throw("ERROR: local_path param $OPT{local_path} must end with a '/'");
} else {
	print "local_path: $OPT{local_path}\n";
}

my $source_names = _prompt_option(-name => 'source_name ranges', -format => "[eg. APOSLE_cohort1-APOSLE_cohort10,APOSLE_single1]");
my @source_name_list = _parse_id_list($source_names);
my @sample_name_list;

my $excl_result_file_regex = "filterList";

my %sample_info;

for my $source_name (@source_name_list){
	
	#Get sample_name list for source_name
	my ($source_obj) = modules::Adaptors::Source->search('external_source_name' => $source_name);
	if(!defined $source_obj){
		modules::Exception->throw("Source $source_name not in database");
	}
	my $source_type = $source_obj->source_type;
	my ($source_group_obj) = modules::Adaptors::Source_Group->search('source_id' => $source_obj->id);
	my @samples = modules::Adaptors::Sample->search('source_group_id'=>$source_group_obj->id);
    
    my $rundir_base = $cluster_config->read($source_type,'base_directories','base_run_directory');
    my $resultdir = $cluster_config->read($source_type,'base_directories','base_results_directory');
    
    my $local_source_dir = $OPT{local_path}.$source_name;
    my $mkdir_command = "mkdir $local_source_dir";
    print $mkdir_command."\n";
   	if($run){$sys_call->run($mkdir_command);}
   
   	if($OPT{bam}){
   		my $local_bam_dir = $local_source_dir."/bam";
   		my $mkdir_command = "mkdir $local_bam_dir";
   		print $mkdir_command."\n";
   		if($run){$sys_call->run($mkdir_command);}
   		
   		my $bam_files = $resultdir."/".$source_name."/"."*bam*";
       	
       	if($run){
       		_scp_from_raijin(-download_files=>$bam_files, -dir=>$local_bam_dir, -run=>'1');
       	}else{
       		_scp_from_raijin(-download_files=>$bam_files, -dir=>$local_bam_dir, -run=>'0');
       	}
       	   		
   		next;
   	}
   
    for my $sample_obj ( @samples ) {
    	
    	my $sample_name = $sample_obj->sample_name;
    	my $sample_id = $sample_obj->id;
    	my ($run_obj) = modules::Adaptors::Run->search_latest_run_date($sample_id);
		my $runid = $run_obj->id;
		
		$sample_info{$source_name}{$sample_name}=$sample_obj->external_sample_name;
		
		#The name of the run directory (eg A15_sg1_tumour2_12)
		my $rundir_name = $sample_name . '_' . $runid;
		
		#The full directory (eg /home/matt/runs/A15_sg1_tumour2_12)
		my $workdir = join('/',$rundir_base,$source_name,$sample_name.'_runs',$rundir_name);
		my $summary_files = $workdir."/summary/*";
       	
       	my $local_summary_dir = $local_source_dir."/";

       	if($source_type eq 'human_related'){
       		$local_summary_dir = $local_summary_dir."individual_report/";
    		my $mkdir_command = "mkdir $local_summary_dir";
    		print $mkdir_command."\n";
	  		if($run && !-d $local_summary_dir){$sys_call->run($mkdir_command);}
    	
    		$local_summary_dir = $local_summary_dir.$sample_name."/";
    		$mkdir_command = "mkdir $local_summary_dir";
    		print $mkdir_command."\n";
	   		if($run){$sys_call->run($mkdir_command);}
    	}
       	if($run){
       		_scp_from_raijin(-download_files=>$summary_files, -dir=>$local_summary_dir, -run=>'1');
       	}else{
       		_scp_from_raijin(-download_files=>$summary_files, -dir=>$local_summary_dir, -run=>'0');
       	}
       	
       	my $remove_command = "rm $local_summary_dir*$excl_result_file_regex*";
       	print $remove_command."\n";
       	if($run){$sys_call->run($remove_command);}
         	
    }
    
    # download group summary for human_related source
   	if($source_type eq 'human_related'){

   		my $local_group_summary_dir = $local_source_dir."/"."group_summary_report/";
    	my $mkdir_command = "mkdir $local_group_summary_dir";
    	print $mkdir_command."\n";
   		if($run){$sys_call->run($mkdir_command);}
   		
   		my @files_regex = ("cover","indel","snv");
   		my @group_summary_files;
   		
   		for my $file_regex (@files_regex){
   			my $file = $resultdir."/".$source_name."/source_group_summaries/*".$file_regex."*tsv";
   			push @group_summary_files, $file;
   		}

   		if($run){
       		_scp_from_raijin(-download_files=>join(",",@group_summary_files), -dir=>$local_group_summary_dir, -run=>'1');
       	}else{
       		_scp_from_raijin(-download_files=>join(",",@group_summary_files), -dir=>$local_group_summary_dir, -run=>'0');
       	}
   	}
   	
   	#create sample_info file to map sample ids to external_sample_name
   	
   	my $sample_info_file = $local_source_dir."/".$source_name."_sample_info.tsv";
#   	print Dumper %sample_info;
	open(FILE,">$sample_info_file") || modules::Exception->throw("Can't open file $sample_info_file\n") if($run);
	print "Writing $sample_info_file \n";

	print FILE "Source_name  \tIndividual_ID    \tExternal_name\n" if($run);
		
	for my $sample (keys %{$sample_info{$source_name}}){
		print FILE join("\t",$source_name,$sample,$sample_info{$source_name}{$sample})."\n" if($run);
	}
	
	close(FILE) if($run);
	
}

 
 
sub _scp_from_raijin{
	my %args = @_;
    my @required_args = (
			             -download_files, # list of files in array,
			             -dir			# local dir to save download files
						 );
	
	my $p_file = '/home/vicky/.raijinrc';
	
	open(my $P_FILE, $p_file)
	or modules::Exception->throw("ERROR: Could not open $p_file");
	my $firstLine = <$P_FILE>;
	my $download_list = $args{-download_files};
	my $down_dir = $args{-dir};
	my $download_command = "sshpass -p '".decode_base64($firstLine)."' scp ".$OPT{remote_scp}."\\{".$download_list."\\} $down_dir";
	print "Downloading: sshpass -p 'hidden' scp ".$OPT{remote_scp}."\\{".$download_list."\\} $down_dir\n";
	
	my $return_code = 0;
	if($args{-run}==1){
		$return_code = system($download_command);
		modules::Exception->throw("ERROR: scp exited with error $return_code") if($return_code != 0);
	}
}


# get individual IDs from input string
sub _parse_id_list{
	
	my ($user_ids) = @_;
	my @user_id_entries = split /[,\s]+/, $user_ids;
	my @range_split;
	my @ids = ();
	my $range_start;
	my $range_end;
	my $prefix;
	my $suffix;
	
	# for each entry separated by a comma or space
	foreach my $entry (@user_id_entries){
	
		# if it's a range
		if($entry =~ /-/){
			
			# get the range prefix, start and end
			@range_split = split('-',$entry);
			$range_split[0] =~ /(\D+)/;
			$prefix = $1;
			$range_split[0] =~ /(\d+)/;
			$range_start = $1;
			$range_split[1] =~ /(\d+)/;
			$range_end = $1;
			
			# go through the range and add to the relevant list
			while($range_start <= $range_end){

				$suffix = $range_start;

				push @ids, $prefix.$suffix;
				$range_start++;
			}
		
		# otherwise just add it straight to the relevant list	
		} else{
			push @ids, $entry;
		}
	}
	
	return @ids;
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
# get option from STDIN 
sub _get_option {
    
    my %args = @_;
    
    my $variable = $args{-default};
    my $name = $args{-name};
    my @possibilities = ();
    if (defined $args{-possibilities}) {
		@possibilities = @{$args{-possibilities}};
    }
    my $format_str = defined $args{-format}?$args{-format}:'';
    my $not_empty = defined $args{-required}?1:0;
    
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
		my $user_input = <STDIN>;
		chomp $user_input;
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
