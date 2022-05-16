#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use File::Basename;
use modules::Adaptors::Release_File;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"files=s",
	   	"runid=i",
	   	"release_files",
	   	"workdir=s",
	   	"source_name=s",
	   	"copybam=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{files} || !$OPT{source_name});

	   
=pod

=head1 SYNOPSIS

scp_results.pl -files files_to_copy -runid run_id -source_name source_name -copybam bam_to_mdss [options]

Required flags: -files -source_name

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

scp_results.pl -> Wrapper to copy tarballs to destination directory either on local or remote host

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

scp_results.pl 

=cut

my $cluster_config = modules::Pipeline::get_cluster_conf();

my $sys_call = modules::SystemCall->new();

my @files = ();
push @files, $OPT{files};

my $source_name = $OPT{source_name};
my $source_type = modules::Pipeline::get_source_type(-source_name=>$source_name);

if ($OPT{release_files}) {
	if (!$OPT{workdir}) {
		modules::Exception->throw("ERROR: Need workdir arg with release_files arg");
	}
	my $workdir = $OPT{workdir};
	my @release_files_obj = modules::Adaptors::Release_File->search_release_files_run($OPT{runid});
	for my $release_file_obj (@release_files_obj) {
		my $release_file_name = $release_file_obj->file_name;
		next if $release_file_name =~ /bam/; #These are already handled in copy_bam step
		if ($release_file_name =~ /vcf/) {
        	push @files, "$workdir/vcf/$release_file_name";
      	} elsif ($release_file_name =~ /summary/) {
        	push @files, "$workdir/summary/$release_file_name";
      	} else {
        	push @files, "$workdir/$release_file_name";
    	}
	}
}

my $archive_command = $cluster_config->read('common','archive','archive_command');
my $mdss = $cluster_config->read('common','mdss','mdss_flag');

my $copy_flag = my $scp_flag  = 0;
if ($archive_command eq 'cp') {
	$copy_flag = 1;
} elsif ($archive_command eq 'scp') {
	$scp_flag = 1;
}  else {
	modules::Exception->throw("ERROR: archive command must be cp or scp");
}



my $files = join(" ",@files);

if ($OPT{runid}) {
	$files =~ s/RUN/$OPT{runid}/g;
}

@files = split(" ",$files);

#Error check all the files exist
for my $file ( @files ) {
	
    if (!-e $file) {
    	modules::Exception->throw("ERROR: File $file doesn't exist");
    }
}

	
my $destdir = $cluster_config->read($source_type,'base_directories','base_results_directory') . '/'. $source_name;
mkdir $destdir if !-d $destdir;

#This directory should have been created during copy_bam step
if ($copy_flag && !-d $destdir) {
	modules::Exception->throw("Need to define the source directory for output");
} elsif  ($scp_flag && $cluster_config->empty('common','archive','file_server')) {
	modules::Exception->throw("Need to define the fileServer directory for scp'ing");
}

my $copy_command;
		
if ($copy_flag) {
	$copy_command = "cp $files $destdir";			
} elsif ($scp_flag) {
	my $user = `whoami`;
	chomp $user;		
	my $scp_server = $cluster_config->empty('common','archive','file_server');	
	$copy_command = "scp $files $user\@$scp_server:$destdir";
} 
#First we run the copy command
$sys_call->run("$copy_command");

#Then run the mdss command if required
if ($mdss) {
	my $mdss_results_dir = $cluster_config->read('common','mdss','mdss_results');
	my $mdss_mkdir_cmd = "mdss mkdir $mdss_results_dir/$source_type/$source_name";
	$sys_call->run($mdss_mkdir_cmd);
	for my $local_file (@files) {
		my $local_file_base = basename($local_file);
		$sys_call->run("mdss rm -f $mdss_results_dir/$source_type/$source_name/$local_file_base");
	}
	my $mdss_copy_cmd = "mdss put $files $mdss_results_dir/$source_type/$source_name";
	$sys_call->run($mdss_copy_cmd);  
	if ($OPT{copybam}) {
		if (!-e $OPT{copybam}) {
			modules::Exception->throw("ERROR: File $OPT{copybam} doesn't exist");
		}
		$sys_call->run("mdss put $OPT{copybam} $mdss_results_dir/$source_type/$source_name");
	}
}


#Then we remove the local tarball
for my $file (@files) {
	if ( $files =~ /.gz$/) {
		$sys_call->run("rm $files");
	}
}


