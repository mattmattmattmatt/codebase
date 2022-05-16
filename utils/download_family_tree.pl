#! /usr/bin/perl -w

use strict;
use MIME::Base64;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;
use File::Copy;
use File::Basename;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;



GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"local_dir=s",
	   	"remote_scp=s",
	   	"apf_req_id=s",
	   	"raijin_dir=s"
);

	   
pod2usage(1) if ($OPT{help});
pod2usage(1) if ($OPT{help} || !$OPT{apf_req_id} || !$OPT{remote_scp});


=pod

=head1 SYNOPSIS

download_family_tree.pl -local_dir [path to save family tree file] -remote_scp [usr@host:path] -afp_req_id -raijin_dir [raijin path to copy family tree file]

Required flags: -apf_req_id -remote_scp

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

download_family_tree.pl -> downloads APF's family tree

=head1 DESCRIPTION

This script downloads APF's family tree for a given APF request id by wget, scp back to raijin

This script requires the sshpass utility to be installed for automated scp

This script requires an rc file with the pword on the first line encoded as base64 (encode it here: http://www.base64decode.org/)

=head1 AUTHOR

Vicky Cho

=head1 EXAMPLE

download_family_tree.pl -local_path /home/vicky/tmp
-remote_scp vxc221@raijin.nci.org.au -afp_req_id 212 -raijin_dir /g/data/u86/snv_pipeline_runs/v2.1_human_genome/apf_family_tree

=cut

my $local_save_dir = defined $OPT{local_dir}?$OPT{local_dir}:"./";
my $apf_req_id = $OPT{apf_req_id};
my $apf_key = "\"https://databases.apf.edu.au/ApfBioinformaticsRequest/webservice/requestFamilyTree?requestId=".$apf_req_id."&key=UJl6BOLPO0f72cOGGK63M14H4fD3yXYz\"";
my $outfile = $local_save_dir. "/Request" . $apf_req_id. "-FamilyTree.txt";
my $remote_scp = $OPT{remote_scp};

my $sys_call = modules::SystemCall->new();

my $wget_command = "wget $apf_key -O $outfile";
$sys_call->run($wget_command);

my $p_file = '/home/vicky/.raijinrc';
	
open(my $P_FILE, $p_file) or modules::Exception->throw("ERROR: Could not open $p_file");
my $firstLine = <$P_FILE>;
my $raijin_dir = defined $OPT{raijin_dir}?$OPT{raijin_dir}:"/g/data/u86/snv_pipeline_runs/v2.1_human_genome/apf_family_tree/";
my $scp_command = "sshpass -p '".decode_base64($firstLine)."' scp " . $outfile . " " .$remote_scp.":/" . $raijin_dir;

my $return_code = 0;

print "Copying $outfile to $remote_scp:$raijin_dir\n";
$return_code = system($scp_command);
modules::Exception->throw("ERROR: scp exited with error $return_code") if($return_code != 0);

