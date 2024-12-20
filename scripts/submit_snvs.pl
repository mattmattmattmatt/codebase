#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::Adaptors::Sample;
use modules::Adaptors::Run;
use modules::Pipeline;
use File::Basename;
use Pod::Usage;
use Cwd;
use Cwd 'abs_path';
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"sample_name=s",
	   	"runid=i",
	   	"submit"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sample_name} || !$OPT{runid});

	   
=pod

=head1 SYNOPSIS

call_snvs.pl -sample_name tumour_sample_name -runid runid -submit submit_qsubs

Required flags: -sample_name -runid

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

call_snvs.pl -> Submits all the snv calling jobs for each chromosome

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

submit_snvs.pl -sample_name

NOTE: DON'T CALL DIRECTORY; USE resume_run.pl TO DETERMINE IF SAMPLE READY FOR THIS STEP

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

my $runid = $OPT{runid};
my ($run_obj) = modules::Adaptors::Run->search('id' => $runid);

if (!defined $run_obj) {
	modules::Exception->throw("ERROR: No run db object with runid $runid");
}



#Get the xml files to get variables we need first
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $cluster_config = modules::Pipeline::get_cluster_conf();
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);

my @chrs = split(" ",$pipe_config->read($source_type,'annotation_version','chr'));

my $sample_name = $OPT{sample_name};

#No longer true
#if ($sample_name !~ /tumour/) {
#	modules::Exception->throw("ERROR: Only call snvs on tumour samples\n");
#}

my ($source_name) = modules::Pipeline::get_source_name($sample_name);

#Get the sample and confirm it exists
my ($sample) = modules::Adaptors::Sample->search('sample_name' => $sample_name);
if (!defined $sample) {
	modules::Exception->throw("Sample $sample_name not in database");
}

my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory') . '/'. $source_name .'/snvcalls';

my @qsub_commands = ();

my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);

my $gatkout = $run_dir .'/gatk/'.$sample_name;

opendir(DIR,$qsub_dir) || modules::Exception->throw("Can't open directory $qsub_dir");;
my @files = grep {/$sample_name/ && /tmp$/} readdir DIR;


#Submit the individual chr snv calling jobs; these run the samtools command and then the record the step as complete
for my $file ( @files ) { #Bundled into 10 jobs now so don't count by chromosome
    my $qsub_file_tmp = $qsub_dir.'/'.$file;
    (my $qsub_file = $qsub_file_tmp) =~ s/.tmp//;
    
    if ( !-e $qsub_file_tmp ) {
    	modules::Exception->throw("File $qsub_file_tmp doesn't exist");	
    }

    if ($qsub_file_tmp !~ /$sample_name/) {
	next; #avoid submitting all for pedigrees
    }
    
    open(TMP_QSUB,"$qsub_file_tmp") || modules::Exception->throw("ERROR: Can't open new qsub file $qsub_file");
    open(QSUB,">$qsub_file") || modules::Exception->throw("ERROR: Can't open new qsub file $qsub_file");
    
    #Need to update the qsub to have the run ids
    while (<TMP_QSUB>) {
    	if (/RUNID/) {
    		$_ =~ s/RUNID/$runid/;
    	}
		if (/GATKOUT/) {
			$_ =~ s/GATKOUT/$gatkout/g;
		}
	
		$_ =~ s/\.tmp//;
	    print QSUB $_;
    }
    close QSUB;
    close TMP_QSUB;
    #print STDERR "Removing file: $qsub_file_tmp\n" if $OPT{submit};
    system("rm $qsub_file_tmp") if $OPT{submit};
    
    
    push @qsub_commands, "qsub $qsub_file";
}
print STDERR "\nRun the following commands:\n\n" if @qsub_commands;
for my $command ( @qsub_commands ) {
   print STDERR "$command\n";
   system("$command") if $OPT{submit};
}




