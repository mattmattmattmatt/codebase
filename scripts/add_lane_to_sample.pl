#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use modules::PipelineSample;
use modules::Adaptors::Source;
use modules::Adaptors::Source_Group;
use modules::Adaptors::Sample;
use modules::Adaptors::Lane;
use modules::Adaptors::Read_File;
use modules::Adaptors::Run;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"sample_name=s",
		"sequencing_centre=s",
		"readdir=s",
		"read1=s",
		"read2=s",
		"test"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sample_name} || !$OPT{readdir} || !$OPT{read1} || !$OPT{read2});

	   
=pod

=head1 SYNOPSIS

add_lane_to_sample.pl -sample_name sample_name -readdir read_directory -read1 read1_name -read2 read2_name [options]

Required flags: -sample_name -readdir -read1 -read2

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

add_lane_to_sample.pl -> Add lanes to existing sample

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

add_lane_to_sample.pl 

=cut


my $sample_name = $OPT{sample_name};

my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);

if (!defined $sample_obj) {
	modules::Exception->throw("ERROR: Can't find db object for sample $sample_name");
}

my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);


my $cluster_config = modules::Pipeline::get_cluster_conf();
my $threadnum = $cluster_config->read('common','qsub_vars','thread_num');

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

if ($svndir =~ /trunk/) {
	modules::Exception->throw("ERROR: Do not use trunk for svn directory") unless $OPT{test};
}

my $pipe_config = modules::Pipeline::get_pipe_conf();
my %sequencing_centres = map{ $_ => 1 } split(",",$pipe_config->read('common','sequencing_centres'));

my $sequencing_centre = defined $OPT{sequencing_centre}?$OPT{sequencing_centre}:'AGRF';
if (!exists $sequencing_centres{$sequencing_centre}) {
	my $seq_centre_str = join (',',keys %sequencing_centres);
	modules::Exception->throw("ERROR: Sequencing centre must be $seq_centre_str from pipe.xml");
}




my $readdir = $OPT{readdir};
my $read1 = $OPT{read1};
my $read2 = $OPT{read2};

if (!-d $readdir) {
	modules::Exception->throw("ERROR: readdir $readdir doesn't exist");
}

if (!-e "$readdir/$read1") {
	modules::Exception->throw("ERROR: read1 file $readdir/$read1 doesn't exist");
}

if (!-e "$readdir/$read2") {
	modules::Exception->throw("ERROR: read1 file $readdir/$read2 doesn't exist");
}

#Get the highest existing lane
my $highest_lane = 0;

my @lane_objs = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);

for my $lane_obj ( @lane_objs ) {
    if ($lane_obj->lane_number > $highest_lane) {
    	$highest_lane = $lane_obj->lane_number;
    }
}
$highest_lane++;




#Create the read_files
#Acceptable suffices
my %suffices = (
				gz=>1,
				bz2=>1
				);
				
				
#Create the lane
my $lane_name = $sample_name . '_l'.$highest_lane;
my $read1_name = $lane_name . '_r1';
my $read2_name = $lane_name . '_r2';

my $read1_compressed = 0;
my $read2_compressed = 0;				
my $read1_suffix;
my $read2_suffix;
my %commands = ();
my $read1_full_file = "$readdir/$read1";
my $read2_full_file = "$readdir/$read2";

#Create symlinks to the files if they don't exist
if (-e "$readdir/$read1_name") {
	modules::Exception->throw("ERROR: Symlinks already exist to $readdir/$read1_name");
} else {
	system("cd $readdir; ln -s $read1 $read1_name");
}

#Create symlinks to the files if they don't exist
if (-e "$readdir/$read2_name") {
	modules::Exception->throw("ERROR: Symlinks already exist to $readdir/$read2_name");
} else {
	system("cd $readdir; ln -s $read2 $read2_name");
}

for my $suffix ( keys %suffices ) {
	if ($read1 =~ /$suffix$/) {
       	$read1_compressed = 1;
       	$read1_suffix = $suffix;
       	if ($suffix eq 'bz2') {
			push @{$commands{bwa}},"$svndir/scripts/compress.pl -keep -suffix $suffix -threadNum $threadnum -files $read1_full_file";
		}
    }
            
    if ($read2 =~ /$suffix$/) {
       	$read2_compressed = 1;
       	$read2_suffix = $suffix;
       	if ($suffix eq 'bz2') {
			push @{$commands{bwa}},"$svndir/scripts/compress.pl -keep -suffix $suffix -threadNum $threadnum -files $read2_full_file";
		}
    }
}
			


my %lane_info = (
	        		lane_number=>$highest_lane,
	        		lane_name=>$lane_name,
	        		sequencing_centre=>$sequencing_centre,
	        		sample_id=>$sample_obj->id
	        	);
	        	
my $lane_db_id = modules::Adaptors::Lane->insert(\%lane_info);
	        
print STDERR "Create lane $lane_name with id $lane_db_id\n";
				
#Create the read_files
my %read1_info = (
					file_name => $read1,
					is_compressed=>$read1_compressed,
					compression_suffix=>$read1_suffix,
					read_file_number=>1,
					read_directory=>$readdir,
					read_file_name=>$read1_name,
					lane_id=>$lane_db_id
					);
					
my %read2_info = (
					file_name => $read2,
					is_compressed=>$read2_compressed,
					compression_suffix=>$read2_suffix,
					read_file_number=>2,
					read_directory=>$readdir,
					read_file_name=>$read2_name,
					lane_id=>$lane_db_id
					);        

modules::Adaptors::Read_File->insert(\%read1_info);
modules::Adaptors::Read_File->insert(\%read2_info);

print STDERR "Create two read_files $read1_name and $read2_name\n";				

#Finally update the sample to add a new lane
my $current_total_lanes = $sample_obj->total_lanes;
my $new_total_lanes = $current_total_lanes + 1;
$sample_obj->total_lanes($new_total_lanes);
$sample_obj->update();

#Create the qsubs
#By default don't check the quality (assume phred33)
my $quality = defined $OPT{skip_quality}?1:0;
my $source_name = modules::Pipeline::get_source_name($sample_name);
my $outdir = $cluster_config->read($source_type,'base_directories','base_run_directory').'/'.$source_name;
my $qsub_dir = $cluster_config->read($source_type,'base_directories','base_qsub_directory').'/'.$source_name;
my $lanes_outdir = $outdir.'/'.$sample_name.'_lanes';
my $scheduler = $cluster_config->read('common','scheduler');
my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
my $pipeline_sample = modules::PipelineSample->new(-source_name=>$source_name,-cluster_obj=>$cluster_obj,-source_type=>$source_type);
my $bam_file;	
			
if (keys %commands) {
	(undef,$bam_file) = $pipeline_sample->align_lane_qsub(-sample_name=>$sample_name,-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id,-commands=>\%commands);
} else {
	(undef,$bam_file) = $pipeline_sample->align_lane_qsub(-sample_name=>$sample_name,-skip_quality=>$quality, -outdir=>$lanes_outdir,-read1=>$read1_full_file,-read2=>$read2_full_file,-lane_name=>$lane_name,-lane_id=>$lane_db_id);
}
my $steps_xml_dir = $qsub_dir . '/runs/';
if (!-d $steps_xml_dir) {
	modules::Exception->throw("ERROR: Can't locate qsub run directory $steps_xml_dir");
}
	    
my $steps_xml_out = $steps_xml_dir.$sample_name.'.xml';
if ( !-e $steps_xml_out ) {
	modules::Exception->throw("File $steps_xml_out doesn't exist");	
}

(my $delimited_bam = $bam_file) =~ s{/}{\\/}g;
#need to add this BAM file to merge_bam xml
#my $perl_command = "/usr/bin/perl -i.bak -pe 's/bam; &samtools; index OUTBAM.merge_bam.out.bam/bam $delimited_bam; &samtools; index OUTBAM.merge_bam.out.bam/' $steps_xml_out";
my $perl_command = "/usr/bin/perl -i.bak -pe 's/bamfiles /bamfiles $delimited_bam,/' $steps_xml_out";
print "$perl_command\n";
system($perl_command);

