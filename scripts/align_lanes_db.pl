#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::Adaptors::Align_Lane;
use modules::Adaptors::Lane;
use modules::Adaptors::Syscall;
use modules::Adaptors::Sample;
use modules::Pipeline;
use File::Basename;
use Pod::Usage;
use Cwd;
use Cwd 'abs_path';
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"lane_id=i",
	   	"phred_quality=s",
	   	"qsub_bwa=s",
	   	"qsub_sam=s",
	   	"qsub_bam=s",
	   	"qsub_bam_stats=s",
	   	"clean",
	   	"force_pass"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{lane_id} || !$OPT{qsub_bwa} || !$OPT{qsub_bam} || !$OPT{qsub_bam_stats} || !$OPT{qsub_sam});

	   
=pod

=head1 SYNOPSIS

align_lanes_db.pl -force_pass force_pass_regardless_of_alignment_percent -lane_id lane_id -qsub_bam bam_qsub_file -clean clean_up_existing_entries_first -qsub_bwa bwa_qsub_file -qsub_sam sam_qsub_file -phred_quality phred_encoding(default=phred33)

Required flags: -lane_id -qsub_bam -qsub_bwa -qsub_bam_stats -qsub_sam

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

align_lanes_db.pl -> Create align_lane db entries after sorted and indexed bam created

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

align_lanes_db.pl 

=cut

#Reconstruct the command for Syscall
my $command =  abs_path($0)." ";

for my $opt (keys %OPT) {
	$command .= "-$opt $OPT{$opt} ";
}

$command =~ s/ $//;

my $qsub_bwa = $OPT{qsub_bwa};
if ( !-e $qsub_bwa ) {
	modules::Exception->throw("ERROR: File $qsub_bwa doesn't exist");	
}

my $qsub_sam = $OPT{qsub_sam};
if ( !-e $qsub_sam ) {
	modules::Exception->throw("ERROR: File $qsub_sam doesn't exist");	
}

my $qsub_bam = $OPT{qsub_bam};
if ( !-e $qsub_bam ) {
	modules::Exception->throw("ERROR: File $qsub_bam doesn't exist");	
}

my $qsub_bam_stats = $OPT{qsub_bam_stats};
if ( !-e $qsub_bam_stats ) {
	modules::Exception->throw("ERROR: File $qsub_bam_stats doesn't exist");	
}

#Create align_lanes in database; check db entry doesn't exist
my $lane_id = $OPT{lane_id};
my ($lane_obj) = modules::Adaptors::Lane->search(id=>$lane_id);

my $encoding = defined $OPT{phred_quality}?$OPT{phred_quality}:'phred33';

if ($encoding ne 'phred33' && $encoding ne 'phred64') {
	modules::Exception->throw("ERROR: phred_quality must be phred33 or phred64");
}

if (!defined $lane_obj) {
	modules::Exception->throw("ERROR: Cannot retrieve lane from database with lane id $lane_id");
}

my $bam_stats_commands = modules::Pipeline::get_commands($qsub_bam_stats);

my ($read_report_file) = $bam_stats_commands->[0] =~ /output\s+(\S+)/;

if ( !-e $read_report_file ) {
	modules::Exception->throw("File $read_report_file doesn't exist");	
}

#Check a min portion of reads are aligned -> source_type from lane_id
my ($sample_obj) = modules::Adaptors::Sample->search(id=>$lane_obj->sample_id);
my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_obj->sample_name);
my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_obj->sample_name);
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $min_percent = $pipe_config->read($source_type,'min_align_percent');

open(READS,"$read_report_file") || modules::Exception->throw("Can't open file $read_report_file\n");
my $total = my $aligned = 0;
while (<READS>) {
	if (/^total: (\d+)/) {
		$total = $1;
	} elsif (/^aligned: (\d+)/) {
		$aligned = $1;
	}
}

if ($total == 0 || $aligned == 0) {
	modules::Exception->throw("ERROR: Couldn't get read values from file $read_report_file");
}

my $percent_aligned = sprintf("%.2f",$aligned/$total * 100);

if ($sequence_type ne 'targeted') { #targeted can have variable % reads aligning
	if ($percent_aligned < $min_percent && !$OPT{force_pass}) {
		modules::Exception->throw("ERROR: Lane with id $lane_id has only $percent_aligned reads and required minimum $min_percent\n");
	}
}

#Check there are no existing Align_lane
my @align_lane_objs = modules::Adaptors::Align_Lane->search(lane_id=>$lane_id);

if (@align_lane_objs) {
	if ($OPT{clean}) {
		modules::Adaptors::Align_Lane->search(lane_id=>$lane_id)->delete_all;
	} else {
		modules::Exception->throw("ERROR: Lane from database with lane id $lane_id already has Lane_Align entry");
	}
}

print "P $percent_aligned\n";

my %align_lane_info = (
						quality_encoding => $encoding,
						lane_id => $lane_id,
						read_align_percent => $percent_aligned
						);

my $align_lane_db_id = modules::Adaptors::Align_Lane->insert(\%align_lane_info);
print STDERR "Created align_lane with id $align_lane_db_id for lane $lane_id\n";

for my $bam_stats_command (@{$bam_stats_commands}) {
	my $step_name = 'bam_stats';
	
	my %syscall_info = (
						command => $bam_stats_command,
						level => 'lane',
						analysis_step_name=>$step_name,
						align_lane_id=>$align_lane_db_id
						);
	my $syscall_db_id = modules::Adaptors::Syscall->insert(\%syscall_info);
}








#Now create the bam and bwa syscalls
my $bwa_commands = modules::Pipeline::get_commands($qsub_bwa);

for my $bwa_command (@{$bwa_commands}) {
	my $step_name;
	if ($bwa_command =~ /bwa aln/) {
		$step_name = 'bwa_aln';
	} 
	
	my %syscall_info = (
						command => $bwa_command,
						level => 'lane',
						analysis_step_name=>$step_name,
						align_lane_id=>$align_lane_db_id
						);
						
	my $syscall_db_id = modules::Adaptors::Syscall->insert(\%syscall_info);
}

my $sam_commands = modules::Pipeline::get_commands($qsub_sam);

for my $sam_command (@{$sam_commands}) {
	my $step_name;
	if ($sam_command =~ /bwa sampe/) {
		$step_name = 'bwa_sampe';
	} 
	
	my %syscall_info = (
						command => $sam_command,
						level => 'lane',
						analysis_step_name=>$step_name,
						align_lane_id=>$align_lane_db_id
						);
	my $syscall_db_id = modules::Adaptors::Syscall->insert(\%syscall_info);
}

my $bam_commands = modules::Pipeline::get_commands($qsub_bam);

for my $bam_command (@{$bam_commands}) {
	my $step_name;
	if ($bam_command =~ /samtools index/) {
		$step_name = 'bam_index';
	} elsif ($bam_command =~ /samtools view/) {
		$step_name = 'bam_sort';
	} 
	
	my %syscall_info = (
						command => $bam_command,
						level => 'lane',
						analysis_step_name=>$step_name,
						align_lane_id=>$align_lane_db_id
						);
	my $syscall_db_id = modules::Adaptors::Syscall->insert(\%syscall_info);
}



#Finally the align_lanes command
my %syscall_info = (
					command => $command,
					level => 'lane',
					analysis_step_name=>'align_lanes',
					align_lane_id=>$align_lane_db_id
					);

my $syscall_db_id = modules::Adaptors::Syscall->insert(\%syscall_info);
