#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use File::Basename;
use modules::ConfigXML;
use modules::Cluster;
use modules::SystemCall;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"run",
		"single=s"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});

	   
=pod

=head1 SYNOPSIS

sync_pipes.pl -run run_all_commands -single only_run_on_single_source [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

sync_pipes.pl -> Sychronise v2.1 and v2.2

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Dan Andrews

	=head1 EXAMPLE

split_by_chr.pl -input input_file -output_prefix output_stub

=cut


my $svndir = $ENV{'SVNDIR'};
my $cluster_xml = "$svndir/conf/cluster.xml";
my $cluster_config = modules::ConfigXML->new($cluster_xml);
my $scheduler = $cluster_config->read('common','scheduler');
my $cluster_obj = modules::Cluster->new(-svndir=>$svndir,-scheduler=>$scheduler);
$cluster_obj->get_job_list();
my %run_info = ();

my $run_base = '/g/data/u86/variantdb/runs/';
my $rundir_v21_single = $run_base . 'v2.1/human_single';
my $rundir_v22_single = $run_base . 'v2.2/human_single_gatk';
my $rundir_v21_related = $run_base . 'v2.1/human_related';
my $rundir_v22_related = $run_base . 'v2.2/human_related_gatk';

my $v21_rundir = '/g/data/u86/variantdb/runs/v2.1/';
my $v22_rundir = '/g/data/u86/variantdb/runs/v2.2/';
my $v21_qsub = '/g/data/u86/variantdb/qsub/v2.1/';
my $v22_qsub = '/g/data/u86/variantdb/qsub/v2.2/';

opendir(SINGLE21,$rundir_v21_single);
opendir(SINGLE22,$rundir_v22_single);
opendir(RELATED21,$rundir_v21_related);
opendir(RELATED22,$rundir_v22_related);


my @single_21 = grep {/\w/} readdir SINGLE21;
my @single_22 = grep {/\w/} readdir SINGLE22;
my @related_21 = grep {/\w/} readdir RELATED21;
my @related_22 = grep {/\w/} readdir RELATED22;

for my $single21 (sort @single_21) {
	my $job_id = $cluster_obj->check_job_running(-sample_name=>$single21.'_sg1_humansingle1');
	my $lane_job_id = $cluster_obj->check_lanejob_running(-lane_name=>$single21.'_sg1_humansingle1_l1');
	if ($job_id || $lane_job_id) {
		#Don't worry about in process jobs
		next;
	}
	
	my $results_dir = '/g/data/u86/massdata_tmp/variantdb/results/v2.1/human_single/'.$single21;
	if (-d $results_dir) {
		my @files = glob("$results_dir/$single21*tar.gz");
		my @bams = glob("$results_dir/$single21*.bam");
		if (!@bams) {
			print Dumper \@bams;
			modules::Exception->throw("Can't find bam files v21\n");
		}
		$run_info{single}{$single21.'_sg1_humansingle1'}{v21} = $bams[0] if @files;
	}
}

for my $single22 (sort @single_22) {
	my $job_id = $cluster_obj->check_job_running(-sample_name=>$single22.'_sg1_humansingle1');
	my $lane_job_id = $cluster_obj->check_lanejob_running(-lane_name=>$single22.'_sg1_humansingle1_l1');
	if ($job_id || $lane_job_id) {
		#Don't worry about in process jobs
		next;
	}
	my $results_dir = '/g/data/u86/massdata_tmp/variantdb/results/v2.2/human_single_gatk/'.$single22;
	if (-d $results_dir) {
		my @files = glob("$results_dir/$single22*tar.gz");
		my @bams = glob("${v22_rundir}/human_single_gatk/*/${single22}_sg1_humansingle1_runs/*/bam/*.bam");
		if (!@bams) {
        	print Dumper \@bams;
            modules::Exception->throw("Can't find bam files v22 ${v22_rundir}/human_single_gatk/*/${single22}_sg1_humansingle1_runs/*/bam/\n");
        }
		$run_info{single}{$single22.'_sg1_humansingle1'}{v22} = $bams[0] if @files;
	}
}

for my $related21 (sort @related_21) {
	my $results_dir = '/g/data/u86/massdata_tmp/variantdb/results/v2.1/human_related/'.$related21;
	if (-d $results_dir) {
		my @files = glob("$results_dir/$related21*tar.gz");
		for my $file (@files) {	
			my ($short_file) = basename($file);
			(my $sample = $short_file) =~ s/_\d+.tar.gz//;
			my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample);
			my $lane_job_id = $cluster_obj->check_lanejob_running(-lane_name=>$related21.'_l1');
			if ($job_id || $lane_job_id) {
				next;
			}
			(my $bam = $file) =~ s/tar.gz/merge_bam.out.bam/; 
			if (!-e $bam) {
				modules::Exception->throw("Can't find bam file $bam");		
			}
 			$run_info{related}{$sample}{v21} = $bam;
		}
	}
}

for my $related22 (sort @related_22) {
	my $results_dir = '/g/data/u86/massdata_tmp/variantdb/results/v2.2/human_related_gatk/'.$related22;
	if (-d $results_dir) {
		my @files = glob("$results_dir/$related22*tar.gz");
		for my $file (@files) {
 			my ($short_file) = basename($file);
            (my $sample = $short_file) =~ s/_\d+.tar.gz//;
            my $job_id = $cluster_obj->check_job_running(-sample_name=>$sample);
			my $lane_job_id = $cluster_obj->check_lanejob_running(-lane_name=>$related22.'_l1');
			if ($job_id || $lane_job_id) {
				next;
			}
			my @bams = glob("${v22_rundir}/human_related_gatk/*/${sample}_runs/*/bam/*.bam");
			if (!@bams) {
        		print Dumper \@bams;
            	modules::Exception->throw("Can't find bam files v22 ${v22_rundir}/human_related_gatk/*/${sample}_runs/*/bam/\n");
        	}
			$run_info{related}{$sample}{v22} = $bams[0];
     	}
	}
}

#print Dumper \%run_info;
#exit;
my %commands = ();


for my $single (sort keys %{$run_info{single}}) {
	(my $source = $single) =~ s/_sg1_humansingle1//;
	my $cp_stats_cmd;
	my @align_lane_cmds;
	my $dummy_run_cmd;
	next if keys %{$run_info{single}{$single}} == 2;
	if (!exists $run_info{single}{$single}{v21}) {
		print "Add single $single to v2.1\n";
		$cp_stats_cmd = "cp ${v22_rundir}/human_single_gatk/${source}/${single}_lanes/*stats ${v21_rundir}/human_single/${source}/${single}_lanes/";
		@align_lane_cmds = split("\n", `grep align_lanes ${v21_qsub}/human_single/${source}/lanes/*align_lane*qsub | awk '{print \"source /g/data/u86/variantdb/v2.1/conf/export_env.txt /g/data/u86/variantdb/v2.1;\",\$0,\"-force_pass\"}'`);
		$dummy_run_cmd = `grep run_pipe ${v21_qsub}/human_single/${source}/runs/*merge_bam.qsub | awk '{print \"source /g/data/u86/variantdb/v2.1/conf/export_env.txt /g/data/u86/variantdb/v2.1;\",\$0,\"-no_run\"}'`;
	} elsif (!exists $run_info{single}{$single}{v22}) {
        print "Add single $single to v2.2\n";
		$cp_stats_cmd = "cp ${v21_rundir}/human_single/${source}/${single}_lanes/*stats ${v22_rundir}/human_single_gatk/${source}/${single}_lanes/";
		@align_lane_cmds = split("\n", `grep align_lanes ${v22_qsub}/human_single_gatk/${source}/lanes/*align_lane*qsub | awk '{print \"source /g/data/u86/variantdb/v2.2/conf/export_env.txt /g/data/u86/variantdb/v2.2;\",\$0,\"-force_pass\"}'`);
		$dummy_run_cmd = `grep run_pipe ${v22_qsub}/human_single_gatk/${source}/runs/*merge_bam.qsub | awk '{print \"source /g/data/u86/variantdb/v2.2/conf/export_env.txt /g/data/u86/variantdb/v2.2;\",\$0,\"-no_run\"}'`;
	}
	chomp $dummy_run_cmd;
	push @{$commands{$source}}, $cp_stats_cmd, @align_lane_cmds, $dummy_run_cmd ;
}

for my $related (sort keys %{$run_info{related}}) {
    (my $source = $related) =~ s/_sg1_.*$//;
	my $cp_stats_cmd;
    my @align_lane_cmds;
    my $dummy_run_cmd;
	next if keys %{$run_info{related}{$related}} == 2;
	if (!exists $run_info{related}{$related}{v21}) {
		print "Add related $related to v2.1\n";
     	$cp_stats_cmd = "cp ${v22_rundir}/human_related_gatk/${source}/${related}_lanes/${related}*stats ${v21_rundir}/human_related/${source}/${related}_lanes/";
        @align_lane_cmds = split("\n", `grep align_lanes ${v21_qsub}/human_related/${source}/lanes/${related}*align_lane*qsub | awk '{print \"source /g/data/u86/variantdb/v2.1/conf/export_env.txt /g/data/u86/variantdb/v2.1;\",\$0,\"-force_pass\"}'`);
        $dummy_run_cmd = `grep run_pipe ${v21_qsub}/human_related/${source}/runs/${related}*merge_bam.qsub | awk '{print \"source /g/data/u86/variantdb/v2.1/conf/export_env.txt /g/data/u86/variantdb/v2.1;\",\$0,\"-no_run\"}'`;
 	} elsif (!exists $run_info{related}{$related}{v22}) {
        print "Add related $related to v2.2\n";
		$cp_stats_cmd = "cp ${v21_rundir}/human_related/${source}/${related}_lanes/${related}*stats ${v22_rundir}/human_related_gatk/${source}/${related}_lanes/";
        @align_lane_cmds = split("\n", `grep align_lanes ${v22_qsub}/human_related_gatk/${source}/lanes/${related}*align_lane*qsub | awk '{print \"source /g/data/u86/variantdb/v2.2/conf/export_env.txt /g/data/u86/variantdb/v2.2;\",\$0,\"-force_pass\"}'`);
        $dummy_run_cmd = `grep run_pipe ${v22_qsub}/human_related_gatk/${source}/runs/${related}*merge_bam.qsub | awk '{print \"source /g/data/u86/variantdb/v2.2/conf/export_env.txt /g/data/u86/variantdb/v2.2;\",\$0,\"-no_run\"}'`;
   	}
	chomp $dummy_run_cmd;
    push @{$commands{$source}}, $cp_stats_cmd, @align_lane_cmds, $dummy_run_cmd;
}

my $sys_call = modules::SystemCall->new();

for my $source (sort keys %commands) {
	if ($OPT{single}) {
		next unless $OPT{single} =~ /$source/;
	}
	print "Source $source commands...\n\n";
	for my $command (@{$commands{$source}}) {
		print "$command\n";
		if ($command =~ /^source/) {
			`$command` if $OPT{run}
		} else {	
			$sys_call->run($command) if $OPT{run};
		}
	}
	print "\n";
}

my %link_commands;
for my $single (sort keys %{$run_info{single}}) {
    (my $source = $single) =~ s/_sg1_humansingle1//;
	next if keys %{$run_info{single}{$single}} == 2;
    if (!exists $run_info{single}{$single}{v21}) {
		my $bam = $run_info{single}{$single}{v22};
        my @rundirs = split("\n",`ls ${v21_rundir}/human_single/${source}/${single}_runs`);
        next if !@rundirs;
        if (@rundirs == 2) {
        	print "ERROR: Multiple run directories for ${v21_rundir}/human_single/${source}/${single}_runs\n";
            exit;
       	}
        my $rundir = $rundirs[0];
        push @{$link_commands{$source}},"ln -s $bam ${v21_rundir}/human_single/${source}/${single}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam", "ln -s ${bam}.bai ${v21_rundir}/human_single/${source}/${single}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam.bai" unless -e "${v21_rundir}/human_single/${source}/${single}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam";
	} elsif (!exists $run_info{single}{$single}{v22}) {
		my $bam = $run_info{single}{$single}{v21};
        my @rundirs = split("\n",`ls ${v22_rundir}/human_single_gatk/${source}/${single}_runs`);
        next if !@rundirs;
		if (@rundirs == 2) {
        	print "ERROR: Multiple run directories for ${v22_rundir}/human_single_gatk/${source}/${single}_runs\n";
            exit;
       	}
        my $rundir = $rundirs[0];
        push @{$link_commands{$source}},"ln -s $bam ${v22_rundir}/human_single_gatk/${source}/${single}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam", "ln -s ${bam}.bai ${v22_rundir}/human_single_gatk/${source}/${single}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam.bai" unless -e "${v22_rundir}/human_single_gatk/${source}/${single}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam";
	}
}

for my $related (sort keys %{$run_info{related}}) {
	(my $source = $related) =~ s/_sg1_.*$//;
    next if keys %{$run_info{related}{$related}} == 2;
    if (!exists $run_info{related}{$related}{v21}) {
    	my $bam = $run_info{related}{$related}{v22};
        my @rundirs = split("\n",`ls ${v21_rundir}/human_related/${source}/${related}_runs`);
        next if !@rundirs;
        if (@rundirs == 2) {
        	print "ERROR: Multiple run directories for ${v21_rundir}/human_related/${source}/${related}_runs\n";
            exit;
      	}
        my $rundir = $rundirs[0];
        push @{$link_commands{$source}},"ln -s $bam ${v21_rundir}/human_related/${source}/${related}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam", "ln -s ${bam}.bai ${v21_rundir}/human_related/${source}/${related}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam.bai" unless -e "${v21_rundir}/human_related/${source}/${related}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam";
  	} elsif (!exists $run_info{related}{$related}{v22}) {
    	my $bam = $run_info{related}{$related}{v21};
        my @rundirs = split("\n",`ls ${v22_rundir}/human_related_gatk/${source}/${related}_runs`);
        next if !@rundirs;
        if (@rundirs == 2) {
        	print "ERROR: Multiple run directories for ${v22_rundir}/human_related_gatk/${source}/${related}_runs\n";
            exit;
        }
        my $rundir = $rundirs[0];
        push @{$link_commands{$source}},"ln -s $bam ${v22_rundir}/human_related_gatk/${source}/${related}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam", "ln -s ${bam}.bai ${v22_rundir}/human_related_gatk/${source}/${related}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam.bai" unless -e "${v22_rundir}/human_related_gatk/${source}/${related}_runs/${rundir}/bam/${rundir}.merge_bam.out.bam";
   	}
}



#print Dumper \%link_commands;
	
for my $source (sort keys %link_commands) {
  	if ($OPT{single}) {
		next unless $OPT{single} =~ /$source/;
	}
  	print "Source $source commands...\n\n";
    for my $command (@{$link_commands{$source}}) {
    	print "$command\n";
        $sys_call->run($command) if $OPT{run};
  	}
    print "\n";
}










