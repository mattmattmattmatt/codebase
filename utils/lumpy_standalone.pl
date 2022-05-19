#! /usr/bin/perl -w

use strict;
use modules::VariantXML;
use modules::SystemCall;
use modules::Pipeline;
use modules::Lumpy;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;

use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "bam=s",
		   "outvcf=s",
		   "outfile=s",
		"tmpdir=s"
	    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{bam} || !$OPT{outvcf} || !$OPT{outfile} || !$OPT{tmpdir});


=pod

=head1 SYNOPSIS

lumpy_standalone.pl -bam bam_input -outvcf lumpy_out_vcf -outfile lumpy_outfile -tmpdir lumpy_tmpdir

Required flags: -bam -outvcf -outfile -tmpdir

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

lumpy.pl -> wrapper to run lumpy and parse results for 4 sv types

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./lumpy.pl

=cut

# Put command line options into the right places
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $clus_conf = modules::Pipeline::get_cluster_conf();
my $runid = $OPT{runid};
my $source_type = 'human_single_gatk';
my @sv_types = qw(del dup inv ins tra);

my $lumpy_executable = $pipe_conf->read($source_type,'binaries','lumpy','binary');
my $lumpy_exclude = $pipe_conf->read($source_type,'binaries','lumpy','exclude');
my $vcf_out = $OPT{outvcf};
my $file_out = $OPT{outfile};
my $bam = $OPT{bam};

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	print "ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n";
	exit;
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $tmpdir = $OPT{tmpdir};

if (!-e $tmpdir) {
	mkdir($tmpdir);
}


if (-e $tmpdir) {
	print "rm -Rf $tmpdir\n";
	system("rm -Rf $tmpdir");
}

my $lumpy = modules::Lumpy->new(-input_bam      => $bam,
								-executable_path => $lumpy_executable,
								-lumpy_exclude => $lumpy_exclude,
								-output_vcf => $vcf_out,
								-tmpdir => $tmpdir
								);
								
$lumpy->run();
my $results = $lumpy->parse_result();

my %sv_data = ();
open(FILE,">$file_out") || modules::Exception->throw("Can't open file to write $file_out\n");
my $breakpoint_file = $file_out . '.breakpoint';
open(BREAK,">$breakpoint_file") || modules::Exception->throw("Can't open file $breakpoint_file\n");

#Need separate file for translocations as we can't overlap
open(TRA,">$file_out.tra") || modules::Exception->throw("Can't open file to write $file_out.tra"); 


for my $sv_type ( keys %{$results} ) {
	
	for my $coord_str1 (sort {my ($a_chr,$a_coord) = $a =~ /([0-9XYM]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9XYM]+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys %{$results->{$sv_type}}) {
		for my $coord_str2 ( keys %{$results->{$sv_type}{$coord_str1}}) {
			my ($chr1,$coord1) = split(':',$coord_str1);
			my ($chr2,$coord2) = split(':',$coord_str2);
			my $anno_str = join("^^^",'SVCALLER=lumpy',
										'SVTYPE='.$sv_type,
										'ID='.$results->{$sv_type}{$coord_str1}{$coord_str2}{id},
										'PE='.$results->{$sv_type}{$coord_str1}{$coord_str2}{pe},
									  'SR='.$results->{$sv_type}{$coord_str1}{$coord_str2}{sr},
									  'LEN='.$results->{$sv_type}{$coord_str1}{$coord_str2}{length}
									 );
			$chr1 =~ s/chr//;
			$chr2 =~ s/chr//;
			
			my $second_coord = $sv_type eq 'tra'?$coord1:$coord2;
			my $sv_key = $chr1 .':' . $coord1;
			
			$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{id} = $results->{$sv_type}{$coord_str1}{$coord_str2}{id};
			$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{pe} = $results->{$sv_type}{$coord_str1}{$coord_str2}{pe};
			$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{sr} = $results->{$sv_type}{$coord_str1}{$coord_str2}{sr};
			$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{len} = $results->{$sv_type}{$coord_str1}{$coord_str2}{length};
			
			#Extra things for translocations
			if ($sv_type eq 'tra') {
				$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{chr2_pair} = $chr2;
				$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{coord2} = $coord2;
			} elsif ($sv_type eq 'ins') {
				$sv_data{$sv_key}{$second_coord}{lumpy}{$sv_type}{bases} = $results->{$sv_type}{$coord_str1}{$coord_str2}{seq};
			}
			
			if ($sv_type eq 'tra') {
				print TRA join ("\t",
							$chr1,
							$coord1,
							$chr2,
							$coord2,
							$anno_str
							) ."\n";
			} else {
					print FILE join("\t",
							$chr1,
							$coord1,
							$coord2,
							$anno_str
							) ."\n";
			}
			
			print BREAK join("\t", 
							$chr1,
							$coord1,
							$coord1,
							$anno_str .  '_1'
								) . "\n";
			
			if ($coord1 != $coord2) {				
				print BREAK join("\t", 
							$chr2,
							$coord2,
							$coord2,
							$anno_str .  '_2'
								) . "\n";
			}
			
		}
	}
}

