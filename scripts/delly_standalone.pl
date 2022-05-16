#! /usr/bin/perl -w

use strict;
use modules::VariantXML;
use modules::SystemCall;
use modules::Pipeline;
use modules::Delly;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Env qw($ENSEMBL_REGISTRY);
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "bam=s",
		   "outvcf=s",
		   "outfile=s",
		   "fast"
		    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{bam} || !$OPT{outvcf} || !$OPT{outfile});


=pod

=head1 SYNOPSIS

delly_standalone.pl -bam bam_input -outvcf delly_out_vcf -outfile delly_out_file -fast run_with_faster_params

Required flags: -bam -outvcf -outfile

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

delly.pl -> wrapper to run delly and parse results for 5 sv types

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./delly.pl

=cut

# Put command line options into the right places
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $clus_conf = modules::Pipeline::get_cluster_conf();
my $source_type = 'human_single_gatk';
my $genome_fasta = $clus_conf->read($source_type,'svn','fasta_file');

my $delly_executable = $pipe_conf->read($source_type,'binaries','delly','binary');
my $delly_exclude = $pipe_conf->read($source_type,'binaries','delly','exclude');
my $bcftools_executable = $pipe_conf->read($source_type,'binaries','bcftools','binary');
my $vcf_out_base = $OPT{outvcf};
my $file_out = $OPT{outfile};
my $bam = $OPT{bam};
my @sv_types = qw(del dup inv ins tra);


my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	print "ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n";
	exit;
} else {
	$svndir = $ENV{'SVNDIR'};
}

my %sv_data = ();
open(FILE,">$file_out") || modules::Exception->throw("Can't open file to write $file_out\n");
my $breakpoint_file = $file_out . '.breakpoint';
open(BREAK,">$breakpoint_file") || modules::Exception->throw("Can't open file $breakpoint_file\n");

#Need separate file for translocations as we can't overlap
open(TRA,">$file_out.tra") || modules::Exception->throw("Can't open file to write $file_out.tra"); 

my $fast = defined $OPT{fast}?1:0;

for my $sv_type ( @sv_types ) {
	(my $bcf_out = $vcf_out_base) =~ s/.vcf/.$sv_type.bcf/;
	my $call_args = $pipe_conf->read($source_type,'binaries','delly','call','args_'.$sv_type);
    my $delly = modules::Delly->new(-input_bam      => $bam,
									-executable_path => $delly_executable,
									-delly_exclude => $delly_exclude,
									-output_bcf => $bcf_out,
									-call_args => $call_args,
									-sv_type => $sv_type,
									-ref_fasta => $genome_fasta,								
									-bcftools=>$bcftools_executable,
									-fast => $fast				
									);
									
	(my $vcf_final = $bcf_out) =~ s/bcf/vcf/;
	if (-e $vcf_final && -s $vcf_final) {
		#Needed if doesn't complete in 96 hours; make sure final vcf exists and isn't empty
		print "Skip $sv_type as file $vcf_final already exists\n";		
	} else {
		$delly->run();
	}									
	
	my $results;
	if (-e $bcf_out) {
		$results = $delly->parse_result($bcf_out);
	} else {
		modules::Exception->warning("Warning: no bcf for type $sv_type");
	}
	
	for my $coord_str1 (sort {my ($a_chr,$a_coord) = $a =~ /([0-9XYM]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9XYM]+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord }  keys %{$results->{$sv_type}}) {
		for my $coord_str2 ( keys %{$results->{$sv_type}{$coord_str1}}) {
			my ($chr1,$coord1) = split(':',$coord_str1);
			my ($chr2,$coord2) = split(':',$coord_str2);
			my $anno_str = join("^^^",'SVCALLER=delly',
										'SVTYPE='.$sv_type,
										'ID='.$results->{$sv_type}{$coord_str1}{$coord_str2}{id},
										'PE='.$results->{$sv_type}{$coord_str1}{$coord_str2}{pe},
									  'SR='.$results->{$sv_type}{$coord_str1}{$coord_str2}{sr},
									  'LEN='.$results->{$sv_type}{$coord_str1}{$coord_str2}{length},
									  'QUAL='.$results->{$sv_type}{$coord_str1}{$coord_str2}{qual}
									 );
			$anno_str .= '^^^ins='.$results->{$sv_type}{$coord_str1}{$coord_str2}{seq} if $sv_type eq 'ins';						 
			
			my $second_coord = $sv_type eq 'tra'?$coord1:$coord2;
			$chr1 =~ s/chr//;
			$chr2 =~ s/chr//;
			my $sv_key = $chr1 .':' . $coord1;
			
			
			$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{id} = $results->{$sv_type}{$coord_str1}{$coord_str2}{id};
			$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{pe} = $results->{$sv_type}{$coord_str1}{$coord_str2}{pe};
			$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{sr} = $results->{$sv_type}{$coord_str1}{$coord_str2}{sr};
			$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{len} = $results->{$sv_type}{$coord_str1}{$coord_str2}{length};
			$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{qual} = $results->{$sv_type}{$coord_str1}{$coord_str2}{qual};
			
			#Extra things for translocations
			if ($sv_type eq 'tra') {
				$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{chr2_pair} = $chr2;
				$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{coord2} = $coord2;
			} elsif ($sv_type eq 'ins') {
				$sv_data{$sv_key}{$second_coord}{delly}{$sv_type}{bases} = $results->{$sv_type}{$coord_str1}{$coord_str2}{seq}; #only event type with bases
			}
			
			
			
			
			if ($sv_type eq 'tra') {
				print TRA join ("\t",
							$chr1,
							$coord1,
							$chr2,
							$coord2,
							$anno_str
							) ."\n";
			}  else {
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


