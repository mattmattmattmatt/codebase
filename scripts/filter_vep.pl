#! /usr/bin/perl -w

use strict;
use modules::VEP;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::BulkInsert;
use modules::VariantXML;
use modules::SystemCall;
use modules::Pipeline;
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
		   "overlap_outfile_snv=s",
		   "ref_file_snv=s",
		   "ref_file_indel=s",
		   "vep_filter_name=s",
		   "tmpdir=s",
		   "runid=i",
		   "writeDB=i",
		   "chr=s"
	    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{overlap_outfile_snv} || !$OPT{ref_file_snv} || !$OPT{ref_file_indel} || !$OPT{runid} || !$OPT{tmpdir});


=pod

=head1 SYNOPSIS

filter_vep.pl -ref_file_snv <overlap_exon_infile_snv> -ref_file_indel <overlap_exon_infile_indel> -overlap_outfile_snv <overlap_outfile_snv> -chr chr_to_run_on -vep_filter_name <filter_name for vep,default filter_vep> -writeDB 1|0

Required flags: -runid -overlap_outfile_snv -ref_file_snv -ref_file_indel -tmpdir

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_vep.pl -> Script to drive variant_effect_predictor.pl on all variants

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_vep.pl

=cut

# Put command line options into the right places
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $clus_conf = modules::Pipeline::get_cluster_conf();
my $runid = $OPT{runid};
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);


my $vep_file = $OPT{overlap_outfile_snv};
my $dbvar_file_snv = $OPT{ref_file_snv};
my $dbvar_file_indel = $OPT{ref_file_indel};

my $vep_executable = $pipe_conf->read($source_type,'binaries','variant_predictor','binary');
my $vep_db_dir = $clus_conf->read($source_type,'svn','conf_dir') . '/vep_index';
my $working_dir = $OPT{tmpdir};
my $vep_filter_name = defined $OPT{vep_filter_name}?$OPT{vep_filter_name}:'filter_vep';
my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;
my $organism = $pipe_conf->read($source_type,'organism');

my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);




my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	print "ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n";
	exit;
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $vep_input_snv = $working_dir . '/vep.filter_vep.in.snv';
my $vep_input_indel = $working_dir . '/vep.filter_vep.in.indel';

my $vep_output_snv = $working_dir . '/vep.filter_vep.out.snv';
my $vep_output_indel = $working_dir . '/vep.filter_vep.out.indel';
    
open(my $VEP_INPUT_SNV, ">$vep_input_snv")
    or modules::Exception->throw("Unable to open file for writing input format [$vep_input_snv]");

open(my $VEP_INPUT_INDEL, ">$vep_input_indel")
    or modules::Exception->throw("Unable to open file for writing input format [$vep_input_indel]");

#Connect to database

my $run = modules::Adaptors::Run->search('id' => $runid);

my $snv_found = my $indel_found = 0;

for my $snv_file (($dbvar_file_snv)) {

	if (!-e $snv_file) {
		modules::Exception->warning("Unable to open input file [$snv_file]");
		next;
	}
	
	$snv_found = 1;
	open(SNV, $snv_file) or modules::Exception->throw("Unable to open input file [$snv_file]");

	#Build up the string for input from the db_variants file
	while (<SNV>){
	    chomp;
	
		#1	55776	55776	SNV^^^C->T^^^Q202^^^TUM81^^^D34^^^het^^^16/34^^^MQ40^^^TBS:,TTt,,,.,T,,.,,tT,.T.Ttt.TtTT.tT,,NBS:.,.,.,.,.,..,,,,,....,,,.,.^^^41:40:41:41:40:41:41:41:41:41:34:34:41:41:40:38:34:39:41:41:41:38:32:37:41:41:36:39:39:37:28:26:34:35
	
	    my ($chr, $start, $end, $attributes) = split /\t/;
		my (undef,$base_change) = split('\^\^\^',$attributes);
		my ($ref_base,$var_base) = $base_change =~ /(.*)->(.*)/;
	
	
	    if (!defined $chr || ! defined $start || ! defined $end || ! defined $attributes){
			modules::Exception->warning("Did not parse input line successfully [$_]");
	    }
	
	    unless ($var_base =~ /[ATGC]/) {
			modules::Exception->warning("Did not parse attributes successfully [$attributes]");
	    }
	
		my $base_str = $ref_base . '/' . $var_base;
	    print $VEP_INPUT_SNV join("\t", $chr, $start, $end, $base_str, '+' ) . "\n";
	}
	close SNV;
}



for my $indel_file (($dbvar_file_indel)) {

	if (!-e $indel_file) {
		modules::Exception->warning("Unable to open input file [$indel_file]");
		next;
	}
	$indel_found = 1;
	open(INDEL, $indel_file) or modules::Exception->throw("Unable to open input file [$indel_file]");
		
	#Build up the string for input from the db_variants file
	while (<INDEL>){
	    chomp;
	
		#1	55776	55776 	DEL^^^AAC^^^Q51.5^^^REF_NORM60^^^SOMATIC^^^D31^^^other^^^1/31^^^MQ5^^^TBS:.$.-3AAC............,..,............^],^^^NBS.,.......,....,.......,.,.^^^20:2:5:4:2:5:5:37:5:5:5:5:41:41:35:5:5:5:5:39:5:40:5:5:32:25:24:24:24:24:26
	    my ($chr, $start, $end, $attributes) = split /\t/;
	    my $ref_base = my $var_base;
	    if ($attributes =~ /DEL/) {
			(undef,$ref_base) = split('\^\^\^',$attributes);
			$var_base = '-';
	    } elsif ($attributes =~ /INS/) {
			(undef,$var_base) = split('\^\^\^',$attributes);
			$ref_base = '-';
	    	#$end--;
	    	$start++;
	    } else {
	    	modules::Exception->throw("ERROR: Must be DEL or INS for line $_");
	    }
	
	
	    if (!defined $chr || ! defined $start || ! defined $end || ! defined $attributes){
			modules::Exception->warning("Did not parse input line successfully [$_]");
	    }
	
		my $base_str = $ref_base . '/' . $var_base;
	    print $VEP_INPUT_INDEL join("\t", $chr, $start, $end, $base_str, '+' ) . "\n";
	}
	close INDEL;
}


close($VEP_INPUT_SNV);
close($VEP_INPUT_INDEL);

my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');



# Run vep and parse result
if ($snv_found) {
	my %mutant_snv_data;
	my @vep_inserts_snv = ();
	# Work through results by chromosome. Write vep matches to database and intermediate files.
	open(my $VEP_SNV, ">$vep_file.snv") or modules::Exception->throw("Unable to open file for writing");
	my $vep_snv
	    = modules::VEP->new(-input_file      => $vep_input_snv,
					-executable_path => $vep_executable,
					-db_dir          => $vep_db_dir,
					-working_dir     => $working_dir,
					-source_type		  => $source_type,
					-exon => 0,
					-output_file => $vep_output_snv
	);
	
	if (!$vep_snv->run) {
		modules::Exception->throw("ERROR: VEP failed to run");
	}
	
	my $results_snv = $vep_snv->parse_result;
	my %snv_chrs;
	foreach my $result (@$results_snv){
	    $snv_chrs{$result->[0]}++;
	}
	foreach my $chr (sort keys %snv_chrs) {
	    foreach my $result (@$results_snv){
			next unless $result->[0] eq $chr;
		
			my ($snv_chr,
			    $snv_start,
			    $snv_end,
				$var_type,
				$se_type,
			    $var_base,
			 	$rs, 
			 	$gmaf, 
			 	$domain, 
			 	$pubmed, 
			 	$clinical, 
			 	$exon_str) = @$result;
		
			
			my $info_field = join('^^^',$rs,$gmaf,$domain,$pubmed,$clinical,$exon_str);
			
			my %attributes = ();
			
			if ($rs ne 'N/A') {
				$attributes{rs} = $rs; 
			} 
			if ($gmaf ne 'N/A') {
				$attributes{gmaf} = $gmaf; 
			} 
			if ($domain ne 'N/A') {
				$attributes{domain} = $domain; 
			} 
			if ($pubmed ne 'N/A') {
				$attributes{pubmed} = $pubmed; 
			} 
			if ($clinical ne 'N/A') {
				$attributes{clinsig} = $clinical; 
			} 
			if ($exon_str ne 'N/A') {
				$attributes{exon} = $exon_str; 
			} 
			
			my $attribute;
			if (keys %attributes) {
				$attribute = join (";", map {"$_=$attributes{$_}"} keys %attributes);
			} else {
				$attribute = "N/A";
			}
			
					                    
			print $VEP_SNV join("\t", 
			     			$snv_chr, 
			     			$snv_start, 
			     			$snv_end,
			     			'') 
							. join('^^^', 
			       					$info_field
			       				) . "\n";
			
			my $mutant_snv_key =  $snv_chr.":".$snv_start;
			$mutant_snv_data{$mutant_snv_key}{$snv_start}{$var_base}{snv_filter_vep_string} = $attribute;
			$mutant_snv_data{$mutant_snv_key}{$snv_start}{$var_base}{pass} = 1;
			
			
	    }
	}
	
	my $file_name = $sample_name . '_' . $runid . '.filter_vep.snv.xml';
	$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_snv_data, -chr=>'all');
	my $full_xml_file = $run_dir . '/conf/'.$file_name;
	$var_xml->split_xml(-file_name=>$full_xml_file);
	
	
	close($VEP_SNV);
	
}

if ($indel_found) {
	
	my %mutant_indel_data = ();
	
	my $vep_indel
	    = modules::VEP->new(-input_file      => $vep_input_indel,
					-executable_path => $vep_executable,
					-db_dir          => $vep_db_dir,
					-working_dir     => $working_dir,
					-source_type		  => $source_type,
					-exon => 0,
					-output_file => $vep_output_indel
					);
	
	if (!$vep_indel->run) {
		modules::Exception->throw("ERROR: VEP failed to run");
	}
	
	my $results_indel = $vep_indel->parse_result;
	
	
	# Get a list of chromosomes included in the results
	
	
	my %indel_chrs;
	foreach my $result (@$results_indel){
	    $indel_chrs{$result->[0]}++;
	}
	my @vep_inserts_indel = ();

	#Now the same for indel
	open(my $VEP_INDEL, ">$vep_file.indel") 
	    or modules::Exception->throw("Unable to open file for writing");
	
	
	
	foreach my $chr (sort keys %indel_chrs){
	
	    my %chr_indels;
	
		my $file_name = $sample_name . '_' . $runid . '.db_variants.indel.'.$chr.'.xml';
		$var_xml->load_xml_keys(-file_name=>$file_name);
		my $indels = $var_xml->search_xml_keys(-chr=>$chr,-var_type=>'indel');
		
		for my $indel_obj (@{$indels}) {
			my $var_type;
			if ($indel_obj->{var_bases} =~ /\+/) {
				$var_type = 'INS';
			} elsif ($indel_obj->{var_bases} =~ /-/) {
				$var_type = 'DEL';
			} else {
				print Dumper $indel_obj;
				modules::Exception->throw("ERROR: Can't determine if ins or del");
			}
			$chr_indels{$indel_obj->{start_coord}}{$indel_obj->{end_coord}}{$var_type} = $indel_obj->{chr};
		}
			

	    foreach my $result (@$results_indel){
			next unless $result->[0] eq $chr;
	
		
			my ($var_chr,
			    $var_start,
			    $var_end,
			    $var_type,
			    $se_type,
			    $var_base,
			 	$rs, 
			 	$gmaf, 
			 	$domain, 
			 	$pubmed, 
			 	$clinical, 
			 	$exon_str) = @$result;
		
			
			my $info_field = join('^^^',$rs,$gmaf,$domain,$pubmed,$clinical,$exon_str);
			
			my %attributes = ();
			
			if ($rs ne 'N/A') {
				$attributes{rs} = $rs; 
			} 
			if ($gmaf ne 'N/A') {
				$attributes{gmaf} = $gmaf; 
			} 
			if ($domain ne 'N/A') {
				$attributes{domain} = $domain; 
			} 
			if ($pubmed ne 'N/A') {
				$attributes{pubmed} = $pubmed; 
			} 
			if ($clinical ne 'N/A') {
				$attributes{clinsig} = $clinical; 
			} 
			if ($exon_str ne 'N/A') {
				$attributes{exon} = $exon_str; 
			} 
			
			
			my $attribute;
			if (keys %attributes) {
				$attribute = join (";", map { "$_=$attributes{$_}" } keys %attributes);
			} else {
				$attribute = "N/A";
			}
			
		    my $var_id;
		    if ($var_type eq 'INS') {
		    	$var_start--;
		    	$var_end--;
		    } 	                    
			
			$var_id = $chr_indels{$var_start}{$var_end}{$var_type};  	
			if ($var_id !~ /\w/) {
				modules::Exception->throw("ERROR: Couldn't retrieve var_id for $chr $var_start $var_end");
			}		       
					                    
			print $VEP_INDEL join("\t", 
			     			$var_chr, 
			     			$var_start, 
			     			$var_end,
			     			'') 
							. join('^^^', 
			       					$info_field
			       				) . "\n";
		    
	    	my $mutant_indel_key =  $var_chr.":".$var_start;
			$mutant_indel_data{$mutant_indel_key}{$var_end}{$var_base}{indel_filter_vep_string} = $attribute;
			$mutant_indel_data{$mutant_indel_key}{$var_end}{$var_base}{pass} = 1;
		   
	    }
	}
	
	my $file_name = $sample_name . '_' . $runid . '.filter_vep.indel.xml';
	$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_indel_data, -chr=>'all');
	my $full_xml_file = $run_dir . '/conf/'.$file_name;
	$var_xml->split_xml(-file_name=>$full_xml_file);

	close($VEP_INDEL);

}