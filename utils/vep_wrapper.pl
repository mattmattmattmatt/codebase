#! /usr/bin/perl -w

use strict;
use modules::Exception;
use modules::VEP;;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use Data::Dumper;
use modules::Pipeline;
use Bio::EnsEMBL::Registry;
use Env qw($ENSEMBL_REGISTRY);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"coord_str=s",
		"var_base=s",
		"ref_base=s",
		"strand=s",
		"all",
		"vep_in=s",
		"vep_bin=s",
		"parse_file=s"
    	);

	   
pod2usage(1) if ($OPT{help} || (!$OPT{coord_str} || !$OPT{var_base}) && !$OPT{vep_in});


=pod

=head1 SYNOPSIS

./vep_wrapper.pl -coord_str chr:coord -var_base variant_base -vep_bin pass_in_vep_bin -ref_base ref_base(default=use_API) -strand strand(default=+) -all run_on_all -indel run_for_indels(default=SNVs) -vep_in pass_in_batch_field -parse_file parse_file_name(don't run)

Required flags: (-coord_str && -var_base) || -vep_in

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

vep_wrapper.pl -> Get vep info from variant outside pipeline

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./sample_list.pl -sample_type normal_cell_line

=cut

my $source_type = 'human_related_gatk';
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $clus_conf = modules::Pipeline::get_cluster_conf();

my $vep_executable;
if (defined $OPT{vep_bin}) {
	$vep_executable = $OPT{vep_bin};
} else {
	$vep_executable = $pipe_conf->read($source_type,'binaries','variant_predictor','binary');
}
print "VEP $vep_executable\n";
my $vep_db_dir = $clus_conf->read($source_type,'svn','conf_dir') . '/vep_index';

my $organism = $pipe_conf->read($source_type,'organism');

my $workdir = `pwd`;
chomp $workdir;

my $exon = defined $OPT{all}?0:1;

my $vep_args;

if ($exon) {
	$vep_args = $pipe_conf->read($source_type,'binaries','variant_predictor','args_exon');
} else {
	$vep_args = $pipe_conf->read($source_type,'binaries','variant_predictor','args_all');
}

my $vep_input = 'vep.input'. $$;
open(INPUT,">$vep_input") || modules::Exception->throw("Can't open file for input");

my $coord_str;

my %ref_base_map;

if ($OPT{coord_str}) {
	$coord_str = $OPT{coord_str};
	my $start_coord = my $end_coord = my $chr;
	
	if ($coord_str =~ /([0-9XY]+):(\d+)\-(\d+)/) {
		$chr = $1;
		$start_coord = $2;
		$end_coord = $3;
	} elsif ($coord_str =~ /([0-9XY]+):(\d+)/) {
		$chr = $1;
		$start_coord = $2;
		$end_coord = $2;
	} else  {
		modules::Exception->throw("ERROR: coord_str must be chr:coord");
	}
	
	my $ref_base;
	if (defined $OPT{ref_base}) {
		$ref_base = $OPT{ref_base};
	} else {
		my $registry = 'Bio::EnsEMBL::Registry';
		$registry->load_all($ENSEMBL_REGISTRY);
		my $slice_adaptor = $registry->get_adaptor($organism, 'core', 'Slice');
		$ref_base = $slice_adaptor->fetch_by_region('chromosome',$chr, $start_coord, $end_coord)->seq();
	}
	$ref_base_map{"$chr:$start_coord-$end_coord"} = $ref_base;
	my $var_base = $OPT{var_base};
	my $strand = defined $OPT{strand}?$OPT{strand}:'+';
	
	print INPUT join("\t",
					$chr,
					$start_coord,
					$end_coord,
					$ref_base .'/'.$var_base,
					$strand
					) . "\n";
	
	close INPUT;
	
} else {
	#my $registry = 'Bio::EnsEMBL::Registry';
	#$registry->load_all($ENSEMBL_REGISTRY);
	#my $slice_adaptor = $registry->get_adaptor($organism, 'core', 'Slice');
	open(FILE,"$OPT{vep_in}") || modules::Exception->throw("Can't open file to write $OPT{vep_in}\n");
	while (<FILE>) {
	    chomp;
	    next unless $_ =~ /^[0-9XY]/;
	    my ($chr,$start_coord,$end_coord,$ref_base,$var_base) = split;
    	#my $ref_base = $slice_adaptor->fetch_by_region('chromosome',$chr, $start_coord, $end_coord)->seq();
		my $strand = defined $OPT{strand}?$OPT{strand}:'+';
		$ref_base_map{"$chr:$start_coord-$end_coord"} = $ref_base;
		print INPUT join("\t",
					$chr,
					$start_coord,
					$end_coord,
					$ref_base .'/'.$var_base,
					$strand
					) . "\n";
	    
	}
}
	
my $outfile;

if ($OPT{parse_file}) {
	$outfile = $OPT{parse_file};
	if ( !-e $outfile ) {
		modules::Exception->throw("File $outfile doesn't exist");
	}
} else {
	$outfile = 'vep.out'.$$;
}


my $vep = modules::VEP->new(-input_file      => $vep_input,
							-executable_path => $vep_executable,
	  						-db_dir          => $vep_db_dir,
							-working_dir     => $workdir,
							-exon => $exon,
							-source_type => $source_type,
							-output_file => $outfile
							);
$vep->run unless $OPT{parse_file};
	
my $results = $vep->parse_result;

foreach my $result (@$results){
	if ($exon) {
		my ($snv_chr,
		    $snv_start,
		    $snv_end,
		    $var_type,
			$se_type,
		    $var_base,
		    $snv_gene,
		    $snv_exon,
		    $snv_ref_aa,
		    $snv_var_aa,
		    $snv_aa_pos,
		    $poly_predict,
		    $poly_score,
		    $sift_predict,
		    $sift_score,
	    	    $cadd_phred) = @$result;


		my $ref_base = $ref_base_map{"$snv_chr:$snv_start-$snv_end"};
		my $aa_string = $snv_ref_aa.$snv_aa_pos.$snv_var_aa;
		print join("\t",
					$snv_chr,
					$snv_start,
					$snv_end,
					$ref_base,
					$var_base,
					$aa_string,
					$snv_gene,
					$snv_exon,
					$poly_predict,
	    			$poly_score,
	    			$sift_predict,
	    			$sift_score,
				$cadd_phred
					) ."\n";
	    
	}  else {
		my ($snv_chr,
		    $snv_start,
		    $snv_end,
			$var_type,
			$vep_category,
		    $var_base,
		 	$rs, 
		 	$gmaf, 
		 	$domain, 
		 	$pubmed, 
		 	$clinical, 
		 	$exon_str,
		 	$snv_gene,
		 	$snv_trans,
		 	$cadd_phred,
	    	    $gnomad_af,
	    	    $gene_name) = @$result;
		 my $ref_base = $ref_base_map{"$snv_chr:$snv_start-$snv_end"};
		 print join("\t",
					$snv_chr,
					$snv_start,
					$snv_end,
					$ref_base,
					$var_base,
					$rs, 
		 			$gmaf, 
		 			$domain, 
		 			$pubmed, 
		 			$clinical, 
		 			$exon_str,
		 			$snv_gene,
		 			$snv_trans,
					$cadd_phred,
					$vep_category,
	    	    $gnomad_af,
	    	    $gene_name
					) ."\n";
		
		
	}
	
		
	
} 










