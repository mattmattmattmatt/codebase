#! /usr/bin/perl -w

use strict;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;;
use modules::Adaptors::Filter;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::Source_Group;
use modules::VariantXML;
use modules::Overlap;
use modules::Exception;
use modules::Pipeline;
use modules::Adaptors::BulkInsert;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Env qw($ENSEMBL_REGISTRY);
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "writeDB=i",
		   "coord_file_snv=s",
		   "coord_file_indel=s",
		   "coord_file_sv=s",
		   "overlap_outfile_snv=s",
		   "overlap_outfile_indel=s",
		   "overlap_outfile_sv=s",
		   "ref_file_snv=s",
		   "ref_file_indel=s",
		   "ref_file_sv=s",
		   "coord_filtername=s",
		   "runid=i",
		   "chr=s",
		   "indel_type=s",
		   "match=s",
		   "snv",
		   "exact",
		   "indel",
		   "sv",
		   "add_attribute"
		   );
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{coord_filtername} || (!$OPT{coord_file_snv} && !$OPT{coord_file_indel} && !$OPT{coord_file_sv}) || (!$OPT{overlap_outfile_snv} && !$OPT{overlap_outfile_indel} && !$OPT{overlap_outfile_sv}) || (!$OPT{ref_file_snv} && !$OPT{ref_file_indel} && !$OPT{ref_file_sv}) );

	   
=pod

=head1 SYNOPSIS

overlap_general.pl -snv overlap_snvs -indel overlap_indels -sv overlap_svs -add_attribute add_matching_value_to_snvs_filters -exact only_match_exact_coord_matches(for_deletions) -indel_type indel_enum_to_search(default=INS,DEL) -match determine_if_match_means_pass_or_fail(options=PASS/FAIL;default=FAIL) -runid runid -overlap_outfile_snv snv_output_overlap_file -overlap_outfile_indel indel_output_overlap_file -ref_file_snv snv_input_overlap_file -ref_file_indel indel_input_overlap_file -coord_file_snv formatted_snv_file -coord_file_indel formatted_indel_file -chr chrname -writeDB write_results_to_DB(default=0)[options]

Required flags: -runid -coord_filtername -runid (-overlap_outfile_snv || -overlap_outfile_indel || -overlap_outfile_sv) (-ref_file_snv || -ref_file_indel || -ref_file_sv) (-coord_file_snv || -coord_file_indel || -coord_file_sv)

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

overlap_general.pl -> Script to run a generic overlap

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./overlap_general.pl 

=cut

if (!$OPT{snv} && !$OPT{indel} && !$OPT{sv}) {
	modules::Exception->throw("ERROR: Need to overlap either -snv or -indel or -sv");
}

my $coord_filter = $OPT{coord_filtername};

if ($OPT{indel} && $coord_filter =~ /filter_dbsnp/) {
	modules::Exception->throw("ERROR: Use overlap_indels_bytype.pl to overlap indels and dbsnp");
}

if ($OPT{indel} && $coord_filter =~ /filter_common/) {
	modules::Exception->throw("ERROR: Use overlap_indels_bytype.pl to overlap indels and dbsnp");
}

if ($OPT{indel} && $coord_filter =~ /filter_exac/) {
	modules::Exception->throw("ERROR: Use overlap_indels_bytype.pl to overlap indels and dbsnp");
}

if ($OPT{indel} && $coord_filter =~ /filter_gnomad/) {
	modules::Exception->throw("ERROR: Use overlap_indels_bytype.pl to overlap indels and dbsnp");
}

if ($OPT{indel} && $coord_filter =~ /filter_clinvar/) {
	modules::Exception->throw("ERROR: Use overlap_indels_bytype.pl to overlap indels and dbsnp");
}

#Do all the common operations first
my $runid = $OPT{runid};
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $organism = $pipe_conf->read($source_type,'organism');

#Now get the peptide info from ensembl to see if we need to combine any variants
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($ENSEMBL_REGISTRY);
my $gene_adaptor = $registry->get_adaptor($organism, "Core", "Gene");

my $chromosome = defined $OPT{chr}?$OPT{chr}:'';

my $last_chr = 0;
if (defined $OPT{chr} && $OPT{chr} eq 'Y') {
	$last_chr = 1;
}

my $match = defined $OPT{match}?$OPT{match}:'FAIL';
my $filterpass;
if ($match ne 'PASS' && $match ne 'FAIL') {
	modules::Exception->throw("Match parameter must be set to PASS or FAIL (default=FAIL)");
} else {
	$filterpass = $match eq 'FAIL'?0:1;
}

my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);


#Special handling for these two filters
my $filter_exon = 0;
if ($coord_filter eq 'filter_exon') {
	$filter_exon = 1;
}
	
	
my $filter_splice = 0;
if ($coord_filter eq 'filter_splicesite') {
	$filter_splice = 1;
}	

my $filter_dbsnp = 0;
if ($coord_filter =~ /filter_dbsnp/) {
	$filter_dbsnp = 1;
}

my $filter_exac = 0;
if ($coord_filter =~ /filter_exac/) {
	$filter_exac = 1;
}

my $filter_gnomad = 0;
if ($coord_filter =~ /filter_gnomad/) {
	$filter_gnomad = 1;
}

my $filter_clinvar = 0;
if ($coord_filter =~ /filter_clinvar/) {
	$filter_clinvar = 1;
}

my $filter_regulatory_custom = 0;
if ($coord_filter =~ /filter_regulatory_custom/) {
        $filter_regulatory_custom = 1;
}

my $filter_mirna = 0;
if ($coord_filter =~ /filter_mirna/) {
        $filter_mirna = 1;
}

my $filter_common = 0;
my $ignore_string;
if ($coord_filter =~ /filter_common/) {
	#here we get the bioreps from the database to pass as -ignore arguments to overlap
	my $source_group_name = modules::Pipeline::get_source_group_name($sample_name);
	my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$source_group_name);
	
	if (!defined $source_group_obj) {
		modules::Exception->throw("ERROR: Can't get source group object from name $source_group_name");
	}
	my @sample_objs = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);
	if (!@sample_objs) {
		modules::Exception->throw("ERROR: Can't get sample objects from $source_group_name");
	}
	
	my @sample_names = ();
	for my $sample_obj ( @sample_objs ) {
	    push @sample_names, $sample_obj->sample_name;
	}
	$ignore_string = join(",",@sample_names);
	$filter_common = 1;
}

my %mutant_snv_data;
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');

if ($OPT{snv}) {

	if (!defined $OPT{coord_file_snv} || !defined $OPT{ref_file_snv} || !defined $OPT{overlap_outfile_snv}) {
		modules::Exception->throw("ERROR: For snv need to define -coord_file_snv, -ref_file_snv, and -overlap_outfile_snv");
	}

	my $coord_file = $OPT{coord_file_snv};
	
	#print "IGNORE STRING: $ignore_string\n";
	
		my $overlap_args;
	
	my $overlap = modules::Overlap->new();
	my $outfile_match = $OPT{overlap_outfile_snv};
	my $ref_file =  $OPT{ref_file_snv};
	
	#print "REF $ref_file OUT $outfile_match\n";
	my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
	#Build up the arguments for the overlap call
	my %args = (
				-ref=>$ref_file,
				-coord=>$coord_file,
				);
	

	if ($OPT{exact}) {
		$args{-exact} = 1;
	}
	
	if ($OPT{chr}) {
		$args{-chr} = $chromosome;
		$args{-append} = $outfile_match;
	} else {
		$args{-output} = $outfile_match;
	}
	
	#Here we have to get bioreps from the database and pass them to overlap object
	if ($filter_common && length($ignore_string)) {
		$args{-ignore} = $ignore_string;
	}

	#Call for matching
	($overlap_match,$match_count) = $overlap->overlap(%args);
	
	#Now call for nonmatching; need to change some of the arguments
	(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;
	
	if ($OPT{chr}) {
		$args{-append} = $outfile_nomatch;
	} else {
		$args{-output} = $outfile_nomatch;
	}
	
	$args{-fail} = 1;
	($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);	
	
	# Fetch SNPs from database
	
	my %search_params = ('run_id' => $runid);
	
	if ($OPT{chr}){
		$search_params{'chr'} = $chromosome;	
	}
	
	my @snvs = my @sorted_snvs = ();
	my $snvs;
	if ($OPT{chr}) {
		my $file_name = $sample_name . '_' . $runid . '.db_variants.snv.'.$OPT{chr}.'.xml';
		my $full_file = $var_xml->outdir().'/'.$file_name;
		
		if (-e $full_file) {
			$var_xml->load_xml_keys(-file_name=>$file_name);
			$snvs = $var_xml->search_xml_keys(-chr=>$OPT{chr},-var_type=>'snv');
		}
	} else {
		my $file_name = $sample_name . '_' . $runid . '.db_variants.snv.xml';
		$var_xml->load_xml_keys(-file_name=>$file_name);
		$snvs = $var_xml->search_xml_keys(-var_type=>'snv');
	}
	if ($snvs && @{$snvs}) {
		@sorted_snvs = sort {$a->{start_coord} <=> $b->{start_coord}} @{$snvs};
	}
	

	my @snv_filter_inserts;
	
	foreach my $newsnv  (@sorted_snvs) {
		my %snv_filter = ();
		my $chr = $newsnv->{chr};
		my $coord = $newsnv->{start_coord};
		
		#Make sure to set the attribute if it's dbsnp as this is where we record allele frequency
		if ($filter_dbsnp) {
			$snv_filter{'attribute'} = 'N/A';
		}

		if ($filter_exon) {
			$snv_filter{'attribute'} = 'N/A';
		}
	
		
		if (exists $overlap_match->{PASS}{$chr}{$coord}) {
			#If the snv is matched
			$snv_filter{'filtermatch'} = 1;
			#Set the pass based on the filter context
			$snv_filter{'filterpass'} = $filterpass == 1?1:0;
			
			
			#Check if allele freq is available for snvs
			if ($filter_dbsnp) {		
				#Get the rs id and allele frequency for dbsnp
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	my ($ref_base,$var_base) = $basic_entry =~ /([ACGT])->([ACGT])/;
				    	if (!defined $ref_base) {
				    		print STDERR "ERROR: $basic_entry\n";
				    	}
				    	
				    	#Want to report the highest ref allele freq with multiple entries
				    	my %matches = ();
				    	my %rs = ();
				    	my $rs;
				    	my $allele_freq = my $rs_found = 0;
				    	for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
					    	if ($match_entry =~ /(rs\d+).*([ACGT])\(([0-9\.]+)\)->([ACGT])\(([0-9\.]+)\)/) {
					    		if ($ref_base eq $2 && $var_base eq $4) {
									$matches{$3} = $5;
									$rs{$3} = $1;
									$allele_freq = 1;
					    		}
					    	} elsif ($match_entry =~ /(rs\d+)/) {
					    		#If there are no allele freqs we still want the rs number
					    		$rs = $1;
					    		$rs_found = 1;
					    	} else {
					    		modules::Exception->throw("ERROR: SNV match has no rs number $chr $coord $match_entry");
					    	}
				    	}
					    if ($allele_freq) {
					    	#Get the highest ref_allele value for multiple matches; we use allele freq in determining pass
					    	my ($match_ref_allele) = reverse sort {$a<=>$b} keys %matches; 
					    	$snv_filter{'attribute'} = 'rs='.$rs{$match_ref_allele} . ';ref_allele_freq=' . $match_ref_allele . ';var_allele_freq=' . $matches{$match_ref_allele};
					    } elsif ($rs_found) {
					    	$snv_filter{'attribute'} = 'rs='.$rs;
					    }
				    }
				#}
			} elsif ($filter_exon) {
				#Get the amino acid length and relative position of the mutation for filter exon
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
					    	if ($match_entry =~ /ENS/) {
					    		my $ens_gene_name;
					    	    if ($match_entry =~ /(ENSG\d+)/) {
					    			($ens_gene_name) = $1;
					    	    } elsif ($match_entry =~ /(ENSMUSG\d+)/) {
					    	    	($ens_gene_name) = $1;
					    	    } else {
									modules::Exception->throw("ERROR: Can't get ensembl gene name from $match_entry");
								}
								my $gene = $gene_adaptor->fetch_by_stable_id($ens_gene_name);
								next unless $gene;
								my $canonical_transcript = $gene->canonical_transcript();
								my $transcript_length = 0;
								my @exon_objs = @{$canonical_transcript->get_all_Exons()};
							    foreach my $exon (@exon_objs) {
							    	if ($exon->coding_region_start($canonical_transcript)) {
							    		$transcript_length += abs( $exon->coding_region_end($canonical_transcript) - $exon->coding_region_start($canonical_transcript)) + 1;
							    	}
							    }
							    
								my $strand = $canonical_transcript->strand();
								my $transcriptMapper = $canonical_transcript->get_TranscriptMapper();
								my @pepCoords = $transcriptMapper->genomic2pep($coord,$coord, $strand);
								if (@pepCoords > 1) {
									modules::Exception->throw("ERROR: SNV at $chr $coord has more than one peptide coordinate");
								}
								my $aa_length = int($transcript_length/3) - 1; # Account for stop codon
								my $aa_pos = $pepCoords[0]->start; 
								$snv_filter{'attribute'} = 'aa_pos='.$aa_pos.';aa_len='.$aa_length;
					    	} else {
								my $snv_string = $chr.':'.$coord;
								modules::Exception->throw("ERROR: Exon entry at $snv_string with match $match_entry doesn't contain ensembl gene name");
					    	}
				    	}
				    }
				#}
			} elsif ($filter_splice) {
				#Get the splice distance for filter_splice
				my $splice_distance;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
							($splice_distance) = $match_entry =~ /(\d+)$/;
				    	}
				    }
				#}
				$snv_filter{'attribute'} = 'splice_dist='.$splice_distance; 
				
			} elsif ($filter_exac) {
				#Get the splice distance for filter_splice
				my $exac_match;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
							($exac_match = $match_entry) =~ s/^\d+://; #Get rid of leading coord used to make entry unique
				    	}
				    }
				#}
				$snv_filter{'attribute'} = 'exac='.$exac_match; 
				
			} elsif ($filter_gnomad) {
				#Get the splice distance for filter_splice
				my $gnomad_match;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
							($gnomad_match = $match_entry) =~ s/^\d+://; #Get rid of leading coord used to make entry unique
				    	}
				    }
				#}
				$snv_filter{'attribute'} = 'gnomad='.$gnomad_match; 
				
			} elsif ($filter_clinvar) {
				my $clinvar_match;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
							$clinvar_match = $match_entry;
					}
				    }
				#}
				$snv_filter{'attribute'} = 'clinvar='.$clinvar_match; 
				
			} elsif ($filter_regulatory_custom) {
				my $regulatory_custom_match;
                                #for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				        for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
				                         ($regulatory_custom_match = $match_entry) =~ s/^\d+://; #Get rid of leading coord used to make entry unique
				        }
				    }
				#}
				$snv_filter{'attribute'} = 'regulatory_custom='.$regulatory_custom_match;

			} elsif ($filter_mirna) {
				my $mirna_match;
                                #for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				        for my $match_entry (sort @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}}) {
				                         $mirna_match = $match_entry; 
				        }
				    }
				#}
				$snv_filter{'attribute'} = 'mirna='.$mirna_match;

			} elsif ($OPT{add_attribute}) {
				#Here we want to add the matching field as an attibute
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$coord}} ) {
				    for my $basic_entry ( sort keys %{$overlap_match->{PASS}{$chr}{$coord}{$coord}} ) {
				    	my @matches = @{$overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}};
						if (@matches > 1) {			
							my $snv_str = $chr." ".$coord;
							modules::Exception->throw("ERROR: Shouldn't have more than 1 match for $snv_str\n");
						} else {
							my $match = $overlap_match->{PASS}{$chr}{$coord}{$coord}{$basic_entry}->[0];
							$match =~ s/-DUP//;
							$snv_filter{'attribute'} = $match;
						}
				    }
				#}
			}
		} elsif (exists $overlap_nomatch->{FAIL}{$chr}{$coord}) {
			#If the snv didn't match
			$snv_filter{'filtermatch'} = 0;
			$snv_filter{'filterpass'} = $filterpass == 0?1:0;
			
		} else {
			my $snv_str = $chr." ".$coord;
			
			print modules::Exception->throw("ERROR: Can't find snv $snv_str in the fail or pass list\n");
		}
		
		if ($snv_filter{'filterpass'} == 1 || $filter_dbsnp) {
			my $mutant_snv_key =  $chr.":".$coord;
			my $var_base = $newsnv->{var_bases};
			my $final_name = 'snv_'.$coord_filter.'_string';
			my $attribute;
			
			if ($filter_dbsnp && $snv_filter{'attribute'} eq 'N/A') {
				$attribute = 'novel';
			} else {
				$attribute = $snv_filter{'attribute'};
			}
			
			$mutant_snv_data{$mutant_snv_key}{$coord}{$var_base}{$final_name} = $attribute;
			
			if ($snv_filter{'filterpass'} == 1) {
				$mutant_snv_data{$mutant_snv_key}{$coord}{$var_base}{pass} = 1;
			} else {
				$mutant_snv_data{$mutant_snv_key}{$coord}{$var_base}{pass} = 0;
			}
		}
		
	}
	
	my $file_name = $sample_name . '_' . $runid . '.'.$coord_filter.'.snv.xml';
	if ($OPT{chr}) {
		$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_snv_data, -chr=>$OPT{chr});
		if ($last_chr) { #Chr is Y
			my $full_xml_file = $run_dir . '/conf/'.$file_name;
			$var_xml->split_xml(-file_name=>$full_xml_file);
		}
	}  else {
		$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_snv_data, -chr=>'all');
		my $full_xml_file = $run_dir . '/conf/'.$file_name;
		$var_xml->split_xml(-file_name=>$full_xml_file);
	}
	
	print STDERR "Total match count: $match_count\n";
	print STDERR "Total nomatch count: $nomatch_count\n";

}


my %mutant_indel_data;


#mostly same as snv loop above except no dbsnp handling
if ($OPT{indel}) {	
	if (!defined $OPT{coord_file_indel} || !defined $OPT{ref_file_indel} || !defined $OPT{overlap_outfile_indel}) {
		modules::Exception->throw("ERROR: For indel need to define -coord_file_indel, -ref_file_indel, and -overlap_outfile_indel");
	}
	
	my $coord_file = $OPT{coord_file_indel};
	my $indel_type = defined $OPT{indel_type}?$OPT{indel_type}:'INS,DEL';
	
	my $overlap_args;
	
	my $overlap = modules::Overlap->new();
	my $outfile_match = $OPT{overlap_outfile_indel};
	my $ref_file =  $OPT{ref_file_indel};
	
	#print "REF $ref_file OUT $outfile_match\n";
	my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
	#Build up the arguments for the overlap call
	my %args = (
				-ref=>$ref_file,
				-coord=>$coord_file,
				);
	
	#Used for deletion matching in filter_common
	if ($OPT{exact}) {
		$args{-exact} = 1;
	}
	
	if ($OPT{chr}) {
		$args{-chr} = $chromosome;
		$args{-append} = $outfile_match
	} else {
		$args{-output} = $outfile_match;
	}
	
	#Call for matching
	($overlap_match,$match_count) = $overlap->overlap(%args);
	
	#Now call for nonmatching; need to change some of the arguments
	(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;
	
	if ($OPT{chr}) {
		$args{-append} = $outfile_nomatch;
	} else {
		$args{-output} = $outfile_nomatch;
	}
	
	$args{-fail} = 1;
	($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);
		
	# Fetch INDELS from database
	
	my %search_params = ('run_id' => $runid);
	
	if ($OPT{chr}){
		$search_params{'chr'} = $chromosome;	
	}
	
	my @db_indels;
	my @sorted_indels;

	my $indels;
	if ($OPT{chr}) {
		my $file_name = $sample_name . '_' . $runid . '.db_variants.indel.'.$OPT{chr}.'.xml';
		my $full_file = $var_xml->outdir().'/'.$file_name;
 		if (-e $full_file) {
 			$var_xml->load_xml_keys(-file_name=>$file_name);
 			$indels = $var_xml->search_xml_keys(-chr=>$OPT{chr},-var_type=>'indel');
		}
	} else {
		my $file_name = $sample_name . '_' . $runid . '.db_variants.indel.xml';
		$var_xml->load_xml_keys(-file_name=>$file_name);
		$indels = $var_xml->search_xml_keys(-var_type=>'indel');
	}
		
	if ($indels && @{$indels}) {
		@sorted_indels = sort {$a->{start_coord} <=> $b->{start_coord}} @{$indels};
	}
	
	
	my @indel_filter_inserts;
	
	foreach my $newindel  (@sorted_indels) {
		my %indel_filter = ();
		
			
		my $chr = $newindel->{chr};
		my $start_coord = $newindel->{start_coord};
		my $end_coord = $newindel->{end_coord};	
		
		if (exists $overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}) {
			
			#With other inputs we can't check the allele so go by start and end coord only
			$indel_filter{'filtermatch'} = 1;
			$indel_filter{'filterpass'} = $filterpass == 1?1:0;
			if ($filter_exon) {
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$start}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}} ) {
				    	for my $match_entry (@{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}{$basic_entry}}) {
					    	if ($match_entry =~ /ENS/) {
					    		my $ens_gene_name;
					    	    if ($match_entry =~ /(ENSG\d+)/) {
					    			($ens_gene_name) = $1;
					    	    } elsif ($match_entry =~ /(ENSMUSG\d+)/) {
					    	    	($ens_gene_name) = $1;
					    	    } else {
									modules::Exception->throw("ERROR: Can't get ensembl gene name from $match_entry");
								}
								my $gene = $gene_adaptor->fetch_by_stable_id($ens_gene_name);
								next unless $gene;
								my $canonical_transcript = $gene->canonical_transcript();
								my $transcript_length = 0;
								my @exon_objs = @{$canonical_transcript->get_all_Exons()};
							    foreach my $exon (@exon_objs) {
							    	if ($exon->coding_region_start($canonical_transcript)) {
							    		$transcript_length += abs( $exon->coding_region_end($canonical_transcript) - $exon->coding_region_start($canonical_transcript)) + 1;
							    	}
							    }
							    
								my $strand = $canonical_transcript->strand();
								my $transcriptMapper = $canonical_transcript->get_TranscriptMapper();
								my @pepCoords = $transcriptMapper->genomic2pep($start_coord,$start_coord, $strand);
								if (@pepCoords > 1) {
									modules::Exception->throw("ERROR: Indel at $chr $start_coord has more than one peptide coordinate");
								}
								my $aa_length = $transcript_length/3 - 1; # Account for stop codon
								my $aa_pos = $pepCoords[0]->start; 
								$indel_filter{'attribute'} = 'aa_pos='.$aa_pos.';aa_len='.$aa_length;
					    	} else {
					    		my $indel_string = $chr.':'.$start_coord;
					    		modules::Exception->throw("ERROR: Exon entry at $indel_string with match $match_entry doesn't contain ENSEMBL gene id");
					    	}
				    	}
				    }
				#}
			} elsif ($filter_splice) {
				my $splice_distance;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$start}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}} ) {
				    	for my $match_entry (@{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}{$basic_entry}}) {
							($splice_distance) = $match_entry =~ /(\d+)$/;
				    	}
				    	
				    }
				#}
				$indel_filter{'attribute'} = 'splice_dist='.$splice_distance; 
				
			} elsif ($filter_regulatory_custom) {
				my $regulatory_custom_match;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$start}} ) {
				    for my $basic_entry (keys %{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}}){
					for my $match_entry (@{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}{$basic_entry}}) {
						($regulatory_custom_match = $match_entry) =~ s/^\d+://; #Get rid of leading coord used to make entry unique
					}
				    }
				#}
				$indel_filter{'attribute'} = 'regulatory_custom='.$regulatory_custom_match;				

			} elsif ($filter_mirna) {
				my $mirna_match;
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$start}} ) {
				    for my $basic_entry (keys %{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}}){
					for my $match_entry (@{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}{$basic_entry}}) {
						$mirna_match = $match_entry					}
				    }
				#}
				$indel_filter{'attribute'} = 'mirna='.$mirna_match;				

			} elsif ($OPT{add_attribute}) {
				#Here we want to add the matching field as an attibute
				#for my $end ( keys %{$overlap_match->{PASS}{$chr}{$start}} ) {
				    for my $basic_entry ( keys %{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}} ) {
				    	my @matches = @{$overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}{$basic_entry}};
						if (@matches > 1) {			
							my $indel_str = $chr." ".$start_coord;
							modules::Exception->throw("ERROR: Shouldn't have more than 1 match for $indel_str\n");
						} else {
							$indel_filter{'attribute'} = $overlap_match->{PASS}{$chr}{$start_coord}{$end_coord}{$basic_entry}->[0];
						}
				    }
				#}
			}
		} elsif (exists $overlap_nomatch->{FAIL}{$chr}{$start_coord}{$end_coord}) {
			#If the indel didn't match
			$indel_filter{'filtermatch'} = 0;
			$indel_filter{'filterpass'} = $filterpass == 0?1:0;
		} else {
			#This should never happen as the indel should either match or not match
			if ($coord_filter =~ /filter_common/) {
				next;
			}
			my $indel_str = $chr." ".$start_coord;
			print  modules::Exception->throw("ERROR: Can't find indel $indel_str in the fail or pass list\n");
		}
			
		if ($indel_filter{'filterpass'} == 1) {
			my $mutant_indel_key =  $chr.":".$start_coord;
			my $var_base = $newindel->{var_bases};
			my $final_name = 'indel_'.$coord_filter.'_string';
			$mutant_indel_data{$mutant_indel_key}{$end_coord}{$var_base}{$final_name} = $indel_filter{'attribute'};
		}
		
	}
	
	my $file_name = $sample_name . '_' . $runid . '.'.$coord_filter.'.indel.xml';
	if ($OPT{chr}) {
		$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_indel_data, -chr=>$OPT{chr});
		if ($last_chr) {
			my $full_xml_file = $run_dir . '/conf/'.$file_name;
			$var_xml->split_xml(-file_name=>$full_xml_file);
		}
	}  else {
		$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_indel_data, -chr=>'all');
		my $full_xml_file = $run_dir . '/conf/'.$file_name;			
		$var_xml->split_xml(-file_name=>$full_xml_file);
	}
	
	print STDERR "Total match count: $match_count\n";
	print STDERR "Total nomatch count: $nomatch_count\n";
}


#Here we overlap to gene coords both for entire sv region and breakpoints
if ($OPT{sv}) {	
	if (!defined $OPT{coord_file_sv} || !defined $OPT{ref_file_sv} || !defined $OPT{overlap_outfile_sv}) {
		modules::Exception->throw("ERROR: For indel need to define -coord_file_sv, -ref_file_sv, and -overlap_outfile_sv");
	}
	my @sv_callers = split(",",$pipe_conf->read('common','sv_callers'));
	
	my $coord_file = $OPT{coord_file_sv};
	
	my $overlap_args;
	
	my $overlap = modules::Overlap->new();
	my $overlap_breakpoint = modules::Overlap->new();
	my $outfile_match_base = $OPT{overlap_outfile_sv};
	my $ref_file_base =  $OPT{ref_file_sv};
	
	my %sv_data_all;
	my %sv_data_breakpoint;
	
	for my $sv_caller ( @sv_callers ) {
        (my $ref_file = $ref_file_base) =~ s/SVCALLER/$sv_caller/;
        (my $outfile_match = $outfile_match_base) =~ s/SVCALLER/$sv_caller/;
	
		#print "REF $ref_file OUT $outfile_match\n";
		my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
		#Build up the arguments for the overlap call
		my %args = (
					-ref=>$ref_file,
					-coord=>$coord_file,
					);
	
		

		if ($OPT{chr}) {
			$args{-chr} = $chromosome;
			$args{-append} = $outfile_match
		} else {
			$args{-output} = $outfile_match;
		}		

		#print Dumper \%args;
		
		
		#Call for matching -> filter_sv_exon only does breakpoint so skip
		($overlap_match,$match_count) = $overlap->overlap(%args);
		
		#Now call for nonmatching; need to change some of the arguments
		(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;

		if ($OPT{chr}) {
			$args{-chr} = $chromosome;
			$args{-append} = $outfile_nomatch
		} else {
			$args{-output} = $outfile_nomatch;
		}	
		
		
		$args{-fail} = 1;
		
		($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);
		
		my @sorted_svs;
		my $file_name_sv = $sample_name . '_' . $runid . '.'.$sv_caller.'.sv.xml';
			
	
	
#		if ( !-e $file_name ) {
#			modules::Exception->throw("File $file_name doesn't exist");	
#		}
		$var_xml->load_xml_keys(-file_name=>$file_name_sv,-sv=>1);
		my $svs = $var_xml->search_xml_keys(-var_type=>'sv');

		if ($svs && @{$svs}) {
			@sorted_svs = sort {$a->{chr} cmp $b->{chr} || $a->{start_coord} <=> $b->{start_coord}} @{$svs};
		}
	
		my $filter_str = 'sv_'.$coord_filter.'_string';
	
		
		#normal case with no breakpoint checking needed
		for my $sv (@sorted_svs) {
			my $newsv_chr = my $newsv_start = my $newsv_end  = my $newsv_type = my $newsv_caller;
			$newsv_chr = $sv->{chr};
			$newsv_start = $sv->{start_coord};
			$newsv_end = $sv->{end_coord};
			$newsv_type = $sv->{sv_type};
			$newsv_type =~ s/sv//;
			$newsv_caller = $sv->{sv_caller};
			$newsv_caller =~ s/sv//;
			
			next unless $sv_caller eq $newsv_caller;			
			
			if ($OPT{chr}) {
				next unless $newsv_chr eq $OPT{chr};
			}
			
			my $sv_key = $newsv_chr . ':' . $newsv_start;
			
			#print "1) SV_TYPE $sv_key $newsv_end $sv_caller $newsv_type\n";
			
		
			if (exists $overlap_match->{PASS}{$newsv_chr} && exists $overlap_match->{PASS}{$newsv_chr}{$newsv_start} && exists $overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end}) {
				for my $sv_entry ( keys %{$overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end}} ) {
					my %matches = ();
					for my $match_entry (sort @{$overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end}{$sv_entry}}) {
						$matches{$match_entry}++;
					}
					my $match_str = join(",",keys %matches);
					$sv_data_all{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{pass} = 1;
					$sv_data_all{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{$filter_str} = $match_str;				
					
				}
			} 
			
			
			
			
		}
		
		
		if ($coord_filter eq 'filter_gene' || $coord_filter eq 'filter_sv_exon') {
			#Here we need to do breakpoint again
			for my $sv_caller ( @sv_callers ) {
		        (my $ref_file_bp = $ref_file_base) =~ s/SVCALLER/$sv_caller/;
		        (my $outfile_match_bp = $outfile_match_base) =~ s/SVCALLER/$sv_caller/;
			
				#print "REF $ref_file OUT $outfile_match\n";
				my ($overlap_match_bp,$match_count_bp,$overlap_nomatch_bp,$nomatch_count_bp);
				#Build up the arguments for the overlap call
				my %args = (
							-ref=>$ref_file_bp . '.breakpoint',
							-coord=>$coord_file,
							);
				(my $out_match_bp = $outfile_match_bp) =~ s/.match/.breakpoint.match/;

				if ($OPT{chr}) {
					$args{-append} = $out_match_bp;
					$args{-chr} = $chromosome;
				} else {
					$args{-output} = $out_match_bp;
				}
			
				
				#Call for matching
				($overlap_match_bp,$match_count_bp) = $overlap_breakpoint->overlap(%args);
				
				
				#Now call for nonmatching; need to change some of the arguments
				(my $outfile_nomatch_bp = $out_match_bp) =~ s/match/nomatch/;
				if ($OPT{chr}) {
					$args{-append} = $outfile_nomatch_bp;
					$args{-chr} = $chromosome;
				} else {
					$args{-output} = $outfile_nomatch_bp;
				}
				
				$args{-fail} = 1;
				($overlap_nomatch_bp,$nomatch_count_bp) = $overlap_breakpoint->overlap(%args);
			
				for my $sv (@sorted_svs) {
					my $newsv_chr = my $newsv_start = my $newsv_end  = my $newsv_type;
					$newsv_chr = $sv->{chr};
					$newsv_start = $sv->{start_coord};
					$newsv_end = $sv->{end_coord};
					$newsv_type = $sv->{sv_type};
					$newsv_type =~ s/sv//;
					
					my $sv_key_start = $newsv_chr . ':' . $newsv_start;
					my $sv_key_end = $newsv_chr . ':' . $newsv_end;
					
					
					#print "3) SV_TYPE $sv_key $newsv_end $sv_caller $newsv_type\n";
					
					#Check for the start first
					
					if (exists $overlap_match_bp->{PASS}{$newsv_chr}{$newsv_start}{$newsv_start}) {
						for my $sv_entry ( keys %{$overlap_match_bp->{PASS}{$newsv_chr}{$newsv_start}{$newsv_start}} ) {
							
							my %matches = ();
							for my $match_entry (sort @{$overlap_match_bp->{PASS}{$newsv_chr}{$newsv_start}{$newsv_start}{$sv_entry}}) {
								$matches{$match_entry}++;
							}
							my $match_str = join(",",keys %matches);
							$sv_data_breakpoint{$sv_key_start}{$newsv_start}{$sv_caller}{$newsv_type}{pass} = 1;
							$sv_data_breakpoint{$sv_key_start}{$newsv_start}{$sv_caller}{$newsv_type}{$filter_str} = $match_str;	
						}
					} 
					
					#Then check for the end
					
					if (exists $overlap_match_bp->{PASS}{$newsv_chr}{$newsv_end}{$newsv_end}) {
						for my $sv_entry ( keys %{$overlap_match_bp->{PASS}{$newsv_chr}{$newsv_end}{$newsv_end}} ) {
							
							my %matches = ();
							for my $match_entry (sort @{$overlap_match_bp->{PASS}{$newsv_chr}{$newsv_end}{$newsv_end}{$sv_entry}}) {
								$matches{$match_entry}++;
							}
							my $match_str = join(",",keys %matches);
							$sv_data_breakpoint{$sv_key_end}{$newsv_end}{$sv_caller}{$newsv_type}{pass} = 1;
							$sv_data_breakpoint{$sv_key_end}{$newsv_end}{$sv_caller}{$newsv_type}{$filter_str} = $match_str;	
						}
					} 
					
				}
	
		
			}
		} 
		
		my $file_name = $sample_name . '_' . $runid . '.'.$coord_filter.'.sv.xml';
		my $file_name_bp = $sample_name . '_' . $runid . '.'.$coord_filter.'.bp.sv.xml';
		
		if ($OPT{chr}) {
			$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%sv_data_all, -chr=>$OPT{chr},-sv=>1);
			$var_xml->create_var_xml(-file_name=>$file_name_bp,-data=>\%sv_data_breakpoint, -chr=>$OPT{chr},-sv=>1);
			
			if ($last_chr) {
				my $full_xml_file = $run_dir . '/conf/'.$file_name;			
				$var_xml->split_xml(-file_name=>$full_xml_file);
				my $full_xml_file_bp = $run_dir . '/conf/'.$file_name_bp;			
				$var_xml->split_xml(-file_name=>$full_xml_file_bp);
			}
			
		} else {
			$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%sv_data_all, -chr=>'all',-sv=>1);
			my $full_xml_file = $run_dir . '/conf/'.$file_name;			
			$var_xml->split_xml(-file_name=>$full_xml_file);
			
			$var_xml->create_var_xml(-file_name=>$file_name_bp,-data=>\%sv_data_breakpoint, -chr=>'all',-sv=>1);
			my $full_xml_file_bp = $run_dir . '/conf/'.$file_name_bp;			
			$var_xml->split_xml(-file_name=>$full_xml_file_bp);
			
		}
		
		
			
	}
}


