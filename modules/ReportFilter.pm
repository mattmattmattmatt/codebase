package modules::ReportFilter;

use strict;
use modules::Exception;
use modules::Utils;
use modules::ConfigXML;
use Data::Dumper;
use File::Basename;
my $EXON_NS_FILTER = 'filter_exon_ns';
my $EXON_FILTER = 'filter_exon';
my $DBSNP_FILTER = 'filter_dbsnp_snv';
my $DBSNP_FILTER_INDEL = 'filter_dbsnp_indel';
my $COSMIC_FILTER = 'filter_cosmic';
my $SPLICE_FILTER = 'filter_splicesite';
my $VEP_FILTER = 'filter_vep';
my $EXAC_FILTER = 'filter_exac_snv';
my $EXAC_FILTER_INDEL = 'filter_exac_indel';
my $GNOMAD_FILTER = 'filter_gnomad_snv';
my $GNOMAD_FILTER_INDEL = 'filter_gnomad_indel';
my $CLINVAR_FILTER = 'filter_clinvar_snv';
my $CLINVAR_FILTER_INDEL = 'filter_clinvar_indel';
my $REGULATORY_CUSTOM_FILTER = 'filter_regulatory_custom';
my $FILTER_MIRNA = 'filter_mirna';
my $FILTER_DGV = 'filter_dgv';

my $PUBMED_URL='http://www.ncbi.nlm.nih.gov/pubmed/';


sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;
	
    return $self;
}

#gets the relevant info from the snv_filter or variant_filter
sub get_filter_info {
	my %args = @_;

	#print Dumper \%args;

    my @required_args = (
			             -var_type, 
						 -filter_name
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $var_type = $args{-var_type}; #snv or indel or sv
    my $filter_obj;
      
    my $filter_name = $args{-filter_name}; #filter name
    
	
	my $attr = "N/A";
	my $attr_flag = 0;
	
	my $genome = defined $args{-genome}?1:0;

	if ($genome) {
		if (defined $args{-attr}) {
			$attr_flag = 1;
			$attr = $args{-attr};
		}
	} else {
		$filter_obj = $args{-filter_obj}; #FilteredSNV or FilteredVariant object
	    if ($var_type eq 'snv') {
	    	if (!$filter_obj->isa("modules::FilteredSNV")) {
	    		modules::Exception->throw("ERROR: Filter object isn't a filtered snv");
	    	}
	    } elsif ($var_type eq 'indel') {
	    	if (!$filter_obj->isa("modules::FilteredVariant")) {
	    		modules::Exception->throw("ERROR: Filter object isn't a filtered variant");
	    	}
	    } else {
	    	modules::Exception->throw("ERROR: var_type must be snv or indel");
	    }
		
		if ($var_type eq 'snv' && $filter_obj->snv_filter_hash->{$filter_name} ne 'nodata' && $filter_obj->snv_filter_hash->{$filter_name}->filtermatch == '1' && $filter_obj->snv_filter_hash->{$filter_name}->attribute) {
			$attr = $filter_obj->snv_filter_hash->{$filter_name}->attribute;
			$attr_flag = 1;
		} elsif ($var_type eq 'indel' && $filter_obj->variant_filter_hash->{$filter_name} ne 'nodata' && $filter_obj->variant_filter_hash->{$filter_name}->filtermatch == '1' && $filter_obj->variant_filter_hash->{$filter_name}->attribute) {
			$attr = $filter_obj->variant_filter_hash->{$filter_name}->attribute;
			$attr_flag = 1;
		} 
		
	}
	
	
	#All possible return values with defaults
	my $aa_change = "N/A";
	my $polyphen_prediction = "N/A";
	my $polyphen_score = "N/A";
	my $sift_prediction = "N/A";
	my $sift_score = "N/A";
	my $cadd_phred = "N/A";
	my $gmaf_1000_genomes = "N/A";
	my $clinical_significance = "N/A";
	my $protein_domains = "N/A";
	my $pubmed = "N/A";
	my $exon_intron_count = "N/A";
	my $dbsnp_var_allele_freq = "N/A";
	my $dbsnp_match = "N/A";
	my $cosmic_coord = 'N/A';
	my $exon_type = "SYN";
	my $splice_exon_type = "N/A";
	my $aa_position = "N/A";
	my $aa_length = "N/A";
	my $known_variation = "N/A";
	my $exac_str = "N/A";
	my $gnomad_str = "N/A";
	my $clinvar_str = "N/A";
	my $regulatory_custom_str = "N/A";
	my $dgv_freq = "N/A";
	my $mirna_str = "N/A";
	
	if ($filter_name eq $EXON_NS_FILTER) {
		if ($attr_flag) {
	    	if ($attr =~ /Stop/) {
	    		$exon_type = "NONSENSE";
	    	} else {
				$exon_type = "MISSENSE";
			}
			#aa change is required
			if ($attr =~ /aa_change\=([^\;]+)/) {
				$aa_change = $1;
				if ($attr =~ /combined\=([^\;]+)/) {
					$aa_change .= " (COMBINED:$1)";	
				}
			} #else {
				#modules::Exception->throw("ERROR: snv filter_exon_ns doesn't have the aa change");
			#}
				
						
			#polyphen and sift may be absent
			if ($attr =~ /poly_pred\=([^\;]+)/) {
				$polyphen_prediction = $1;
			}
			if ($attr =~ /poly_score\=([^\;]+)/) {
				$polyphen_score = $1;
			}
			if ($attr =~ /sift_pred\=([^\;]+)/) {
				$sift_prediction = $1;
			}
			if ($attr =~ /sift_score\=([^\;]+)/) {
				$sift_score = $1;
			}
			if ($attr =~ /cadd_phred\=([^\;]+)/) {
				$cadd_phred = $1;
			}
		}
		return ($exon_type,$aa_change,$polyphen_score,$polyphen_prediction,$sift_score,$sift_prediction,$cadd_phred);
	} elsif ($filter_name eq $EXON_FILTER) {
		if ($attr_flag) {
			if ($attr =~ /aa_pos\=([^\;]+)/) {
				$aa_position = $1;
			}
			if ($attr =~ /aa_len\=([^\;]+)/) {
				$aa_length = $1;
			}
		}
		return ($aa_position,$aa_length);		
	} elsif ($filter_name eq $SPLICE_FILTER) { 
		if ($attr_flag) {
			if ($attr =~ /splice_dist\=([^\;]+)/) {
				$splice_exon_type = "SPLICE (".$1.")";
			}							
		}
		return ($splice_exon_type);		
	} elsif ($filter_name eq $EXAC_FILTER || $filter_name eq $EXAC_FILTER_INDEL) { 
		if ($attr_flag) {
			if ($attr =~ /exac\=([^\;]+)/) {
				$exac_str = $1;
			}							
		}
		return ($exac_str);		
	} elsif ($filter_name eq $GNOMAD_FILTER || $filter_name eq $GNOMAD_FILTER_INDEL) { 
		if ($attr_flag) {
			if ($attr =~ /gnomad\=([^\;]+)/) {
				$gnomad_str = $1;
			}							
		}
		return ($gnomad_str);		
	} elsif ($filter_name eq $CLINVAR_FILTER || $filter_name eq $CLINVAR_FILTER_INDEL) {
		if ($attr_flag) {
			if ($attr =~ /clinvar\=(.*)/) {
				$clinvar_str = $1;
			}
		}
		return ($clinvar_str);
	} elsif ($filter_name eq $REGULATORY_CUSTOM_FILTER) {
		if ($attr_flag) {
			if ($attr =~ /regulatory_custom\=(.*)/) {
				$regulatory_custom_str = $1;
			}
		}
		return ($regulatory_custom_str);
	} elsif ($filter_name eq $FILTER_MIRNA) {
		if ($attr_flag) {
			if ($attr =~ /mirna\=(.*)/) {
				$mirna_str = $1;
			}
		}
		return ($mirna_str);
	} elsif ($filter_name eq $VEP_FILTER) {
		
		if ($attr_flag) {
			if ($attr =~ /pubmed\=([^\;]+)/) {
				$pubmed = '';
				my $pubmed_tmp = $1;
				#Join the urls
				my @pubmed_nums = split(",",$pubmed_tmp);
				for my $pubmed_num (@pubmed_nums) {
					$pubmed .= $PUBMED_URL . $pubmed_num . ' ';
				}
				$pubmed =~ s/ $//;
				
			}
			if ($attr =~ /exon\=([^\;]+)/) {
				$exon_intron_count = $1;
			}
			if ($attr =~ /gmaf\=([^\;]+)/) {
				$gmaf_1000_genomes = $1;
			}
			if ($attr =~ /clinsig\=([^\;]+)/) {
				$clinical_significance = $1;
			}
			if ($attr =~ /domain\=([^\;]+)/) {
				$protein_domains = $1;
			}
			if ($attr =~ /rs\=([^\;]+)/) {
				$known_variation = $1;
			}
			
		}
		return ($pubmed,$exon_intron_count,$gmaf_1000_genomes,$clinical_significance,$protein_domains,$known_variation);		
	} elsif ($filter_name eq $DBSNP_FILTER) {
		if ($attr_flag) {
			if ($attr =~ /rs=(rs\d+);ref_allele_freq\=([0-9\.]+);var_allele_freq\=([0-9\.]+)/) {
				$dbsnp_match = 'dbsnp:'. $1 . ':R' . $2 . ':V'. $3;
				$dbsnp_var_allele_freq = $3;
			} elsif ($attr =~ /rs=(rs\d+)/) {
				$dbsnp_match = 'dbsnp:'.$1.':No_freq';
				$dbsnp_var_allele_freq = 'No_freq';
			} elsif ($attr !~ /novel/) {
				$dbsnp_match = 'DBSNP_Diff_Allele';
			}
		}
		return ($dbsnp_match,$dbsnp_var_allele_freq);
		
	} elsif ($filter_name eq $DBSNP_FILTER_INDEL) {
		if ($attr_flag) {
			if ($attr =~ /rs=(rs\d+);ref_allele_freq\=([0-9\.]+);var_allele_freq\=([0-9\.]+)/) {
				$dbsnp_match = 'dbsnp:'. $1 . ':R' . $2 . ':V'. $3;
				$dbsnp_var_allele_freq = $3;
			} elsif ($attr =~ /rs=(rs\d+)/) {
				$dbsnp_match = 'dbsnp:'.$1.':No_freq';
				$dbsnp_var_allele_freq = 'No_freq';
			} elsif ($attr !~ /novel/) {
                                $dbsnp_match = 'DBSNP_Diff_Allele';
                        }
		}
		return ($dbsnp_match,$dbsnp_var_allele_freq);
		
	} elsif ($filter_name eq $FILTER_DGV) {
		
		if ($attr_flag) {
			if ($attr =~ /([0-9\.]+)%/) {
				my @matches = ($attr =~ /([0-9\.]+)%/g);
				my $match_count = @matches;
				my $sum = 0;
				for my $match ( @matches ) {
				    $sum += $match/100;
				}
				
				$dgv_freq = sprintf("%.3f",$sum/$match_count)
			} 
		}
		return ($dgv_freq);
		
	}elsif ($filter_name eq $COSMIC_FILTER) {
		if ($attr_flag) {
			my $cancer_type = defined $args{-cancer_type}?uc($args{-cancer_type}):'melanoma';
	    	if ($attr =~ /$cancer_type/i) {
	    		$cosmic_coord = "$cancer_type:".$attr;
	    	} else {
	    		$cosmic_coord = "OTHERCANCER:".$attr;
	    	}
		}	
		return ($cosmic_coord);
	} else {
		modules::Exception->throw("ERROR: No use case for filter $filter_name");
	}

	
}

return 1;
