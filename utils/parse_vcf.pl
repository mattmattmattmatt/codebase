#! /usr/bin/perl -w

use strict;
use modules::Exception;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Set::IntSpan;
use vars qw(%OPT);
use Statistics::Descriptive::Full;

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vcf=s",
	   "vcf_cutoff=i",
	   "variant_type=s",
	   "bed=s",
	   "out=s",
	   "keep_zyg",
	   "sample_name=s",
	   "genotype_qual=i",
	   "keep_allele_freq",
	   "mean_var_freq"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf});

=pod

=head1 SYNOPSIS

parse_vcf.pl -vcf vcf_file -vcf_cutoff vcf_cutoff_score(default=40) -sample_name only_get_vars_for_this_sample -genotype_quality min_genotype_quality_cutoff -bed only_include_vars_in_bed -variant_type report_only_this_variant_type(del, ins, or snv) -out outfile(default=vcf_name.out) -keep_zyg keep_original_zyg_info -mean_var_freq report_mean_var_freq_across_all_var_samples

Required flags: -vcf 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_vcf.pl -> standalone vcf parser

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

=cut

my $vcf = $OPT{vcf};
if ( !-e $vcf ) {
	modules::Exception->throw("File $vcf doesn't exist");	
}

my $geno_qual = defined $OPT{genotype_qual}?$OPT{genotype_qual}:0;

my $out = defined $OPT{out}?$OPT{out}:$vcf.'.txt';
my $vcf_cutoff = defined $OPT{vcf_cutoff}?$OPT{vcf_cutoff}:40;

my $variant_type = 'all';
if (defined $OPT{variant_type}) {
	if ($OPT{variant_type} ne 'SNV' && $OPT{variant_type} ne 'DEL' && $OPT{variant_type} ne 'INS') {
		modules::Exception->throw("ERROR: variant_type must be INS, DEL, or SNV");
	}
	$variant_type = $OPT{variant_type};
} 


my $unhandled_count = 0;

my %bed_regions = ();

my $sample_name = defined $OPT{sample_name}?$OPT{sample_name}:0;

my $vcf_data;
if ($sample_name) {
	$vcf_data = &parse_vcf(-sample_name=>$sample_name);
} else {
	$vcf_data = &parse_vcf();
}


my $line = 0;

if ($OPT{bed}) {
	open(BED,"$OPT{bed}") || modules::Exception->throw("Can't open file $OPT{bed}\n");
	while (<BED>) {
		if ($line%1000==0) {
			print "L $line\n";
		}
		my ($chr,$start,$end) = split;
		if (!exists $bed_regions{$chr}) {
			my $set = Set::IntSpan->new("$start-$end");
			$bed_regions{$chr} = $set;
		} else {
			$bed_regions{$chr} = $bed_regions{$chr}->union(Set::IntSpan->new("$start-$end"));
		}
		$line++;
	}
}

#print Dumper \%bed_regions;
#print Dumper $vcf_data;

open(FILE,">$out") || modules::Exception->throw("Can't open file to write $out\n");
my @vcf_order;
	
for my $var_key (@vcf_order) {	
	if ($OPT{variant_type}) {
		next unless $OPT{variant_type} eq $vcf_data->{$var_key}{type};
	}
	my ($chr,$coord,$event) = split(':',$var_key);
	my ($start,$end) = split('-',$coord);
	if ($OPT{bed}) {
		if (exists $bed_regions{$chr} && $bed_regions{$chr}->intersect(Set::IntSpan->new("$start-$end"))->empty()) {
			next;
		}
	}
	$chr =~ s/chr//;

	my $var_count = 0;
	if (exists $vcf_data->{$var_key} && defined $vcf_data->{$var_key}{var_count}) {
		$var_count = $vcf_data->{$var_key}{var_count};
	}	

	my $zyg_count = 0;
        if (exists $vcf_data->{$var_key} && defined $vcf_data->{$var_key}{zyg_count}) {
                $zyg_count = $vcf_data->{$var_key}{zyg_count};
        }


	my $vcf_str  = $vcf_data->{$var_key}{type}.';'.$var_key.';Q='.$vcf_data->{$var_key}{qual}. ';AC='.$var_count . ';ZC='.$zyg_count;
	if ($OPT{keep_allele_freq}) {
		$vcf_str .= ";ALLELE=".$vcf_data->{$var_key}{allele};
	} 
	
	if ($OPT{mean_var_freq}) {
		$vcf_str .= ";MEANAF=".$vcf_data->{$var_key}{mean_af}.";MEDAF=".$vcf_data->{$var_key}{median_af}.";VAR_READ_COUNTS=".$vcf_data->{$var_key}{var_read_counts};
	}
	
	
	if ($vcf_data->{$var_key}{type} eq 'INS') {
		$vcf_str .= ";REF=".$vcf_data->{$var_key}{ref_base};
	}
	
	
	if ($OPT{keep_zyg}) {
		my $zyg_str = join("\t",@{$vcf_data->{$var_key}{zyg}});
		print FILE join("\t", 
							$chr,
							$start,
							$end,
							$vcf_str,
							$zyg_str
							)."\n";
	} else {
		print FILE join("\t", 
							$chr,
							$start,
							$end,
							$vcf_str
							)."\n";
		
	}
}


sub parse_vcf {
       
    my %vcf_data = ();
    

    open(VCF,$vcf) || modules::Exception->throw("Can't open file $vcf\n");
    
   	#Genotype quality index
   	my $gq_index = 0;
   	my $sample_index = 0;

	#Other indices strelka needs
	my ($tir_index,$tar_index,$tor_index,$a_index,$c_index,$g_index,$t_index);
	my $strelka_indel = 0;
	my $strelka_snv = 0;
	
	
	
	



    while (<VCF>) {
    	
    	
    	
    	if ($sample_name) {
			unless (/#CHROM/) {
	    		next if /^#/;				
			}    		
    	} else {
	    	next if /^#/;
    	}
    	next unless /\w/;
    	chomp;
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,@alleles) = split;
    	my $line = $_;
   		my @fields = split;
        my @gt_fields = split(':',$gt_fields) if $gt_fields;
		my @details = split(';',$rest);	

    	if (/CHROM/ && $sample_name) {
    		#Get the sample_index
    		for ( my $count = 0 ; $count < @fields ; $count++ ) {
    		    if ($fields[$count] eq $sample_name) {
    		    	$sample_index = $count;
    		    }
    		}
    		
    		if (!$sample_index) {
    			modules::Exception->throw("ERROR: Can't find sample $sample_name in line $_");
    		}
    		next;
    	}
    	
    	#Special handling for strelka AF calculations -> only do once
    	if ($gt_fields =~ /:TIR/ && $strelka_indel == 0) {
    		for ( my $count = 0 ; $count < @gt_fields ; $count++ ) {
    		    if ($gt_fields[$count] eq 'TIR') {
    		    	$tir_index = $count;
    		    } elsif ($gt_fields[$count] eq 'TAR') {
    		    	$tar_index = $count;
    		    } elsif ($gt_fields[$count] eq 'TOR') {
    		    	$tor_index = $count;
    		    } 
    		}
    		$strelka_indel = 1;
    	}
    		
    	if ($gt_fields =~ /:AU/ && $strelka_snv == 0) {
    		for ( my $count = 0 ; $count < @gt_fields ; $count++ ) {
    			if ($gt_fields[$count] eq 'AU') {
    			    $a_index = $count;
    			} elsif ($gt_fields[$count] eq 'CU') {
    			   	$c_index = $count;
    			} elsif ($gt_fields[$count] eq 'GU') {
    			   	$g_index = $count;
    			} elsif ($gt_fields[$count] eq 'TU') {
    			   	$t_index = $count;
    			} 
    		}
    		$strelka_snv = 1;
    	}
    	#exit;
    	
    	
    	if ($chr eq 'MT') {
			$chr = 'M';
		}
    	

    	if ($qual !~ /\d/ && $qual ne '.') {
    		modules::Exception->throw("ERROR: Error with qual $qual format at line $_");
    	}

		if ($ref =~ /N/ || $var_str =~ /N/) {
			next;
		}

		if ($var_str eq '.') {
			next;
		}


		my @vars;
		if ($sample_name) {
			#Get only the relevant genotype
	    	my $var_index = 0;
			my ($zyg_str) = split(':',$fields[$sample_index]);
    		next if $zyg_str eq '0/0';
			next if $zyg_str eq '0|0';
			next if $zyg_str eq '.|.';
    		next if $zyg_str eq './.';
    		my @nums = split('/',$zyg_str);
    		#Get the index from the gt string (eg 0/1 or 0/2 etc)
    		for my $num (@nums) {
    			if ($num != 0) {
    				$var_index = $num;
    			}
    		}
    		#print "sample index $sample_index $zyg_str VI $var_index\n";
    		
    		if (!$var_index) {
    			modules::Exception->throw("Can't get var_index for $_");
    		}
    		
    		#Only add the relevant genotype
    		my @tmp_vars = split(",",$var_str);
    		
    		push @vars,$tmp_vars[$var_index-1];
    		

    		if ($geno_qual && !$gq_index) {
    			#only do this once with the first non-header line; get the index of the GQ field
    			for ( my $count = 0 ; $count < @gt_fields ; $count++ ) {
    			    if ($gt_fields[$count] eq 'GQ') {
    			    	$gq_index = $count;
    			    }
    			}
    			if (!$gq_index) {
    				modules::Exception->throw("ERROR: Can't get GQ index for $gt_fields");
    			}
    		}
		} else {
			@vars = split(",",$var_str);
		}

				

		my @ac_fields =  ();
		my $total_ac = 0;

		for my $detail (@details) {
			if ($detail =~ /^AC=/) {
				$detail =~ s/AC=//;
				@ac_fields = split(",",$detail);
				for my $count (@ac_fields) {
					$total_ac += $count;
				} 
			}
			
			if ($qual eq '.' && $detail =~ /SomaticEVS=([0-9\.]+)/) {
				$qual = $1;
			}
			
			if ($qual eq '.' && $detail =~ /TLOD=([0-9\.]+)/) {
				$qual = $1;
			}
			
		}

		my $zyg_num = 1;	
		for my $var ( @vars ) {
			next if $var eq '*'; #Due to upstream deletion
			next if $var eq $ref; #Happens with Tapestri occasionally
			my ($var_key,$var_type,$ref_base) = _get_variant_key(-type=>'vcf',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var);
			next unless $var_key =~ /(\d+)-(\d+):.+/;
			my ($start,$end) = $var_key =~ /(\d+)-(\d+)/;
			#print "Here $var_key\n";

			#If fails quality test; only apply if quality score available
#			if ($qual ne '.' && $qual <= $vcf_cutoff) {
#				next;
#			}

			if ($sample_name && $geno_qual) {
				#Here we additionally screen on GQ
				my @gt_fields_data = split(':',$fields[$sample_index]);
				my $gq = $gt_fields_data[$gq_index];
				
				if ($gq < $geno_qual) {
					#Make sure it's not a het or hom decision -> this is still a variant
					my @nums = split('/',$gt_fields_data[0]);
					if ($nums[0] eq $nums[1]) {
						next;
					}
				}
				
			}

			$vcf_data{$var_key}{var_count} = shift @ac_fields;
			$vcf_data{$var_key}{zyg_count} = $zyg_num;
			
			if ($qual eq '.') {
				$vcf_data{$var_key}{qual} = "N/A";				
			} else {
				$vcf_data{$var_key}{qual} = $qual;
			}
			$vcf_data{$var_key}{type} = $var_type;
			$vcf_data{$var_key}{ref_base} = $ref_base;
			
			
			if ($OPT{keep_zyg}) {
				if ($sample_name) {
					$vcf_data{$var_key}{zyg} = \split(':',$fields[$sample_index]);
				} else {
					$vcf_data{$var_key}{zyg} = \@alleles;
				}
			}

			#Get the number of variant samples
			if ($OPT{keep_allele_freq}) {
				
				
				
				my ($allele_total) = $rest =~ /AN=(\d+)/;
				$allele_total = 0 if !defined $allele_total;
				
				
				$vcf_data{$var_key}{allele} = $total_ac . '/' . $allele_total;
			}
			
			#Get the frequency across all the samples containing the variant
			if ($OPT{mean_var_freq}) {
				
				my $stats = Statistics::Descriptive::Full->new();
				
				my @var_counts;
				
				#Strelka indels
				if ($gt_fields =~ /:TIR/) {
					for my $allele_str (@alleles) {
						my @gt_fields_data = split(':',$allele_str);
	    				my ($tir) = split(",",$gt_fields_data[$tir_index]);
	    				my ($tar) = split(",",$gt_fields_data[$tar_index]);
	    				my ($tor) = split(",",$gt_fields_data[$tor_index]);
	    				my ($sum) = $tir+$tar+$tor;
	    				my $sample_af = sprintf("%.4f",$tir/$sum);
	    				push @var_counts, "$tir/$sum";
	    				$stats->add_data($sample_af) if $sample_af != 0;
					}
				} elsif ($gt_fields =~ /:AU/) {
					#Strelka SNVs
					
					for my $allele_str (@alleles) {
						my @gt_fields_data = split(':',$allele_str);
	    				my $ref_count = my $var_count;
	    				if ($ref eq 'A') {
	    					($ref_count) = split (",",$gt_fields_data[$a_index]);
	    				} elsif ($ref eq 'C') {
	    					($ref_count) = split (",",$gt_fields_data[$c_index]);
	    				} elsif ($ref eq 'G') {
	    					($ref_count) = split (",",$gt_fields_data[$g_index]);
	    				}  elsif ($ref eq 'T') {
	    					($ref_count) = split (",",$gt_fields_data[$t_index]);
	    				}  
	    				
	    				if ($var eq 'A') {
	    					($var_count) = split (",",$gt_fields_data[$a_index]);
	    				} elsif ($var eq 'C') {
	    					($var_count) = split (",",$gt_fields_data[$c_index]);
	    				} elsif ($var eq 'G') {
	    					($var_count) = split (",",$gt_fields_data[$g_index]);
	    				}  elsif ($var eq 'T') {
	    					($var_count) = split (",",$gt_fields_data[$t_index]);
	    				}  
	    				my $sum = $var_count+$ref_count;
	    				my $sample_af = sprintf("%.4f",$var_count/$sum);
	    				push @var_counts, "$var_count/$sum";
	    				$stats->add_data($sample_af) if $sample_af != 0;
					}
				} else {
					for my $allele_str (@alleles) {
						my @gt_fields_data = split(':',$allele_str);
						#print Dumper \@gt_fields_data;
						my $var_count_str = my $sample_af = 0;
						
						#Some strange GTs don't have multiple values listed....
						if (@gt_fields_data > 1) {
	    					my @ads = split(',',$gt_fields_data[1]);
	    					my $sum = 0;
	    					for my $element (@ads) {
	    						$sum += $element if $element =~ /\d/;
							}	
							
							my $var_count;
							if (defined $ads[$zyg_num] && $ads[$zyg_num] =~ /\d/) {
								$var_count = $ads[$zyg_num];
							} else {
								$var_count = 0;	
							}
	    					if ($sum == 0) {
	    						$sample_af = 0;
	    					} elsif ($var_count == 0) {
	    						$sample_af = 0;
	    					} else {
	    						$sample_af = sprintf("%.4f",$ads[$zyg_num]/$sum)
	    					}
	    					$var_count_str =  "$var_count/$sum";
						}
	    				
	    				push @var_counts,$var_count_str;
	    				
	    				$stats->add_data($sample_af);
					}
					
				}
				my $mean = defined $stats->mean?$stats->mean:0;
				my $median = defined $stats->median?$stats->median:0;
				
				$vcf_data{$var_key}{mean_af} = sprintf("%.3f",$mean);
				$vcf_data{$var_key}{median_af} = sprintf("%.3f",$median);
				$vcf_data{$var_key}{var_read_counts} = join(",",@var_counts);
			}
			
			push @vcf_order, $var_key;
			$zyg_num++;
		}
		
    }
    return \%vcf_data;
}

#Gets variant key from either vcf or vep; standardises naming for loading into data structure
sub _get_variant_key {
	 my @args = @_;
	
	 my %args = @args;


    my @required_args = (
    					-chrom,
    					-first,
    					-ref_seq,
    					-var_seq,
    					-type
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
    
    my $ref = $args{-ref_seq};
    my $var = $args{-var_seq};
    my $first_coord = $args{-first};
    my $chr = $args{-chrom};
    my $type = $args{-type};
    
    my $start_coord = my $end_coord = my $bases;
    my $length_ref = length($ref);
    my $length_var = length($var);
    my $var_type;
    my $ref_base = "N/A";
    
    if ($type eq 'vcf') {
		if ($length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord + 1;
			$end_coord = $start_coord + $length_ref - $length_var - 1;				
			my $del_length = $length_ref - $length_var;
			#print "VCF R $ref L $del_length\n";
			
			$bases = '-'. substr($ref,1,$del_length);
		} elsif ($length_ref < $length_var) {
			#Add the ref length and var length difference to the coord 
			#$start_coord = $end_coord = $first_coord + 1;
			$var_type = 'INS';
			$start_coord = $end_coord = $first_coord;
			my $ins_length = $length_var - $length_ref;
			$bases = '+'.substr($var,1,$ins_length);
			$ref_base = substr($var,0,1);
		} elsif ($length_ref == $length_var && $length_ref != 1) {
			#Handling for cases like AT->AC ot ATATA->CTATA; turn first different base into an SNV
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			
			my @ref_bases = split("",$ref);
			my @var_bases = split("",$var);

			if ($ref_bases[0] ne $var_bases[0]) {
				#Simplest case where first base is different
				$bases = $ref_bases[0] . '->' .$var_bases[0];
			} else {
				my $until_snp_count = 0; #How far away from start we are
				my $length = @ref_bases;
				while ($length > 0) {
					if ($ref_bases[$until_snp_count] ne $var_bases[$until_snp_count]) {
						#print "Before $chr : $start_coord - $end_coord : $ref -> $var\n";
						$start_coord += $until_snp_count;
						$end_coord += $until_snp_count;
						$bases = $ref_bases[$until_snp_count] . '->' .$var_bases[$until_snp_count];
						#print "After $chr : $start_coord - $end_coord : $bases\n\n";
						last;
					}
					$until_snp_count++;
					$length--;
				}
				
			}
			
		} else {
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			$bases = $ref . '->' .$var;
		}
    	
    } else {
    	if ($ref eq '-' || $length_ref < $length_var) {
	    	$start_coord = $end_coord = $first_coord - 1;
			$var_type = 'INS';	
			my $ins_length = $length_var - $length_ref;
			$ins_length++ if $ref eq '-';
			$bases = '+'.substr($var,0,$ins_length);
		}  elsif ($var eq '-' || $length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord;
			my $del_length = $length_ref - $length_var;
			$end_coord = $start_coord + $length_ref - $length_var - 1;
			$del_length++ if $var eq '-';
			$end_coord++ if $var eq '-';
			$bases = '-'. substr($ref,0,$del_length);	
		} elsif ($length_ref == $length_var && $length_ref == 1) {
			#single snvs
			$var_type = 'SNV'; 
			$bases = $ref .'->'.$var;
			$start_coord = $end_coord = $first_coord;
		} else {
			modules::Exception->warning("ERROR: Can't identify var_type doesn't match any var type\n");
			#next;
		}
    }
    
	my $var_key = $chr . ':'.$start_coord .'-'.$end_coord .':' .$bases;
	return($var_key,$var_type,$ref_base);
	
}
