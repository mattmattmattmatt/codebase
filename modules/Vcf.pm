package Vcf;

use strict;
use Data::Dumper;
use modules::Exception;

sub new {
    my ($class, @args) = @_;

	my @required_args = (
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }


    my $self = bless {}, $class;
    
    my $organism = defined $args{-organism}?$args{-organism}:'human';
    
    $self->organism($organism);
    
    return $self;
}




#Get the organism
sub organism {
    my ($self, $organism) = @_;

    if (defined $organism) {
		$self->{'organism'} = $organism;
    } elsif (! defined $self->{'organism'}) {
		modules::Exception->throw("organism not set");
    }

    return $self->{'organism'};
}

#Parse a vcf file
sub parse_vcf {
    my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    
    my %vcf_data = ();
    

    open(VCF,$args{-vcf_file}) || modules::Exception->throw("Can't open file $args{-vcf_file}\n");
    
   	#Genotype quality index
   	my $gq_index = 0;
   	my $sample_index = 0;
	my $vcf_line = 1;

	#Other indices strelka needs
	my ($tir_index,$tar_index,$tor_index,$a_index,$c_index,$g_index,$t_index);
	my $strelka_indel = 0;
	my $strelka_snv = 0;
	
	



    while (<VCF>) {
	    next if /^#/;
    	next unless /\w/;
    	chomp;
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,@alleles) = split;
    	my $line = $_;
   		my @fields = split;
        my @gt_fields = split(':',$gt_fields) if $gt_fields;
		my @details = split(';',$rest);	

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


		my @vars = split(",",$var_str);

				

		my @ac_fields =  ();
		my $total_ac = 0;
		my @af_fields = ();


		for my $detail (@details) {
			if ($detail =~ /^AC=/) {
				$detail =~ s/AC=//;
				@ac_fields = split(",",$detail);
				for my $count (@ac_fields) {
					$total_ac += $count;
				} 
			}
			
			#If we have pre-filled AF (MGRB)
			if ($detail =~ /^AF=/) {
				$detail =~ s/AF=//;
				my @tmp_fields = split(",",$detail);
				for my $tmp (@tmp_fields) {
					push @af_fields, sprintf("%.5g", $tmp);
				}
			}
		
			
			
		}

		my $zyg_num = 0;
		
		
			
		for my $var ( @vars ) {
			$zyg_num++;

			next if $var eq '*'; #Due to upstream deletion
			next if $var eq $ref; #Bug with Tapestri / very rare
			my ($var_key,$var_type,$ref_base) = _get_variant_key(-type=>'vcf',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var);
			next unless $var_key =~ /(\d+)-(\d+):.+/;
			my ($start,$end) = $var_key =~ /(\d+)-(\d+)/;

			$vcf_data{$var_key}{var_count} = shift @ac_fields;
			$vcf_data{$var_key}{zyg_count} = $zyg_num;
			
			if ($qual eq '.') {
				$vcf_data{$var_key}{qual} = "N/A";				
			} else {
				$vcf_data{$var_key}{qual} = $qual;
			}
			$vcf_data{$var_key}{type} = $var_type;
			$vcf_data{$var_key}{ref_base} = $ref_base;
			$vcf_data{$var_key}{zyg} = \@alleles;

			#Get the number of variant samples
			my ($allele_total) = $rest =~ /AN=(\d+)/;
			$allele_total = 0 if !defined $allele_total;
			$vcf_data{$var_key}{allele} = $total_ac . '/' . $allele_total;
			
			#Get the frequency across all the samples containing the variant
				
			my @stats = ();
			
			my @var_counts;
			
			if (@alleles) {

			
				if ($gt_fields =~ /:TIR/) {
						for my $allele_str (@alleles) {
							next unless $allele_str =~ /\d/;
							my @gt_fields_data = split(':',$allele_str);
		    				my ($tir) = split(",",$gt_fields_data[$tir_index]);
		    				my ($tar) = split(",",$gt_fields_data[$tar_index]);
		    				my ($tor) = split(",",$gt_fields_data[$tor_index]);
		    				my ($sum) = $tir+$tar+$tor;
		    				next if $sum == 0;
		    				my $sample_af = sprintf("%.4f",$tir/$sum);
		    				push @var_counts, "$tir/$sum";
		    				push @stats,$sample_af if $sample_af != 0;
						}
				} elsif ($gt_fields =~ /:AU/) {
						#Strelka SNVs
						
						for my $allele_str (@alleles) {
							next unless $allele_str =~ /\d/; 
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
		    				next if $sum == 0;
		    				my $sample_af = sprintf("%.4f",$var_count/$sum);
		    				push @var_counts, "$var_count/$sum";
		    				push @stats, $sample_af if $sample_af != 0;
						}
					} else {
							for my $allele_str (@alleles) {
								my @gt_fields_data = split(':',$allele_str);
								#print Dumper \@gt_fields_data;
								my $var_count_str = my $sample_af = 0;
									
								#Some strange GTs don't have multiple values listed....
								if (@gt_fields_data > 1) {
									#Strings like 10,5,2
				    				my @ads = split(',',$gt_fields_data[1]);
				    				#Total allele count
				    				my $sum = 0;
				    				for my $element (@ads) {
				    					$sum += $element if $element =~ /\d/;
									}	
										
									my $var_count;
									#Variant count for that allele
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
				    				
				    				
					    			push @stats,$sample_af if $var_count>1;
					    			push @var_counts,$var_count_str;
								}
				    			
				    		}
				    	
				    	}
				    my $mean = my $median = 0;
				
					if (@stats) {
						($mean,$median) = _mean_median(-numbers=>\@stats);
					}
				
				
				
					$vcf_data{$var_key}{mean_af} = sprintf("%.3f",$mean);
					$vcf_data{$var_key}{median_af} = sprintf("%.3f",$median);
					$vcf_data{$var_key}{var_read_counts} = join(",",@var_counts);
			   				
				
			} elsif (@af_fields) {
				$vcf_data{$var_key}{mean_af} = $af_fields[$zyg_num-1];
				$vcf_data{$var_key}{median_af} = $af_fields[$zyg_num-1];
				$vcf_data{$var_key}{var_read_counts} = 'N/A';
			} else {
				modules::Exception->throw("No AF info available $chr $first_coord\n");
			}
				
			
			push @{$self->{vcf_order}}, $var_key;
		}
		$vcf_line++;
		if ($vcf_line %1000000 == 0) {
			#print join("\t",$chr,$first_coord) . "\n";
		}
    }
    $self->{data}{$args{-vcf_file}} = \%vcf_data;
}

#Write normalised file
sub write_normalised {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_out,
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    
    
    
    open(FILE,">$args{-vcf_out}") || modules::Exception->throw("Can't open file to write $args{-vcf_out}\n");
	
	
	my $vcf_data = $self->get_vcf(-vcf_file=>$args{-vcf_file});
	
	for my $var_key (@{$self->{vcf_order}}) {	
		
		my ($chr,$coord,$event) = split(':',$var_key);
		my ($start,$end) = split('-',$coord);
		
		$chr =~ s/chr//;
	
		my $var_count = 0;
		if (exists $vcf_data->{$var_key} && defined $vcf_data->{$var_key}{var_count}) {
			$var_count = $vcf_data->{$var_key}{var_count};
		}	
	
		my $zyg_count = 0;
	    if (exists $vcf_data->{$var_key} && defined $vcf_data->{$var_key}{zyg_count}) {
	    	$zyg_count = $vcf_data->{$var_key}{zyg_count};
	    }
	
	
		my $vcf_str  = $vcf_data->{$var_key}{type}.';'.$var_key.';Q='.$vcf_data->{$var_key}{qual}. ';AC='.$var_count . ';ZC='.$zyg_count.";ALLELE=".$vcf_data->{$var_key}{allele}.";MEANAF=".$vcf_data->{$var_key}{mean_af}.";MEDAF=".$vcf_data->{$var_key}{median_af}.";VAR_READ_COUNTS=".$vcf_data->{$var_key}{var_read_counts};
		

		if ($vcf_data->{$var_key}{type} eq 'INS') {
			$vcf_str .= ";REF=".$vcf_data->{$var_key}{ref_base};
		}
		
		
		my $zyg_str = join("\t",@{$vcf_data->{$var_key}{zyg}});
		print FILE join("\t", 
							$chr,
							$start,
							$end,
							$vcf_str,
							$zyg_str
							)."\n";
		
	}
	close FILE;
    
}


#Get vcf data for a file
sub get_vcf {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $self->{data}{$args{-vcf_file}} doesn't exist");	
    }
    return $self->{data}{$args{-vcf_file}};
}



#Check that a vcf file is 'complete'
sub check_vcf {
	my ($self,@args) = @_;
	my @required_args = (
			             -vcf_file,
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
	
	my $vcf = $args{-vcf_file};
	if ( !-e $vcf) {
		modules::Exception->throw("File $vcf doesn't exist");
	}
	my $organism = $self->{organism};
	
	if ($organism ne 'human' && $organism ne 'mouse') {
		modules::Exception->throw("ERROR: organism must be human or mouse");
	}
	

	
	

	return 1;
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


#Simple stats to avoid installing Statistics::Descriptive as not standard
sub _mean_median { 
	my @args = @_;
	my %args = @args;
	
	
	my @numbers = @{$args{-numbers}};
	
	my $sum = 0;
	my @sorted_numbers = sort {$a<=>$b} @numbers;
	foreach (@numbers) { $sum += $_; }
	
    my $mean = $sum/@numbers;
    my $median_pos = int(@numbers / 2);

    my $median;
    if ( @numbers % 2 == 1) {
        $median = $sorted_numbers[$median_pos];
    } else {
        my $median_pos2 = $median_pos - 1;
        $median = ($sorted_numbers[$median_pos] + $sorted_numbers[$median_pos2]) / 2;
    }
    
    
  	return ($mean,$median);
}
return 1;

=comment
package modules::Vcf;

use strict;
use Data::Dumper;
use modules::Exception;
use modules::Pipeline;

sub new {
    my ($class, @args) = @_;

	my @required_args = (
			             -sample_name
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }


    my $self = bless {}, $class;
    
    my $sample_name = $args{-sample_name};
    my $sample_type = modules::Pipeline::get_sample_type(-sample_name=>$sample_name);
    my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);
    my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_name);
    my $pipe_conf = modules::Pipeline::get_pipe_conf();
    my $organism = $pipe_conf->read($source_type,'organism');
    
    $self->source_type($source_type);
	$self->sample_type($sample_type);
    $self->organism($organism);
    $self->sequence_type($sequence_type);
    
    return $self;
}

#Get the sample_type
sub sample_type {
    my ($self, $sample_type) = @_;

    if (defined $sample_type) {
		$self->{'sample_type'} = $sample_type;
    } elsif (! defined $self->{'sample_type'}) {
		modules::Exception->throw("sample_type not set");
    }

    return $self->{'sample_type'};
}

#Get the sample_type
sub sequence_type {
    my ($self, $sequence_type) = @_;

    if (defined $sequence_type) {
		$self->{'sequence_type'} = $sequence_type;
    } elsif (! defined $self->{'sequence_type'}) {
		modules::Exception->throw("sequence_type not set");
    }

    return $self->{'sequence_type'};
}

#Get the sample_type
sub source_type {
    my ($self, $source_type) = @_;

    if (defined $source_type) {
		$self->{'source_type'} = $source_type;
    } elsif (! defined $self->{'source_type'}) {
		modules::Exception->throw("source_type not set");
    }

    return $self->{'source_type'};
}

#Get the organism
sub organism {
    my ($self, $organism) = @_;

    if (defined $organism) {
		$self->{'organism'} = $organism;
    } elsif (! defined $self->{'organism'}) {
		modules::Exception->throw("organism not set");
    }

    return $self->{'organism'};
}

#Parse a vcf file
sub parse_vcf {
    my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    open(VCF,$args{-vcf_file}) || modules::Exception->throw("Can't open file $args{-vcf_file}\n");
    my %variant_data = ();
    
   	my $tumour_flag = $args{-tumour_flag};
    
    
    while (<VCF>) {
    	next if /^#/;
    	
    	if ($_ !~ /^chr[0-9XYM]/ && $_ !~ /^[0-9XYM]/) {
    		modules::Exception->throw("ERROR: Error with VCF format at line $_");
    	}
    	
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,$pass,$rest,$gt_fields,$normal_gt,$tumour_gt) = split("\t");
		$chr =~ s/chr//;    
		if (!$qual) {
    		modules::Exception->throw("ERROR: Error with VCF format at line $_");
    	}

		if ($ref =~ /N/i || $var_str =~ /N/i) {
			next;
		}

		my @normal_fields = split(':',$normal_gt);
    	my @tumour_fields = split(':',$tumour_gt) if $tumour_gt;
		
		my $var_type;
		if ($tumour_flag) {
			if ($normal_fields[0] eq '0/0') {
				$var_type =  'SOMATIC';
			} elsif ($normal_fields[0] eq $tumour_fields[0]) {
				$var_type = 'GERMLINE';
			} elsif ($normal_fields[0] eq '0/1' && ($tumour_fields[0] eq '0/0' || $tumour_fields[0] eq '1/1')) {
				$var_type = 'LOH';
			} else {
				$var_type = 'COMPLEX';
			}
		} else {
			#Can't tell at this point when it's not matched
			$var_type = 'UNKNOWN';
		}
		
		
    	#Get the CLR score
    	my $clr_str = 'NO_CLR';
    	if (/CLR=(\d+)/) {
    		
    		if ($normal_fields[0] eq '0/0') {
    			#Check it's the tumour that is variant
    			$clr_str = 'REF_NORM'.$1;
    		} elsif ($tumour_fields[0] eq '0/0') {
	    		#If it's the tumour that matches reference
    			$clr_str = 'REF_TUM'.$1;
    		} else {
				$clr_str = 'REF_NONE'.$1;				
    		}
    	} 
    	
		my @vars = split(",",$var_str);
		
		for my $var ( @vars ) {
			my $variant_type;
			my $start_coord;
			my $end_coord;
			my $bases;
			my $length_ref = length($ref);
			my $length_var = length($var);
			
			if (length($ref) > length($var)) {
				$variant_type = 'DEL';
				$start_coord = $first_coord + 1;
				$end_coord = $start_coord + $length_ref - $length_var - 1;				
				my $del_length = $length_ref - $length_var;
				$bases = substr($ref,1,$del_length);
			} elsif (length($ref) < length($var)) {
				$variant_type = 'INS';
				#Add the ref length and var length difference to the coord -> was a bug
				#$start_coord = $end_coord = $first_coord + 1;
				$start_coord = $end_coord = $first_coord;
				my $ins_length = $length_var - $length_ref;
				$bases = substr($var,1,$ins_length);
			} else {
				$variant_type = 'SNV';
				$start_coord = $end_coord = $first_coord;
				$bases = $ref . '->' .$var;
			}
			my $key = 'Q'.$qual . '^^^'.$clr_str .'^^^'. $var_type;
			#VQSR pass
			if ($pass eq 'PASS') {
				$key .= '^^^PASS';
			}
			#print "$bases $variant_type\n";
			$variant_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $key;
		}
    }
    
    $self->{data}{$args{-vcf_file}} = \%variant_data;
}

#Get vcf data for a file
sub get_vcf {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $self->{data}{$args{-vcf_file}} doesn't exist");	
    }
    return $self->{data}{$args{-vcf_file}};
}

#Apply filters to vcf data and return passed lines only
sub filter_vcf {
	my ($self, @args) = @_;

    my @required_args = (
			             -vcf_file,
			             -snv_depth_file,
			             -indel_depth_file,
			             -tumour_flag
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $self->{data}{$args{-vcf_file}} doesn't exist");	
    }
    
    if ( !-e $args{-snv_depth_file} ) {
    	modules::Exception->throw("File $args{-snv_depth_file} doesn't exist");	
    }
    
    if ( !-e $args{-indel_depth_file} ) {
    	modules::Exception->throw("File $args{-indel_depth_file} doesn't exist");	
    }
    
   	

    my $tumour_flag = $args{-tumour_flag};
    my $gatk = defined $args{-gatk}?1:0;
    my %indel_depth = ();
    my %snv_depth = ();
    
    #Get the depth from the merge_vcf files
    open(SNV,"$args{-snv_depth_file}") || modules::Exception->throw("Can't open file $args{-snv_depth_file}\n");
    while (<SNV>) {
    	my @fields = split("\t");
    	my ($depth) = $fields[3] =~ /D(\d+)$/;
    	$snv_depth{$fields[0]}{$fields[1]} = $depth;
    }
    
    open(INDEL,"$args{-indel_depth_file}") || modules::Exception->throw("Can't open file $args{-indel_depth_file}\n");
    while (<INDEL>) {
    	my @fields = split("\t");
    	my ($depth) = $fields[3] =~ /D(\d+)$/;
	    $indel_depth{$fields[0]}{$fields[1]} = $depth;
    }
    
    
    my $pipe_config = modules::Pipeline::get_pipe_conf();
	my $sample_type = $self->sample_type;
    my $source_type = $self->source_type;
    
    my $min_snv_quality = $pipe_config->read($source_type,'cutoffs','snv_quality_cutoff');
	my $min_indel_quality = $pipe_config->read($source_type,'cutoffs','indel_quality_cutoff');
	my $min_clr = $pipe_config->read($source_type,'cutoffs','clr_cutoff');
	my $min_depth = $pipe_config->read($source_type,'cutoffs','min_variant_cutoff');
	my $max_depth = $pipe_config->read($source_type,'cutoffs','max_variant_cutoff');
	my %snv_types = ();
	if ($pipe_config->exists($source_type,'db_types')) {
		%snv_types = map {$_ => 1} split(",",$pipe_config->read($source_type,'db_types'));
	} else {
		%snv_types = map {$_ => 1} split(",",$pipe_config->read('common','db_types'));	
	}
	
	my $sequence_type = $self->sequence_type;
	
	my %filter_vcf_data = ();
	my $all_vcf_data = $self->get_vcf(-vcf_file => $args{-vcf_file});
	if (exists $all_vcf_data->{SNV}) {
		for my $chr (sort keys %{$all_vcf_data->{SNV}}) {
			for my $start_coord (sort {$a<=>$b} keys %{$all_vcf_data->{SNV}{$chr}}) {
				for my $end_coord (keys %{$all_vcf_data->{SNV}{$chr}{$start_coord}}) {
					for my $bases (keys %{$all_vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}}) {
						my $rest = $all_vcf_data->{SNV}{$chr}{$start_coord}{$end_coord}{$bases};
						my $depth = $snv_depth{$chr}{$start_coord};
						my $new_rest = $rest . '^^^D' . $depth;
						#Here we need to consider CLR score as well
						my ($quality_str,$clr_str,$type) = split('\^\^\^',$rest);
						
						if (!exists $snv_types{$type} && !exists $snv_types{'ALL'}) {
							#Only pass SNVs we're interested in
							next;
						}
						my ($quality) = $quality_str =~ /(\d+)/;

						if ($tumour_flag) {
							
							#Now filter CLR score
							next if $clr_str eq 'NO_CLR';
							
							my ($clr_score) = $clr_str =~ /(\d+)/;
							if ($sequence_type eq 'targeted') {
								if ($depth >= $min_depth && $quality >= $min_snv_quality && $clr_score >= $min_clr ) {
									$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
								}		
							} else {
								if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_snv_quality && $clr_score >= $min_clr ) {
									$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
								}
							}
													
						} elsif ($gatk) {
							if ($rest =~ /PASS/) {
	                        	$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
							}
							
						} else {
							if ($sequence_type eq 'targeted') {
								if ($depth >= $min_depth && $quality >= $min_snv_quality) {
									$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
								}		
							} else {
								if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_snv_quality) {
									$filter_vcf_data{SNV}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
								}
							}
						}
					}
				}
			}
		}
	}
	
	
	my %indel_types = ();
	if ($pipe_config->exists($source_type,'db_types')) {
		%indel_types = map {$_ => 1} split(",",$pipe_config->read($source_type,'db_types'));
	} else {
		%indel_types = map {$_ => 1} split(",",$pipe_config->read('common','db_types'));	
	}
	my @indel_type = qw(DEL INS);
	
	for my $variant_type (@indel_type) {
		if (exists $all_vcf_data->{$variant_type}) {
			for my $chr (sort keys %{$all_vcf_data->{$variant_type}}) {
				for my $start_coord (sort {$a<=>$b} keys %{$all_vcf_data->{$variant_type}{$chr}}) {
					for my $end_coord (keys %{$all_vcf_data->{$variant_type}{$chr}{$start_coord}}) {
						for my $bases (keys %{$all_vcf_data->{$variant_type}{$chr}{$start_coord}{$end_coord}}) {
							my $rest = $all_vcf_data->{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases};
							my $depth = $indel_depth{$chr}{$start_coord};
							my $new_rest = $rest . '^^^D' . $depth;
							#Here we need to consider CLR score as well
							my ($quality_str,$clr_str,$type) = split('\^\^\^',$rest);
					
							if (!exists $indel_types{$type}  && !exists $snv_types{'ALL'}) {
								#Only pass SNVs we're interested in
								next;
							}
							my ($quality) = $quality_str =~ /(\d+)/;

							if ($tumour_flag) {
								
								#Now filter CLR score
								next if $clr_str eq 'NO_CLR';
								
								my ($clr_score) = $clr_str =~ /(\d+)/;
								if ($sequence_type eq 'targeted') {
									if ($depth >= $min_depth && $quality >= $min_indel_quality && $clr_score >= $min_clr ) {
										$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
									}
								} else {		
									if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_indel_quality && $clr_score >= $min_clr ) {
										$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
									}
								}
							} elsif ($gatk) { 
                            	if ($rest =~ /PASS/) {
                                	$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
                                }
							} else {
								#print "$variant_type $chr $start_coord $end_coord if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_indel_quality)\n";
								if ($sequence_type eq 'targeted') {
									if ($depth >= $min_depth && $quality >= $min_indel_quality) {
										$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
									}
								} else {		
									if ($depth >= $min_depth && $depth <= $max_depth && $quality >= $min_indel_quality) {
										$filter_vcf_data{$variant_type}{$chr}{$start_coord}{$end_coord}{$bases} = $new_rest;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	return \%filter_vcf_data;
}

#Check that a vcf file is 'complete'
sub check_vcf {
	my ($self,@args) = @_;
	my @required_args = (
			             -vcf,
			             -chromosome
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
	
	my $chr_args = $args{-chromosome};
	my $vcf = $args{-vcf};

	my $organism = $self->{organism};
	
	my $chr_count;
	if ($organism eq 'human') {
		$chr_count = 23;
	} elsif ($organism eq 'mouse') {
		$chr_count = 20;
	} else {
		modules::Exception->throw("ERROR: organism must be human or mouse");
	}
	

	my %human_chr_length = (
						1 => 249250621,
						2 => 243199373, 	
						3 => 198022430, 
						4 => 191154276, 
						5 => 180915260, 
						6 => 171115067, 
						7 => 159138663, 
						8 => 146364022, 
						9 => 141213431, 
						10 => 135534747, 
						11 => 135006516, 
						12 => 133851895, 
						13 => 115169878, 
						14 => 107349540, 
						15 => 102531392, 
						16 => 90354753, 
						17 => 81195210, 
						18 => 78077248, 
						19 => 59128983, 
						20 => 63025520, 
						21 => 48129895, 
						22 => 51304566, 
						'X' => 155270560, 
						'Y' => 59373566,
						);
	
	my %mouse_chr_length = (
						1 => 195471971,
						2 => 182113224, 	
						3 => 160039680, 
						4 => 156508116, 
						5 => 151834684, 
						6 => 149736546, 
						7 => 145441459, 
						8 => 129401213, 
						9 => 124595110, 
						10 => 130694993, 
						11 => 122082543, 
						12 => 120129022, 
						13 => 120421639, 
						14 => 124902244, 
						15 => 104043685, 
						16 => 98207768, 
						17 => 94987271, 
						18 => 90702639, 
						19 => 61431566, 
						'X' => 171031299, 
						'Y' => 91744698,
							);
	
	return 1 if $chr_args eq 'Y'; #Not covered for females
	
	my $sequence_type = $self->sequence_type;
	
	
	#Here we have to do the calculation for all the chromosomes at once
	open(VCF,$vcf) || modules::Exception->throw("ERROR: Can't open vcf file $vcf to check");
	my %chr_stats = ();
	
	
	while (<VCF>) {
		next unless (/^chr[0-9XYM]/ || /^[0-9XYM]/);
		my ($chr,$coord) = split /[\t\s]+/;
		$chr =~ s/chr//;


		#initialise values
		if (!exists $chr_stats{$chr}) {
			$chr_stats{$chr}{min} = $coord;
			$chr_stats{$chr}{max} = $coord;
		}

		if ($chr_stats{$chr}{max} < $coord) {
			$chr_stats{$chr}{max} = $coord;
		}

		$chr_stats{$chr}{count}++;

	}

	#For targeted can get empty vcfs
	return 1 if $sequence_type eq 'targeted';

	#If no data in files
	if (!keys %chr_stats) {
		modules::Exception->throw("ERROR: Didn't get any data from vcf file");
	}

	my $chr_covered = keys %chr_stats;
	
	
	#Skip rest of tests for targeted cases
	return 1 if $sequence_type eq 'exome';

	#If we're expecting all the chromosomes and we didn't get enough
	if ($chr_covered < $chr_count && $chr_args eq 'all') {
		modules::Exception->throw("ERROR: Only covered $chr_covered chromosomes and expected $chr_count chromosomes");
	}
	
	for my $chr (sort keys %chr_stats) {
		next if $chr eq 'Y' || $chr eq 'X'; #PAR regions make this test fail sometimes
		
		#If only checking single chromosome skip other chromosomes
		if ($chr_args ne 'all') {
			next unless $chr_args eq $chr;
		}
		
		my $chr_length;
		if ($organism eq 'human') {
			$chr_length = $human_chr_length{$chr};
		} elsif ($organism eq 'mouse') {
			$chr_length = $mouse_chr_length{$chr};
		}
		
		my $chr_max_coord = $chr_stats{$chr}{max};
		my $chr_min_coord = $chr_stats{$chr}{min};
		my $chr_total_count = $chr_stats{$chr}{count};
		
		my $end_proportion = $chr_max_coord / $chr_length;
		my $start_proportion = $chr_min_coord / $chr_length;
		my $snv_proportion =  $chr_total_count / $chr_length;
		
		#If the ratio of lines to chr_length is < 1:10000 then we likely have a problem; only for genome cases where we expect lots of cover all over
		if ($snv_proportion < 0.0001 && $sequence_type eq 'genome') {
			print STDERR "WARNING: Vcf has $chr_total_count lines and chr length $chr is $chr_length ($snv_proportion < 0.0001)\n";
			return 0;
		}
		
		#If we're >1% of the way to the start of the chromosome then we likely have a problem; skip chr with huge N start regions
		if ($organism eq 'human') {
			if ($start_proportion > 0.01 && $chr !~ /1[345]/ && $chr !~ /2[12]/) {
				print STDERR "WARNING: Vcf has start coord $chr_min_coord and chr $chr length is $chr_length ($start_proportion > 0.01)\n";
				return 0;
			}
		} elsif ($organism eq 'mouse') {
			#Mouse has huge number of N bases in every chromosome
			if ($start_proportion > 0.1) {
				print STDERR "WARNING: Vcf has start coord $chr_min_coord and chr $chr length is $chr_length ($start_proportion > 0.01)\n";
				return 0;
			}
		}
		
		#If we're not at least 99% of the way to the end of the chromosome then we likely have a problem; for all cases
		if ($end_proportion < 0.99) {
			print STDERR "WARNING: Vcf has last coord $chr_max_coord and chr $chr length is $chr_length ($end_proportion < 0.99)\n";
			return 0; 
		} 
		
		
		#print "$organism $chr\ncover $snv_proportion < 0.0001  start $start_proportion > 0.01  end $end_proportion < 0.99\n\n";
	}
		
	
	

	return 1;
}


return 1;
=cut
