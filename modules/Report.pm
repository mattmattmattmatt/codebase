package modules::Report;

use strict;
use modules::FilteredSNV;
use modules::Adaptors::SNV;
use modules::Adaptors::SNV_Row;
use modules::Adaptors::Variant_Row;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::BulkInsert;
use modules::ReportFilter;
use modules::FilteredVariant;
use modules::PileupLine;
use modules::Exception;
use modules::Utils;
use modules::ConfigXML;
use modules::VariantXML;
use modules::Pipeline;
use Data::Dumper;
use File::Basename;
use Bio::EnsEMBL::Registry;
use Env qw($ENSEMBL_REGISTRY);
my $EXON_NS_FILTER = 'filter_exon_ns';
my $EXON_FILTER = 'filter_exon';
my $DBSNP_FILTER = 'filter_dbsnp_snv';
my $DBSNP_FILTER_INDEL = 'filter_dbsnp_indel';
my $EXAC_FILTER = 'filter_exac_snv';
my $EXAC_FILTER_INDEL = 'filter_exac_indel';
my $GNOMAD_FILTER = 'filter_gnomad_snv';
my $GNOMAD_FILTER_INDEL = 'filter_gnomad_indel';
my $CLINVAR_FILTER = 'filter_clinvar_snv';
my $CLINVAR_FILTER_INDEL = 'filter_clinvar_indel';
my $REGULATORY_CUSTOM_FILTER = 'filter_regulatory_custom';
my $FILTER_MIRNA = 'filter_mirna';
my $COSMIC_FILTER = 'filter_cosmic';
my $SPLICE_FILTER = 'filter_splicesite';
my $VEP_FILTER = 'filter_vep';
my $FILTER_GENE = 'filter_gene';
my $FILTER_GENE_BP = 'filter_gene_breakpoint';
my $FILTER_SV_EXON = 'filter_sv_exon';
my $FILTER_SV_EXON_BP = 'filter_sv_exon_breakpoint';

my $FILTER_DGV = 'filter_dgv';
my $PUBMED_URL='http://www.ncbi.nlm.nih.gov/pubmed/';

sub new {
    my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

    my @required_args = (
			             -run, 
						 -gene_mapper,
						 -sample_type,
						 -source_type,
						 -gene_col_name
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
	
	#Run id
    $self->run($args{'-run'});
    
    my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$self->run->id);
    $self->sample_name($sample_name);
	#Gene name mapping
	$self->gene_mapper($args{'-gene_mapper'});
	#Sample type
	$self->sample_type($args{-sample_type});
	#Gene column name
	$self->gene_col_name($args{-gene_col_name});
	#Source type
	$self->source_type($args{-source_type});
	my $pipe_config = modules::Pipeline::get_pipe_conf();
	#set the organism from source_type
	$self->organism($pipe_config->read($args{-source_type},'organism'));
	
	my $paired_variants = $pipe_config->read($args{-source_type},'paired_variants');
	#Cancer samples; human_related samples aren't counted in this category
	$self->{paired_variants} = $paired_variants;


	if (defined $args{-confdir}) {
		$self->{genome} = 1;
		if (!-d $args{-confdir}) {
			modules::Exception->throw("ERROR: confdir $args{-confdir} doesn't exist");
		}
		$self->confdir($args{-confdir});
		
	} else {
		$self->{genome} = 0;
	}
	
	#Optional extra lists we've filtered out    
	my @filter_files = defined $args{-filter_files}?@{$args{-filter_files}}:();
	$self->filter_files(\@filter_files) if @filter_files;
	
	#Used to differentiate indel from snv reports when different behavior is required
	if (defined $args{-indel}) {
		$self->{variant_type} = 'indel';
	} elsif (defined $args{-sv}) {
		$self->{variant_type} = 'sv';
		$self->{sv_caller} = $args{'-sv_caller'};
	} else {
		$self->{variant_type} = 'snv';
		
		#These are snv specific
		if (defined $args{-polyphen_file}) {
			$self->_parse_polyphen($args{-polyphen_file});
		}	

		#This is for the polyphen entries that don't have scores
		if (defined $args{-polyphen_info_file}) {
			$self->_parse_polyphen_info($args{-polyphen_info_file})
		}
	}
	
#	my $registry = 'Bio::EnsEMBL::Registry';
#	$registry->load_all($ENSEMBL_REGISTRY);
#	my $slice_adaptor = $registry->get_adaptor('human', 'core', 'Slice');
#	$self->{adaptor} = $slice_adaptor;	
	
    return $self;
}



#Parse the polyphen file to add these scores to the columns -> these are the polyphen entries that have been run
sub _parse_polyphen {
	my $self = shift;
	my $polyphen_file = shift;
	my %polyphen = ();
	open(FILE,"$polyphen_file") || modules::Exception->throw("Can't open file $polyphen_file\n");
	while (<FILE>) {
		my ($chr,$coord,undef,$score,$prediction,$rest) = split("\t");
		my @rest_fields = split(':',$rest);
		my $nt_change = $rest_fields[2];
		$polyphen{$chr}{$coord}{$nt_change} = $prediction  . ',' . $score;
	}
	close FILE;
	$self->{polyphen} = \%polyphen; 
}

#Parse the polyphen info file -> this contains info to allow user to run polyphen themselves
sub _parse_polyphen_info {
	my $self = shift;
	my $polyphen_info_file = shift;
	my %polyphen_info = ();
	open(FILE,"$polyphen_info_file") || modules::Exception->throw("Can't open file $polyphen_info_file\n");
	while (<FILE>) {
		chomp;
		my ($chr,$coord,undef,$poly_info) = split("\t");
		$polyphen_info{$chr}{$coord} = $poly_info;
	}
	close FILE;
	$self->{polyphen_info} = \%polyphen_info; 
}

#Subroutine loads up all the filters used for input/reporting/pass,etc; these conditions are based on sample_type
sub load {
	my $self = shift;
	my %args = @_;
	
    	my $sample_type = $self->{sample_type};

	my $config;

    	#Parse the xml
    	if (defined $args{-report_xml}){
		$config = modules::Pipeline->get_report_conf(-report_xml=>$args{-report_xml});
	} else {
		$config = modules::Pipeline->get_report_conf();
	}

	my $variant_type = $self->{variant_type};

	my $conf_header;
	if ($variant_type eq 'snv') {
		$conf_header = 'snv';
	} elsif ($variant_type eq 'sv') {
		$conf_header = 'sv';
	} else {
		$conf_header = 'indel';
	}


	if (!$config->exists($conf_header,'sample_types',$sample_type)) {
		modules::Exception->throw("ERROR: no sample_type config entry for $sample_type");
	}
    

    #Load up the data struct

    
    #Input filter; db entries or files that are input to the reports (currently filter_exon, filter_splicesite)
    
    my @input_filters = ();
    if ($config->exists($conf_header,'sample_types',$sample_type,'input_filters')) {
    	@input_filters = split(",",$config->read($conf_header,'sample_types',$sample_type,'input_filters'));
    } else {
    	@input_filters = split(",",$config->read($conf_header,'common','input_filters'));
    }
   
    if ($self->{genome}) {
		$self->{primary_filters} = \@input_filters;    	
    } else {
	    my @input_filters_db;
	    for my $input_filter ( @input_filters ) {
		    my ($filter_db) = modules::Adaptors::Filter->search(name => $input_filter);
		    if (!defined $filter_db) {
		    	modules::Exception->throw("ERROR: filter $input_filter doesn't exist in db");
		    }
		    push @input_filters_db, $filter_db;
	    	
		}
	    $self->{primary_filters} = \@input_filters_db;
    }
    
    
    #Filter with snv_filter info or info in file
    my @filter_info = ();
    if ($config->exists($conf_header,'sample_types',$sample_type,'filter_info')) {
    	@filter_info = split(",",$config->read($conf_header,'sample_types',$sample_type,'filter_info'));
    } else {
    	@filter_info = split(",",$config->read($conf_header,'common','filter_info'));
    }
    
    if ($self->{genome}) {
    	$self->{filter_info} = \@filter_info;
    } else {
	    my @filter_info_db;
	    for my $filter_info ( @filter_info ) {
			my ($filter_db) = modules::Adaptors::Filter->search(name => $filter_info);
			if (!defined $filter_db) {
		    	modules::Exception->throw("ERROR: filter $filter_info doesn't exist in db");
		    }
			push @filter_info_db, $filter_db;
		}	
	    $self->{filter_info} = \@filter_info_db; 
    }
    
	
	
	    
    #Extra filters to report; don't need db object
    my @extra_filters = ();
    if ($config->exists($conf_header,'sample_types',$sample_type,'extra_filters')) {
    	@extra_filters = split(",",$config->read($conf_header,'sample_types',$sample_type,'extra_filters'));
    } else {
    	@extra_filters = split(",",$config->read($conf_header,'common','extra_filters'));
    }
   
#    my @extra_filters_db;
#    for my $extra_filter ( @extra_filters ) {
#		my ($filter_db) = modules::Adaptors::Filter->search(name => $extra_filter);
#		if (!defined $filter_db) {
#	    	modules::Exception->throw("ERROR: filter $extra_filter doesn't exist in db");
#	    }
#		push @extra_filters_db, $filter_db;
#	}	
    $self->{filter_list_order} = \@extra_filters; 
	
	
	#Report filter data structure
	
	
	my %report_set = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'summary_report')) {
		%report_set = map {$_ => 1} split(",",$config->read($conf_header,'sample_types',$sample_type,'summary_report'));
	} else {
		%report_set = map {$_ => 1} split(",",$config->read($conf_header,'common','summary_report'));	
	}
	
	$self->{report_set} = \%report_set;    

	#Filter pass conditions 
	my $rules;
	if ($config->exists($conf_header,'sample_types',$sample_type,'pass_conditions','rule')) {
		$rules = $config->read($conf_header,'sample_types',$sample_type,'pass_conditions','rule');
	} else {
		$rules = $config->read($conf_header,'common','pass_conditions','rule');		
	}
	my @rules = ();
	#Get the rules	
	if (ref($rules) eq 'ARRAY'){ # Cope with steps being a single step or an array of steps
    	@rules = @$rules;
	} else {
    	@rules = ($rules);	
	}
	
	my %rule_set = ();
	for my $rule ( @rules ) {
	    my ($name,$condition) = split(':',$rule);
	    $rule_set{$name} = $condition;
	    
	}
	
	$self->{filter_set} = \%rule_set;

	my @annotation_headers = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'annotations')) {
		@annotation_headers = split(",",$config->read($conf_header,'sample_types',$sample_type,'annotations'));
	} else {
		@annotation_headers = split(",",$config->read($conf_header,'common','annotations'));
	}

	my @run_headers = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'run_headers')) {
		@run_headers = split(",",$config->read($conf_header,'sample_types',$sample_type,'run_headers'));
	} else {
		@run_headers = split(",",$config->read($conf_header,'common','run_headers'));
	}

	my @extra_headers = ();
	if ($config->exists($conf_header,'sample_types',$sample_type,'extra_headers')) {
		@extra_headers = split(",",$config->read($conf_header,'sample_types',$sample_type,'extra_headers'));
	} else {
		@extra_headers = split(",",$config->read($conf_header,'common','extra_headers'));
	}
	
	$self->{extra_headers} = \@extra_headers;
	$self->{run_headers} = \@run_headers;
	$self->{annotation_headers} = \@annotation_headers;

	#Header lookup for columns
	$self->headers();
}

#Add the additional filter files
sub filter_files {
	my $self = shift;
	my ($filter_files) = @_;
	my %filter_files = ();
	#Handle any additional gene lists we want to filter
	my $count = 1;
	for my $filter_file ( @{$filter_files} ) {
	    if ( !-e $filter_file ) {
	    	modules::Exception->throw("File $filter_file doesn't exist");	
	    }
	    #Default filter name is the base of the filename
	    my $filter_name = basename($filter_file);
	
	    my %filter_list = ();
	    open(FILE,"$filter_file") || die "Can't open file $filter_file\n";
	    while (<FILE>) {
	    	chomp;
	    	my ($entry) = $_ =~ /^(\S+)/;
	    	next unless $entry;
	    	$filter_list{uc($entry)} = 1;
	    }
	    $filter_files{$filter_name} = \%filter_list;
		$count++;
	}
	$self->{filter_files} = \%filter_files;
	return $self->{filter_files};
}


#Add the organism
sub organism {
	 my ($self, $organism) = @_;

    if (defined $organism) {
		$self->{'organism'} = $organism;
    } elsif (! defined $self->{'organism'}) {
		modules::Exception->throw("organism not set");
    }

    return $self->{'organism'};
}

#Add the genemapper
sub gene_mapper {
	 my ($self, $gene_mapper) = @_;

    if (defined $gene_mapper) {
		$self->{'gene_mapper'} = $gene_mapper;
    } elsif (! defined $self->{'gene_mapper'}) {
		modules::Exception->throw("genemapper not set");
    }

    return $self->{'gene_mapper'};
}

#Map all the headers
sub headers {
	my $self = shift;
	my %header_col_lookup;
	
#						
#    foreach my $filter (@{$self->{filter_list_order}}) {
#		push @filter_names, $filter->name;
#    }

	my @filter_names = @{$self->{filter_list_order}};
	my @annotation_headers = @{$self->{annotation_headers}};
	my @run_headers = @{$self->{run_headers}};
	my @extra_headers = @{$self->{extra_headers}};


	my @common_headers = (
						  @run_headers,			      
					      @filter_names,
					      @extra_headers,
					      @annotation_headers
						  );
	

	for (my $i = 0; $i < scalar @common_headers; $i++){
    	$header_col_lookup{$common_headers[$i]} = $i;
    	
	}
	$self->{headers} = \@common_headers;
	$self->{header_lookup} = \%header_col_lookup;
}

#Get the runid
sub run {
    my ($self, $run) = @_;

    if (defined $run) {
		$self->{'run'} = $run;
    } elsif (! defined $self->{'run'}) {
		modules::Exception->throw("run not set");
    }

    return $self->{'run'};
}

#Get the confdir directory for whole genome cases 
sub confdir {
    my ($self, $confdir) = @_;

    if (defined $confdir) {
		$self->{'confdir'} = $confdir;
    } elsif (! defined $self->{'confdir'}) {
		modules::Exception->throw("confdir not set");
    }

    return $self->{'confdir'};
}

#Get the sample_typeid
sub source_type {
    my ($self, $source_type) = @_;

    if (defined $source_type) {
		$self->{'source_type'} = $source_type;
    } elsif (! defined $self->{'source_type'}) {
		modules::Exception->throw("source_type not set");
    }

    return $self->{'source_type'};
}



#Get the sample_typeid
sub sample_type {
    my ($self, $sample_type) = @_;

    if (defined $sample_type) {
		$self->{'sample_type'} = $sample_type;
    } elsif (! defined $self->{'sample_type'}) {
		modules::Exception->throw("sample_type not set");
    }

    return $self->{'sample_type'};
}

#Get the sample_name
sub sample_name {
    my ($self, $sample_name) = @_;

    if (defined $sample_name) {
		$self->{'sample_name'} = $sample_name;
    } elsif (! defined $self->{'sample_name'}) {
		modules::Exception->throw("sample_name not set");
    }

    return $self->{'sample_name'};
}

#Get the gene_col_name
sub gene_col_name {
	my ($self, $gene_col_name) = @_;

    if (defined $gene_col_name) {
		$self->{'gene_col_name'} = $gene_col_name;
    } elsif (! defined $self->{'gene_col_name'}) {
		modules::Exception->throw("gene_col_name not set");
    }

    return $self->{'gene_col_name'};
}

#Generate the file containing the polyphen info
sub generate_polyphen_input {
	my $self = shift;
	my %args = @_;
		
	my %sorted_pass_rows = %{$self->{sorted_pass_rows}};
	#Now print the appropriate fields
	my $divider_count = @{$self->{headers}} - 1;

	my $count = 0;

	#Need this for polyphen	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $ref_base_col_num = $header_col_lookup{'ref_allele'};
	my $var_base_col_num = $header_col_lookup{'var_allele'};
	my $aa_change_col_num = $header_col_lookup{'aa_change'};
	my $uniprot_col_num = $header_col_lookup{'uniprot'};
	
	my %polyphen_data = ();
	

	for my $type (keys %sorted_pass_rows) {
		
		for my $splice ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}} ) {
		    for my $depth ( sort {$b<=>$a} keys %{$sorted_pass_rows{$type}{$splice}} ) {
		    	for my $chr ( sort keys %{$sorted_pass_rows{$type}{$splice}{$depth}} ) {
		    		for my $coord ( sort {$a<=>$b} keys %{$sorted_pass_rows{$type}{$splice}{$depth}{$chr}} ) {
		    			$polyphen_data{$chr}{$coord}{ref} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$ref_base_col_num];
		    			$polyphen_data{$chr}{$coord}{var} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$var_base_col_num];
						$polyphen_data{$chr}{$coord}{aachange} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$aa_change_col_num];
						$polyphen_data{$chr}{$coord}{uniprot} = $sorted_pass_rows{$type}{$splice}{$depth}{$chr}{$coord}->[$uniprot_col_num];
						$count++;
		    		}
				}
			}
		}
	}
    return (\%polyphen_data,$count);
}


#Print the final TSV that goes to the spreadsheet and the match and no match files
sub print_to_files {
    my ($self) = shift;
    my %args = @_;
    
    
    my @required_args = (
    					'-tsv_file',
    					'-pass_file'
			 			);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $output_filename = $args{-tsv_file};
    open(my $OUT, ">$output_filename")
	or modules::Exception->throw("Unable to open output file [$output_filename]");

	my $pass_file = $args{-pass_file};
	open(my $PASS, ">$pass_file")
    or modules::Exception->throw("Unable to open output file [$pass_file]");

	(my $fail_file = $pass_file) =~ s/match/nomatch/;

	open(my $FAIL, ">$fail_file") 
	or modules::Exception->throw("Unable to open output file [$fail_file]");


	# print the headers
    print $OUT join("\t",@{$self->{headers}}) . "\n";
	
	my %sorted_pass_rows = %{$self->{sorted_pass_rows}};
	my %rare_allele_rows = %{$self->{sorted_allele_rows}};
	my %no_freq_rows = %{$self->{sorted_no_freq_rows}};
	
	
	#Now print the appropriate fields
	my $divider_count = @{$self->{headers}} - 1;

	my %header_col_lookup = %{$self->{header_lookup}};
	
	my $variant_type = $self->{variant_type};

	if ($variant_type eq 'snv') {

		
		for my $row_type (qw(PASS RARE_ALLELE NO_FREQ FAIL)) {
			my %row_data;
			if ($row_type eq 'FAIL') {
				my $final_divider = "LOW_PRIORITY LINE\t" . "---\t" x $divider_count;
				chop($final_divider);
				print $OUT "$final_divider\n";
			    
			
				my @fail_rows = @{$self->{fail_rows}};
				for my $field_array ( @fail_rows) {
				    my $line = join("\t",@{$field_array});
				    print $OUT $line,"\n";
				 	my $rest;
				 	if ($self->{paired_variants}) {
				 		$rest = join("^^^",  'LOW_PRIORITY',
				 							$field_array->[$header_col_lookup{ref_allele}]."->".$field_array->[$header_col_lookup{var_allele}],
				 	 						'ref_count:'.$field_array->[$header_col_lookup{ref_allele_count}],
				 	 						'var_count:'.$field_array->[$header_col_lookup{var_allele_count}],
				 	 						'snv_score:'.$field_array->[$header_col_lookup{snv_score}],
				 	 						'clr_score:'.$field_array->[$header_col_lookup{clr_score}],
				 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}],
				 	 						'normal_alleles:'.$field_array->[$header_col_lookup{normal_alleles}],
				 	 						'snv_class:'.$field_array->[$header_col_lookup{snv_class}]
				 	 						
				 	 						);
				 	} else {
				 		$rest = join("^^^",  'LOW_PRIORITY',
				 							$field_array->[$header_col_lookup{ref_allele}]."->".$field_array->[$header_col_lookup{var_allele}],
				 	 						'ref_count:'.$field_array->[$header_col_lookup{ref_allele_count}],
				 	 						'var_count:'.$field_array->[$header_col_lookup{var_allele_count}],
				 	 						'snv_score:'.$field_array->[$header_col_lookup{snv_score}],
				 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}]
				 	 						);
				 	}
				 	$rest =~ s/ /_/g;
				 	print $FAIL join("\t",$field_array->[$header_col_lookup{chr}], $field_array->[$header_col_lookup{coord}], $field_array->[$header_col_lookup{coord}], $rest) ."\n";
				}
			} else {
				my $row_print;
				if ($row_type eq 'RARE_ALLELE') {
					%row_data = %rare_allele_rows;
					$row_print = 'RARE';
				} elsif ($row_type eq 'NO_FREQ') {
					%row_data = %no_freq_rows;
					$row_print = 'NO_FREQ';
				} else {
					%row_data = %sorted_pass_rows;
					$row_print = 'NOVEL';
				}
				
				for my $type (keys %row_data) {
					my $divider = "$type $row_print\t" . "---\t" x $divider_count;
					chop $divider;
					print $OUT "$divider\n";
					
					for my $splice ( sort {$a<=>$b} keys %{$row_data{$type}} ) {
					    for my $depth ( sort {$b<=>$a} keys %{$row_data{$type}{$splice}} ) {
					    	for my $chr ( sort keys %{$row_data{$type}{$splice}{$depth}} ) {
					    		for my $coord ( sort {$a<=>$b} keys %{$row_data{$type}{$splice}{$depth}{$chr}} ) {
					    			print $OUT join("\t",@{$row_data{$type}{$splice}{$depth}{$chr}{$coord}}) . "\n";
					    			#Make a useful key for overlap file  
					    			
					    			my $rest;
					    			if ($self->{paired_variants}) {
					    				$rest = join("^^^",$row_print,
					    								$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele}] .'->'. $row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele}],
					    								'ref_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele_count}],
				 	 									'var_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele_count}],
				 	 									'snv_score:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{snv_score}],
				 	 									'clr_score:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{clr_score}],
				 	 									'final_status:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{final_status}],
				 	 									'normal_alleles:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{normal_alleles}],
				 	 									'snv_class:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{snv_class}],	
					    								);
					    			} else {
					    				$rest = join("^^^",$row_print,
					    								$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele}] .'->'. $row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele}],
					    								'ref_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{ref_allele_count}],
				 	 									'var_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{var_allele_count}],
				 	 									'snv_score:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{snv_score}],
				 	 									'final_status:'.$row_data{$type}{$splice}{$depth}{$chr}{$coord}->[$header_col_lookup{final_status}]
					    								);
					    			}
					    			$rest =~ s/ /_/g;
					    			
					    			print $PASS join("\t",
					    							$chr,
					    							$coord,
					    							$coord,
					    							$rest
					    							)."\n";
								}
							}
						}
					}
				}
		}
		


		}
	
	} elsif ($variant_type eq 'indel') {
		
		for my $row_type (qw(PASS RARE_ALLELE NO_FREQ FAIL)) {
			my %row_data;
			if ($row_type eq 'FAIL') {
				my $final_divider = "LOW_PRIORITY LINE\t" . "---\t" x $divider_count;
				chop($final_divider);
				print $OUT "$final_divider\n";
						
				my @fail_rows = @{$self->{fail_rows}};
				for my $field_array ( @fail_rows) {
					my $line = join("\t",@{$field_array});
				    print $OUT $line,"\n";
				    my $rest;
					if ($self->{paired_variants}) {
						$rest = join("^^^",  'LOW_PRIORITY',
											$field_array->[$header_col_lookup{var_type}],
											$field_array->[$header_col_lookup{ref_allele}]."->".$field_array->[$header_col_lookup{var_allele}],
				 	 						'ref_count:'.$field_array->[$header_col_lookup{ref_allele_count}],
				 	 						'var_count:'.$field_array->[$header_col_lookup{var_allele_count}],
				 	 						'snv_score:'.$field_array->[$header_col_lookup{var_score}],
				 	 						'clr_score:'.$field_array->[$header_col_lookup{clr_score}],
				 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}],
				 	 						'normal_alleles:'.$field_array->[$header_col_lookup{normal_alleles}],
				 	 						'variant_class:'.$field_array->[$header_col_lookup{variant_class}]
				 	 						
										);
					} else {
						$rest = join("^^^",  'LOW_PRIORITY',
											$field_array->[$header_col_lookup{var_type}],
											$field_array->[$header_col_lookup{ref_allele}]."->".$field_array->[$header_col_lookup{var_allele}],
				 	 						'ref_count:'.$field_array->[$header_col_lookup{ref_allele_count}],
				 	 						'var_count:'.$field_array->[$header_col_lookup{var_allele_count}],
				 	 						'snv_score:'.$field_array->[$header_col_lookup{var_score}],
				 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}]
										);
					}
					$rest =~ s/ /_/g;
					print $FAIL join("\t",$field_array->[$header_col_lookup{chr}], $field_array->[$header_col_lookup{start_coord}], $field_array->[$header_col_lookup{end_coord}], $rest) ."\n";
				}
			} else {
				my $row_print;
				if ($row_type eq 'RARE_ALLELE') {
					%row_data = %rare_allele_rows;
					$row_print = 'RARE';
				} elsif ($row_type eq 'NO_FREQ') {
					%row_data = %no_freq_rows;
					$row_print = 'NO_FREQ';
				} else {
					%row_data = %sorted_pass_rows;
					$row_print = 'NOVEL';
				}
				
				for my $type (keys %row_data) {
					my $divider = "$type $row_print\t" . "---\t" x $divider_count;
					chop $divider;
					print $OUT "$divider\n";
					
					for my $splice ( sort {$a<=>$b} keys %{$row_data{$type}} ) {
					    for my $depth ( sort {$b cmp $a} keys %{$row_data{$type}{$splice}} ) {
					    	for my $chr ( sort keys %{$row_data{$type}{$splice}{$depth}} ) {
					    		for my $start_coord ( sort {$a<=>$b} keys %{$row_data{$type}{$splice}{$depth}{$chr}} ) {
					    			for my $end_coord ( sort {$a<=>$b} keys %{$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}} ) {
					    				for my $var_base ( sort keys %{$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}} ) {
					    			
							    			print $OUT join("\t",@{$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}}) . "\n";
							    			#Make a useful key for overlap file  
							    			my $rest;
					    					if ($self->{paired_variants}) {
							    				$rest = join("^^^",$row_print,
							    									$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_type}],
																	$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele}]."->".$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele}],
				 	 												'ref_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele_count}],
				 	 												'var_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele_count}],
				 	 												'snv_score:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_score}],
				 	 												'clr_score:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{clr_score}],
				 	 												'final_status:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{final_status}],
				 	 												'normal_alleles:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{normal_alleles}],
							    									'variant_class:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{variant_class}]
							    								
							    								);
					    					} else {
					    						$rest = join("^^^",$row_print,
							    									$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_type}],
																	$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele}]."->".$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele}],
				 	 												'ref_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{ref_allele_count}],
				 	 												'var_count:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_allele_count}],
				 	 												'snv_score:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{var_score}],
				 	 												'final_status:'.$row_data{$type}{$splice}{$depth}{$chr}{$start_coord}{$end_coord}{$var_base}->[$header_col_lookup{final_status}]
							    								);
					    					}
							    			$rest =~ s/ /_/g;
							    			print $PASS join("\t",
							    							$chr,
							    							$start_coord,
							    							$end_coord,
							    							$rest
							    							)."\n";
					    				}
					    			}
								}
							}
						}
					}
				}
			
			}
		}
		

	} elsif ($variant_type eq 'sv') {
		my $sv_caller = $self->{sv_caller};
		for my $row_type (qw(PASS RARE_ALLELE NO_FREQ FAIL)) {
			my %row_data;
			if ($row_type eq 'FAIL') {
				my $final_divider = "LOW_PRIORITY LINE\t" . "---\t" x $divider_count;
				chop($final_divider);
				print $OUT "\n$final_divider\n\n";
			    
			
				my @fail_rows = @{$self->{fail_rows}};
				for my $field_array ( @fail_rows) {
				    my $line = join("\t",@{$field_array});
				    print $OUT $line,"\n";				 	
				 	my $rest = join("^^^",  'LOW_PRIORITY',
				 							$sv_caller,
				 							$field_array->[$header_col_lookup{sv_type}],
				 	 						'sv_id:'.$field_array->[$header_col_lookup{sv_id}],
				 	 						'final_status:'.$field_array->[$header_col_lookup{final_status}]
				 	 						);
				 	$rest =~ s/ /_/g;
				 	print $FAIL join("\t",$field_array->[$header_col_lookup{chr_1}], $field_array->[$header_col_lookup{coord1}], $field_array->[$header_col_lookup{chr_2}], $field_array->[$header_col_lookup{coord2}], $rest) ."\n";
				}
			} else {
				my $row_print;
				if ($row_type eq 'RARE_ALLELE') {
					%row_data = %rare_allele_rows;
					$row_print = 'RARE';
				} elsif ($row_type eq 'NO_FREQ') {
					%row_data = %no_freq_rows;
					$row_print = 'NO_FREQ';
				} else {
					%row_data = %sorted_pass_rows;
					$row_print = 'NOVEL';
				}
				
				for my $sv_type (keys %{$row_data{$sv_caller}}) {
					my $divider = "$sv_type $row_print\t" . "---\t" x $divider_count;
					chop $divider;
					print $OUT "\n$divider\n\n";
					
			    	for my $chr1 ( sort keys %{$row_data{$sv_caller}{$sv_type}} ) {
			    		for my $chr2 ( sort keys %{$row_data{$sv_caller}{$sv_type}{$chr1}} ) {
				    		for my $start_coord ( sort {$a<=>$b} keys %{$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}} ) {
				    			for my $end_coord ( sort {$a<=>$b} keys %{$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}{$start_coord}} ) {
				    				for my $gene (keys %{$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}{$start_coord}{$end_coord}}) {
						    			print $OUT join("\t",@{$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}{$start_coord}{$end_coord}{$gene}}) . "\n";
						    			#Make a useful key for overlap file  
						    			
						    			
						    			my $rest = join("^^^",$row_print,
						    								$sv_caller,
						    								$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}{$start_coord}{$end_coord}{$gene}->[$header_col_lookup{sv_type}],
						    								$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}{$start_coord}{$end_coord}{$gene}->[$header_col_lookup{sv_id}],
							 								$row_data{$sv_caller}{$sv_type}{$chr1}{$chr2}{$start_coord}{$end_coord}{$gene}->[$header_col_lookup{final_status}],
						    								);
						    			
						    			$rest =~ s/ /_/g;
						    			
						    			print $PASS join("\t",
						    							$chr1,
						    							$chr2,
						    							$start_coord,
						    							$end_coord,
						    							$rest
						    							)."\n";
				    					
				    				}
								}
				    		}
						}
			    	}
				}
			}
		}
	}
	
	close($OUT);
    return 1;
}

#Generate the pass fail status
sub generate_pass_fail {
	my $self = shift;
	my %args = @_;
	my $debug = defined $args{-debug}?1:0;

	my $variant_type = $self->{variant_type};
	my @lines = ();
	
	if ($variant_type eq 'snv') {
	    # Get the snvs to report on
	    my %snv_ids = $self->_get_variants($args{-chr},$args{-start},$args{-end});

		#Get the formatted lines
		if (defined $args{-chr}) {
			@lines = $self->_generate_snv_lines(-snv_ids=>\%snv_ids,-chr=>$args{-chr});		
		} else {
			@lines = $self->_generate_snv_lines(-snv_ids=>\%snv_ids);
		}
	} elsif ($variant_type eq 'sv') {
		# Get the snvs to report on
	    my %sv_ids = $self->_get_variants($args{-chr},$args{-start},$args{-end});

		#Get the formatted lines
		if (defined $args{-chr}) {
			@lines = $self->_generate_sv_lines(-sv_ids=>\%sv_ids,-chr=>$args{-chr});		
		} else {
			@lines = $self->_generate_sv_lines(-sv_ids=>\%sv_ids);
		}
	} else {
		my %indel_ids = $self->_get_variants($args{-chr},$args{-start},$args{-end});
		if (defined $args{-chr}) {
			@lines = $self->_generate_indel_lines(-indel_ids=>\%indel_ids,-chr=>$args{-chr});
		} else {
			@lines = $self->_generate_indel_lines(-indel_ids=>\%indel_ids);
		}
	}
	$self->_sort_pass_fail(-lines=>\@lines);
}

#sort the lines into pass/fail and put into data structures
sub _sort_pass_fail {
	my $self = shift;
	my %args = @_;
	
	my @required_args = (
						'-lines'
			 			);
    foreach my $required_arg (@required_args){
		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    	
    my $variant_type = $self->{variant_type};
    
    my @pass_rows;
    my @rare_allele_rows;
    my @no_freq_rows;
    my %fail_rows;
	my @lines = @{$args{-lines}};
	my %filter_set = %{$self->{filter_set}};
	my %report_filters = %{$self->{report_set}};
	my %header_col_lookup = %{$self->{header_lookup}};
	my %filter_files = keys %{$self->{filter_files}}?%{$self->{filter_files}}:();
	
	#Keep track of genes filtered out for reporting			
	my %filtered_genes;			
	
	my $dbsnp_column_number;
	if ($variant_type eq 'snv') {
		$dbsnp_column_number = $header_col_lookup{filter_dbsnp_snv};
	} elsif ($variant_type eq 'sv') {
		$dbsnp_column_number = $header_col_lookup{filter_dgv};
	} else {
		$dbsnp_column_number = $header_col_lookup{filter_dbsnp_indel};
	}
	
	my $common_column_number;
	if ($self->organism eq 'mouse') {
		if ($variant_type eq 'snv') {
			$common_column_number = $header_col_lookup{filter_common_snv};
		} else {
			$common_column_number = $header_col_lookup{filter_common_indel};
		}
	}
	
	my %filter_count = ();
	my $line_count = @lines;
	
	if ($variant_type eq 'snv') {
		$filter_count{filter_snv}{pass} = $line_count; 	
	} elsif ($variant_type eq 'sv') {
		$filter_count{filter_sv}{pass} = $line_count;
	} else {
		$filter_count{filter_indel}{pass} = $line_count;
	}
	
ROW:	
	for my $line (@lines) {

		my $dbsnp_fail = 0;
		chomp $line;
		my @cols = split ("\t",$line);




		my $gene;
		if (!defined $cols[$header_col_lookup{$self->{gene_col_name}}]) {
			#Cosmic overlaps don't need to overlap genes so don't report this warning
			#modules::Exception->warning("No gene for snv entry with line ($line)");
		} else {
			$gene = $cols[$header_col_lookup{$self->{gene_col_name}}];
		}
	
		#First get the filter pass counts
		for my $report_filter ( keys %report_filters ) {
		    #Keep a tally of the number of records passing each filter
			if (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] eq 'PASS') {
				$filter_count{$report_filter}{pass}++;
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] eq 'FAIL') {
				$filter_count{$report_filter}{fail}++;
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] eq 'N/A') {
		    	next; #N/A for sv type tra and ins
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] eq '""') {
		    	next; #"" for sv type tra and ins
		    }elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] == 1) {
				$filter_count{$report_filter}{pass}++;
		    } elsif (defined  $header_col_lookup{$report_filter} && $cols[$header_col_lookup{$report_filter}] == 0) {
				$filter_count{$report_filter}{fail}++;
		    }
		}
		
		if ($gene) {
			#Handle any additional file filters of gene names
			for my $file_filter ( keys %filter_files ) {
				if (exists $filter_files{$file_filter}->{$gene}) {
					$filtered_genes{$gene}++;
					$filter_count{$file_filter}{fail}++;
				} else {
					$filter_count{$file_filter}{pass}++;
				}
			}
		}	
		
		if ($cols[$dbsnp_column_number] eq 'PASS') {
			$cols[$dbsnp_column_number] = 'NOVEL';
		}
		
		if ($common_column_number) {
			if ($cols[$common_column_number] eq 'PASS') {
				$cols[$common_column_number] = 'NOT_SEEN';
			} elsif ($cols[$common_column_number] eq 'FAIL') {
				$cols[$common_column_number] = 'PREVIOUSLY_SEEN';
			}		
		}
		
		my %failed_reasons = ();
		#my $local_line = join("\t",@cols);
		
		#Entries with a cutoff
	    foreach my $filter (keys %filter_set) {
			
			#Handle the special cases first (currently only snvs)
			#These are all snv specific filters
			if ($filter_set{$filter} =~ /allele_freq=([0-9\.]+)/) {
				my $cutoff = $1;
				my $allele_freq_str;
				if ($variant_type eq 'sv') {
					$allele_freq_str = $cols[$header_col_lookup{dgv_freq}];
				} else {
					$allele_freq_str = $cols[$header_col_lookup{dbsnp_match}];
				}
				
				#special filtering for dbsnp human where rare alleles are still a pass
				if ($cols[$dbsnp_column_number] eq 'NOVEL' || $cols[$dbsnp_column_number] eq 'N/A') {
					#Skip if we're already passed or don't have the data
					next;
				} elsif ($allele_freq_str =~ /:R([0-9\.]+):V([0-9\.]+)/) {
					my $dbsnp_ref_allele_freq = $1;
					my $dbsnp_var_allele_freq = $2;
					
					if ($dbsnp_var_allele_freq < $cutoff) {
						$cols[$dbsnp_column_number] = 'RARE_ALLELE';
					} elsif ($dbsnp_ref_allele_freq < $cutoff) {
						$cols[$dbsnp_column_number] = 'RARE_REF';
						$failed_reasons{filter_dbsnp}++;
					} else {
						$failed_reasons{filter_dbsnp}++;
					}
				} elsif ($allele_freq_str =~ /[0-9]+/ && $variant_type eq 'sv') {
					#Check there is a pop freq (not all types have this)
					my $dgv_freq = $cols[$header_col_lookup{dgv_freq}];
					if ($dgv_freq < $cutoff) {
						$cols[$dbsnp_column_number] = 'RARE_ALLELE';
					} elsif ((1-$dgv_freq) < $cutoff) {
						$cols[$dbsnp_column_number] = 'RARE_REF';
						$failed_reasons{filter_dgv}++;
					} else {
						$failed_reasons{filter_dgv}++;
					}
					#exit;
				} elsif ($allele_freq_str =~ /No/) {
					$cols[$dbsnp_column_number] = 'NO_FREQ';
					$failed_reasons{filter_dbsnp}++;
				} else {
					if ($variant_type eq 'sv') {
						$failed_reasons{filter_dgv}++;
					} else {
						$failed_reasons{filter_dbsnp}++;
					}
				}
								
				
			} elsif ($filter eq 'base_change') {
				#filter out base changes we don't want
				my @changes = split(",",$filter_set{$filter});
				my $current_change = $cols[$header_col_lookup{ref_allele}] . '->' . $cols[$header_col_lookup{var_allele}];
				my $matched = 0;
				
				for my $change (@changes) {
					if ($change =~ /^\!/) {
						$change =~ s/^\!//;
						#If it matches a change we're filtering out
						if ($current_change eq $change) {
							$failed_reasons{base_change}++;
						}
					} else {
						modules::Exception->throw("ERROR: Can only filter out base changes");
					}
				}
			}  elsif ($filter eq 'supporting_reads') {
				my $min_total = $filter_set{$filter};
				my $total_supporting = $cols[$header_col_lookup{split_reads}] + $cols[$header_col_lookup{paired_reads}];
				if ($total_supporting < $min_total) {
					$failed_reasons{'Not enough supporting reads'}++;
				}
			} elsif (!exists $header_col_lookup{$filter}) {
				#Catchall; shouldn't get here as all the special cases should have been handled already
				modules::Exception->throw("ERROR: Problem with pass/fail analysis for filter $filter");
			} elsif ($cols[$header_col_lookup{$filter}] =~ /FAIL/) {
				#If it's a pass/fail field (eg filter_dbsnp_snv)
				$failed_reasons{$filter}++;
			} elsif ($cols[$header_col_lookup{$filter}] =~ /PREVIOUSLY_SEEN/) {
				#If it's filter_common field
				$failed_reasons{$filter}++;
			} elsif ($cols[$header_col_lookup{$filter}] =~ /^\d/ && $cols[$header_col_lookup{$filter}] < $filter_set{$filter}) {
				#If it's a value cutoff (eg read_depth, median_quality_score)
				my $cutoff_fail = "$filter < $filter_set{$filter}";
				$failed_reasons{$cutoff_fail}++;
			} 
	    }
	
		if (keys %failed_reasons) {
			my $reason_str = join(" , ",keys %failed_reasons);
			my $local_line = join("\t",@cols);
			$local_line =~ s/OVERALL_PASS/FAIL REASONS: $reason_str/;
			$fail_rows{$local_line} = 1;
			next ROW;
		}
	
	
		if ($gene) {
			#Additional file cutoffs
			for my $file_filter ( keys %filter_files ) {
				if (exists $filter_files{$file_filter}->{$gene}) {
					$filter_count{$file_filter}{fail}++;
					my $local_line = join("\t",@cols);
					$local_line =~ s/OVERALL_PASS/problem_gene FAIL/;
					$fail_rows{$local_line}  = 1;
					next ROW;
				}
			}  
		}
		
	    my $or_pass = 0;

		#If this field is either NON-SYN or SPLICE then the row is passed
		if ($variant_type eq 'snv' && $cols[$header_col_lookup{snv_exon_type}] ne 'SYN') {
			$or_pass = 1;
		} elsif ($variant_type eq 'indel' || $variant_type eq 'sv') {
			$or_pass = 1;
		}
	
	    unless ($or_pass) {
	    	my $local_line = join("\t",@cols);
			$local_line =~ s/OVERALL_PASS/synonomous FAIL/;
	    	$fail_rows{$local_line} = 1;
			next ROW;
	    }
	    
	    #Finally determine whether it's rare allele or full pass
	    if ($cols[$dbsnp_column_number] eq 'RARE_ALLELE') {
			push @rare_allele_rows, \@cols;    	
	    } elsif ($cols[$dbsnp_column_number] eq 'NO_FREQ') {
	    	push @no_freq_rows, \@cols;
	    } else {
		    push @pass_rows, \@cols;
	    }
	}
		
	#Now generate the data structure used for custom sorting the passed rows
	my %sorted_pass_rows;
	my %sorted_allele_rows;
	my %sorted_no_freq_rows;
	
	if ($variant_type eq 'snv') {
	
		foreach my $pass_row (@pass_rows) {
			my $mutant_depth = $pass_row->[$header_col_lookup{read_depth}];
			
			if ($pass_row->[$header_col_lookup{snv_exon_type}] =~ /SPLICE/) {
		    	$sorted_pass_rows{'SNV'}{1}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{coord}]} = $pass_row;
		    } else {
		    	$sorted_pass_rows{'SNV'}{0}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{coord}]} = $pass_row;
		    }
		}
		
		foreach my $allele_row (@rare_allele_rows) {
			my $mutant_depth = $allele_row->[$header_col_lookup{read_depth}];
			
			if ($allele_row->[$header_col_lookup{snv_exon_type}] =~  /SPLICE/) {
		    	$sorted_allele_rows{'SNV'}{1}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{coord}]} = $allele_row;
		    } else {
				$sorted_allele_rows{'SNV'}{0}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{coord}]} = $allele_row;
		    }
		}
		
		foreach my $no_freq_row (@no_freq_rows) {
			my $mutant_depth = $no_freq_row->[$header_col_lookup{read_depth}];
			
			if ($no_freq_row->[$header_col_lookup{snv_exon_type}] =~  /SPLICE/) {
		    	$sorted_no_freq_rows{'SNV'}{1}{$mutant_depth}{$no_freq_row->[$header_col_lookup{chr}]}{$no_freq_row->[$header_col_lookup{coord}]} = $no_freq_row;
		    } else {
				$sorted_no_freq_rows{'SNV'}{0}{$mutant_depth}{$no_freq_row->[$header_col_lookup{chr}]}{$no_freq_row->[$header_col_lookup{coord}]} = $no_freq_row;
		    }
		}
		
		
	} elsif ($variant_type eq 'indel') {
		foreach my $pass_row (@pass_rows) {
			my $indel_type = $pass_row->[$header_col_lookup{var_type}];
			my $mutant_depth = $pass_row->[$header_col_lookup{read_depth}];
			
			if ($pass_row->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
		    	$sorted_pass_rows{$indel_type}{1}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{start_coord}]}{$pass_row->[$header_col_lookup{end_coord}]}{$pass_row->[$header_col_lookup{var_allele}]} = $pass_row;
		    } else {
		    	$sorted_pass_rows{$indel_type}{0}{$mutant_depth}{$pass_row->[$header_col_lookup{chr}]}{$pass_row->[$header_col_lookup{start_coord}]}{$pass_row->[$header_col_lookup{end_coord}]}{$pass_row->[$header_col_lookup{var_allele}]} = $pass_row;
		    }
		}
		
		foreach my $allele_row (@rare_allele_rows) {
			my $indel_type = $allele_row->[$header_col_lookup{var_type}];
			my $mutant_depth = $allele_row->[$header_col_lookup{read_depth}];
			
			if ($allele_row->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
		    	$sorted_allele_rows{$indel_type}{1}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{start_coord}]}{$allele_row->[$header_col_lookup{end_coord}]}{$allele_row->[$header_col_lookup{var_allele}]} = $allele_row;
		    } else {
		    	$sorted_allele_rows{$indel_type}{0}{$mutant_depth}{$allele_row->[$header_col_lookup{chr}]}{$allele_row->[$header_col_lookup{start_coord}]}{$allele_row->[$header_col_lookup{end_coord}]}{$allele_row->[$header_col_lookup{var_allele}]} = $allele_row;
		    }
		}

		foreach my $no_freq_rows (@no_freq_rows) {
			my $indel_type = $no_freq_rows->[$header_col_lookup{var_type}];
			my $mutant_depth = $no_freq_rows->[$header_col_lookup{read_depth}];
			
			if ($no_freq_rows->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
		    	$sorted_no_freq_rows{$indel_type}{1}{$mutant_depth}{$no_freq_rows->[$header_col_lookup{chr}]}{$no_freq_rows->[$header_col_lookup{start_coord}]}{$no_freq_rows->[$header_col_lookup{end_coord}]}{$no_freq_rows->[$header_col_lookup{var_allele}]} = $no_freq_rows;
		    } else {
		    	$sorted_no_freq_rows{$indel_type}{0}{$mutant_depth}{$no_freq_rows->[$header_col_lookup{chr}]}{$no_freq_rows->[$header_col_lookup{start_coord}]}{$no_freq_rows->[$header_col_lookup{end_coord}]}{$no_freq_rows->[$header_col_lookup{var_allele}]} = $no_freq_rows;
		    }
		}
		
	} elsif ($variant_type eq 'sv') {
		my $sv_caller = $self->{sv_caller};
		foreach my $pass_row (@pass_rows) {
			my $sv_type = $pass_row->[$header_col_lookup{sv_type}];
		    $sorted_pass_rows{$sv_caller}{$sv_type}{$pass_row->[$header_col_lookup{chr_1}]}{$pass_row->[$header_col_lookup{chr_2}]}{$pass_row->[$header_col_lookup{coord1}]}{$pass_row->[$header_col_lookup{coord2}]}{$pass_row->[$header_col_lookup{ensembl}]} = $pass_row;
		}
		
		foreach my $allele_row (@rare_allele_rows) {
			my $sv_type = $allele_row->[$header_col_lookup{sv_type}];
		    $sorted_pass_rows{$sv_caller}{$sv_type}{$allele_row->[$header_col_lookup{chr_1}]}{$allele_row->[$header_col_lookup{chr_2}]}{$allele_row->[$header_col_lookup{coord1}]}{$allele_row->[$header_col_lookup{coord2}]}{$allele_row->[$header_col_lookup{ensembl}]} = $allele_row;
		}
		foreach my $no_freq_rows (@no_freq_rows) {
			my $sv_type = $no_freq_rows->[$header_col_lookup{sv_type}];
		    $sorted_no_freq_rows{$sv_caller}{$sv_type}{$no_freq_rows->[$header_col_lookup{chr_1}]}{$no_freq_rows->[$header_col_lookup{chr_2}]}{$no_freq_rows->[$header_col_lookup{coord1}]}{$no_freq_rows->[$header_col_lookup{coord2}]}{$no_freq_rows->[$header_col_lookup{ensembl}]} = $no_freq_rows;
		}
	}
	
	#print Dumper \%fail_rows;
	#exit;
	my @fail_rows = ();
	for my $fail_row (sort {my @afields = split("\t",$a); my @bfields = split ("\t",$b); $afields[0] cmp $bfields[0] || $afields[1] <=> $bfields[1]} keys %fail_rows) {
		my @cols = split ("\t",$fail_row);
		push @fail_rows, \@cols;
	}
	
	my @all_rows = (@pass_rows,@rare_allele_rows,@fail_rows);
	$self->{all_rows} = \@all_rows;
	$self->{pass_rows} = \@pass_rows;
	$self->{allele_rows} = \@rare_allele_rows;
	$self->{no_freq_rows} = \@no_freq_rows;
	$self->{fail_rows} = \@fail_rows;
	$self->{filter_count} = \%filter_count;
	$self->{sorted_pass_rows} = \%sorted_pass_rows;
	$self->{sorted_allele_rows} = \%sorted_allele_rows;
	$self->{sorted_no_freq_rows} = \%sorted_no_freq_rows;
	$self->{filtered_genes} = \%filtered_genes;
}

#Get the formatted indel lines
sub _generate_indel_lines {
	my $self = shift;
	my %args = @_;

	my @required_args = (
						'-indel_ids'
			 			);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $indel_ids = $args{-indel_ids};
	my $gene_mapper = $self->{gene_mapper};
	
	#Needed for indel_row entries
	my %indel_lookup = ();
	
	my @filter_names;
	if ($self->{genome}) {
		@filter_names = @{$self->{filter_info}};
	} else {
	    foreach my $filter (@{$self->{filter_info}}){
			push @filter_names, $filter->name;
	    }
		
	}

	my $debug_chr;
	if (defined $args{-chr}) {
		$debug_chr = $args{'-chr'};
	}


	my @lines = ();
	
	my %chr_xml_data;
	my $current_chr = -1;
	
	foreach my $indel_id (sort {$a<=>$b} keys %{$indel_ids}) {
		
		
		my $indel;
		my $filt_indel;
		my $var_base = my $start_coord = my $end_coord = my $xml_var_base = my $chr;
		
		if ($self->{genome}) {
			($chr,$start_coord,$end_coord,$var_base) = split(':',$indel_ids->{$indel_id});
			$xml_var_base = $var_base;
			$xml_var_base =~ s/^\+/i/;
			$xml_var_base =~ s/^\-/d/;
		} else {
			($indel) = modules::Adaptors::Variant->search(id => $indel_id);
			
			$filt_indel
			    = modules::FilteredVariant->new('variant' => $indel,
							    'filter_list_order' => $self->{filter_info});
			
			$var_base = $indel->var_base;
			$start_coord = $indel->start_coord;
			$end_coord = $indel->end_coord;		
			$chr = $indel->chr;						
		}
		
		
		if ($debug_chr) {
			next unless $chr eq $debug_chr;
		}
		
		
		if ($self->{genome}) {
			#Here need to update the db_var xml for the new chromosome
			if ($current_chr ne $chr) {
				
				my $dbvar_file_name = $self->{confdir} . '/' . $self->{sample_name} . '_' . $self->run->id . '.db_variants.indel.'.$chr.'.xml';
				if ( !-e $dbvar_file_name ) {
					modules::Exception->throw("File $dbvar_file_name doesn't exist");	
				}
				%chr_xml_data = (); #Remove old chr entries
				$chr_xml_data{'db_variants'} = modules::ConfigXML->new("$dbvar_file_name");
				for my $info_filter (@filter_names) {
					
					my $file_name = $self->{confdir} . '/' . $self->{sample_name} . '_' . $self->run->id . '.'.$info_filter.'.indel.'.$chr.'.xml';
					if ( !-e $file_name ) {
						modules::Exception->warning("Skip file xml $file_name doesn't exist");
					} else {
						$chr_xml_data{$info_filter} = modules::ConfigXML->new($file_name);
					}	
				}
				$current_chr = $chr;
			}
		}
		
		my $tumour_other_alleles = "N/A";
		my $ref_allele_count = 0;
		my $var_allele_count = 0;
		my $var_bases;
		my $normal_alleles;
		my $other_alleles = "N/A";
		
		my $ref_base;
	
		#Set ref_base for pileup lookup
	    if ($self->{genome}) {
			$ref_base = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","ref_bases");
		} else {
			$ref_base = $indel->ref_base;
		}	
		
		if ($self->{paired_variants}) {
			my $tumour_base_string;
			my $normal_base_string;

			if ($self->{genome}) {
				$tumour_base_string = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","tumour_base_string");
				$normal_base_string = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","normal_base_string");
			} else {
				$tumour_base_string = $indel->tumour_base_string;
				$normal_base_string = $indel->normal_base_string;	
			}

			my $pl_tumour = modules::PileupLine->new();
			$pl_tumour->base_string_indel($tumour_base_string);
	
			my $var_type;
			if ($self->{genome}){
				$var_type = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","indel_type");
			}else{
				$var_type = $indel->var_type
			}	
	
                        my $indel_length = length($var_base) - 1;
                        $var_bases = $var_base;
                        $var_bases =~ s/([+-])/$1$indel_length/; #change +AA to +2AA for example			
	
			my @tumour_bases = @{$pl_tumour->base_frequencies_indel->bases};
			my @tumour_counts = @{$pl_tumour->base_frequencies_indel->counts};
			my @tumour_remaining_alleles = ();
			
			for ( my $count = 0 ; $count < @tumour_bases ; $count++ ) {
				if ($tumour_bases[$count] eq 'REF') {
					$ref_allele_count = $tumour_counts[$count];
				} elsif ($tumour_bases[$count] eq $var_bases) {
					$var_allele_count = $tumour_counts[$count];
				} else {
					push @tumour_remaining_alleles,"$tumour_bases[$count]:$tumour_counts[$count]";
				}
			}
		
			if (@tumour_remaining_alleles) {
				$tumour_other_alleles = join('_', @tumour_remaining_alleles);
			}
		
			my $pl_normal = modules::PileupLine->new();
			$pl_normal->base_string_indel($normal_base_string);
		
			my @normal_bases = @{$pl_normal->base_frequencies_indel->bases};
			my @normal_counts = @{$pl_normal->base_frequencies_indel->counts};
			my @normal_alleles = ();
		
			for ( my $count = 0 ; $count < @normal_bases ; $count++ ) {
				push @normal_alleles, "$normal_bases[$count]:$normal_counts[$count]"
			}
		
			$normal_alleles = join("_", @normal_alleles);
		} else {
			my $pl_normal = modules::PileupLine->new();
			
			my $normal_base_string;
			if ($self->{genome}) {
				$normal_base_string = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","normal_base_string");
			} else {
				$normal_base_string = $indel->normal_base_string;
			}
			
			$pl_normal->base_string_indel($normal_base_string);	
			
			
			my $indel_length = length($var_base) - 1;
			$var_bases = $var_base;
			$var_bases =~ s/([+-])/$1$indel_length/; #change +AA to +2AA for example
	
			my @normal_bases = @{$pl_normal->base_frequencies_indel->bases};
			my @normal_counts = @{$pl_normal->base_frequencies_indel->counts};
			my @normal_remaining_alleles = ();
			
			for ( my $count = 0 ; $count < @normal_bases ; $count++ ) {
				if ($normal_bases[$count] eq 'REF') {
					$ref_allele_count = $normal_counts[$count];
				} elsif ($normal_bases[$count] eq $var_bases) {
					$var_allele_count = $normal_counts[$count];
				} else {
					push @normal_remaining_alleles,"$normal_bases[$count]:$normal_counts[$count]";
				}
			}
		
			if (@normal_remaining_alleles) {
				$other_alleles = join('_', @normal_remaining_alleles);
			}
		}
		
		my %filter_values;
		my $exon_overlap = 'EXON';
		my $aa_position = "N/A";
		my $aa_length = "N/A";			
		my $gmaf_1000_genomes = "N/A";
		my $clinical_significance = "N/A";
		my $protein_domains = "N/A";
		my $pubmed = "N/A";
		my $exon_intron_count = "N/A";
		#Get the allele freq if relevant
		my $dbsnp_match = "N/A";
		my $dbsnp_var_allele_freq = "N/A";
		my $cosmic_coord = 'N/A';
		my $splice_exon_type = "N/A";
		my $known_variation = "N/A";
		my $filter_exac_indel = "N/A";
		my $filter_gnomad_indel = "N/A";
		my $filter_clinvar_indel = "N/A";
		my $filter_regulatory_custom = "N/A";
		my $filter_mirna = "N/A";
				
		my $read_allele_freq;
		if ($self->{genome}) {
			my $ref_freq = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","ref_base_freq");
			my $var_freq = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","var_base_freq");
			$read_allele_freq = "reads:R".$ref_freq.':V'.$var_freq;
		} else {
			$read_allele_freq = "reads:R".$indel->ref_base_freq.':V'.$indel->var_base_freq;
		}

		#Get the variant filter info
		foreach my $filter_name (@filter_names) {
			next unless exists $chr_xml_data{$filter_name};
			if ($self->{genome}) {
				my $lookup_name = 'indel_'.$filter_name.'_string';
				my $attr_string;
				if ($chr_xml_data{$filter_name}->exists("chr$chr","s$start_coord","e$end_coord","v$xml_var_base",$lookup_name)) {
					$attr_string = $chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base",$lookup_name);
				}
				if ($filter_name eq $EXON_FILTER) {
					($aa_position,$aa_length) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $SPLICE_FILTER) { 
					($splice_exon_type) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $VEP_FILTER) {
					($pubmed,$exon_intron_count,$gmaf_1000_genomes,$clinical_significance,$protein_domains,$known_variation) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
					if ($gmaf_1000_genomes =~ /\(([ACTG\-]+)\)/) {
						#Cover the rare allele cases
						if ($1 eq $ref_base && $xml_var_base =~ /^d/) {
							$gmaf_1000_genomes =~ s/\)/:REF\)/;
						} elsif ($1 eq '-' && $xml_var_base =~ /^i/) {
							$gmaf_1000_genomes =~ s/\)/:REF\)/;
						}
					}
				} elsif ($filter_name eq $DBSNP_FILTER_INDEL) {
					($dbsnp_match, $dbsnp_var_allele_freq) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
					if ($attr_string eq 'novel') {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} elsif ($filter_name eq $EXAC_FILTER_INDEL) {
					($filter_exac_indel) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $GNOMAD_FILTER_INDEL) {
					($filter_gnomad_indel) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $CLINVAR_FILTER_INDEL) {
					($filter_clinvar_indel) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1); 				
				} elsif ($filter_name eq $REGULATORY_CUSTOM_FILTER) {
					($filter_regulatory_custom) =modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);  
				} elsif ($filter_name eq $FILTER_MIRNA) {
					($filter_mirna) =modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);  
				} elsif ($filter_name eq $COSMIC_FILTER) {
					my $cancer_type = modules::Pipeline::get_cancer_type(-run_id=>$self->run->id);
					($cosmic_coord) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-attr=>$attr_string, -genome=>1);
				} else {
					
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$start_coord","e$end_coord","v$var_base","pass") && $chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$end_coord","v$var_base","pass") == 1) {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				}
				
			} else {
				if ($filter_name eq $EXON_FILTER) {
					($aa_position,$aa_length) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-filter_obj=>$filt_indel);
				} elsif ($filter_name eq $SPLICE_FILTER) { 
					($splice_exon_type) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-filter_obj=>$filt_indel);
				} elsif ($filter_name eq $VEP_FILTER) {
					($pubmed,$exon_intron_count,$gmaf_1000_genomes,$clinical_significance,$protein_domains,$known_variation) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-filter_obj=>$filt_indel);
				} elsif ($filter_name eq $DBSNP_FILTER_INDEL) {
					($dbsnp_match, $dbsnp_var_allele_freq) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-filter_obj=>$filt_indel);
					if ($filt_indel->variant_filter_hash->{$filter_name} ne 'nodata' && $filt_indel->variant_filter_hash->{$filter_name}->filterpass == 1) {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} elsif ($filter_name eq $COSMIC_FILTER) {
					my $cancer_type = modules::Pipeline::get_cancer_type(-run_id=>$self->run->id);
					($cosmic_coord) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'indel',-filter_obj=>$filt_indel,-cancer_type=>$cancer_type);
				} else {
					if ($filt_indel->variant_filter_hash->{$filter_name} ne 'nodata') {
						if ($filt_indel->variant_filter_hash->{$filter_name}->filterpass == 1) {
							$filter_values{$filter_name} = 'PASS';
						} else {
							$filter_values{$filter_name} = 'FAIL';
						}
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				
				}
			
				
			
			}
			
			
		}
		
		#Overwrite SYN cases that are actually splice site variants
		if ($splice_exon_type ne 'N/A') {
			$exon_overlap = $splice_exon_type;
		}
				
		my $gene_info = $gene_mapper->get_gene_from_coord(-chr => $chr, -coord => $start_coord);
	
		#If the start coord doesn't overlap use the end coord
		if (!exists $gene_info->{'ensembl'}) {
			$gene_info = $gene_mapper->get_gene_from_coord(-chr => $chr, -coord => $end_coord);
		}			
	
		#Store the annotation_header values for reporting
		my %annotation_headers = ();
	
		$annotation_headers{ensembl} = join(",",keys %{$gene_info->{ensembl}});
		
		if ($self->organism eq 'mouse') {
			$annotation_headers{'gene(mgi)'} = join(",",keys %{$gene_info->{'gene(mgi)'}});
		} else {
			$annotation_headers{'gene(hgnc)'} = join(",",keys %{$gene_info->{'gene(hgnc)'}});
		}
		for my $annotation_header (@{$self->{annotation_headers}}) {
			$annotation_headers{$annotation_header} = join(",",keys %{$gene_info->{$annotation_header}});
		}
		
	
		my $var_type;
		my $var_allele;
		my $var_length;
		my $ref_allele;
		
		if ($self->{genome}) {
			if ($xml_var_base =~ /^i/) {
				$var_type = 'INS';
				$var_allele = $var_base;
				$ref_allele = 'N/A';
				$var_length = length($var_allele) - 1;
			} elsif ($xml_var_base =~ /^d/) {
				$var_type = 'DEL';
				($ref_allele = $var_base) =~ s/[-+]//;
				$var_allele =  $var_base;
				$var_length = $end_coord - $start_coord + 1;
			} else {
				my $indel_str = $chr.':'. $start_coord.'-' .$end_coord;
				modules::Exception->throw("ERROR: Can't get var_type for indel $indel_str vartype $var_type");
			}
		} else {
			if ($indel->var_type eq 'INS') {
				$var_type = 'INS';
				$var_allele = '+' . $indel->affected_bases;
				$ref_allele = 'N/A';
				$var_length = length($var_allele) - 1;
			} elsif ($indel->var_type eq 'DEL') {
				$var_type = 'DEL';
				$ref_allele = $indel->affected_bases;
				$var_allele =  '-'. $indel->affected_bases;
				$var_length = $end_coord - $start_coord + 1;
			} else {
				my $indel_str = $chr.':'. $start_coord.'-' .$end_coord;
				modules::Exception->throw("ERROR: Can't get var_type for indel $indel_str");
			}
			$indel_lookup{$chr}{$start_coord}{$end_coord}{$var_type}{$var_allele} = $indel_id;
			
		}
		
		
		#Flag the special cases with no gene info to avoid perl warnings and to add 'no gene overlap'; cosmic overlaps may
		my $gene_overlap = 1;
		if (!keys %{$gene_info->{'ensembl'}}) {
			$gene_overlap = 0;
		}
		
		my @entries = ();
		
		#Iterate over the headers; the variables will be the same name except in the special cases; this allows us to have custom reports based on the headers
		for my $header (@{$self->{headers}}) {
			my $variable = '$'.$header;
			my $value;
			
			if (exists $filter_values{$header}) {
				$value = $filter_values{$header};
			} elsif ($header eq 'chr') {
				$value = $chr;
			} elsif ($header eq 'start_coord') {
				$value = $start_coord;
			} elsif ($header eq 'end_coord') {
				$value = $end_coord;
			} elsif ($header eq 'var_score') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","indel_score");
				} else {
					$value = $indel->var_score;											
				}
			} elsif ($header eq 'clr_score') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","clr_score");
				} else {				
					$value = $indel->clr_score;
				}
			} elsif ($header eq 'median_quality_score') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","median_quality_score");
				} else {
					$value = $indel->median_quality_score;
				}
			} elsif ($header eq 'read_depth') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","read_depth");
				} else {
					$value = $indel->read_depth;
				}
			}  elsif ($header eq 'variant_class') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$start_coord","e$end_coord","v$xml_var_base","variant_class");
				} else {
					$value = $indel->variant_class;
				}			
			} elsif ($header eq 'final_status') {
		    	#Don't know this yet so use a default of pass
		    	$value = 'OVERALL_PASS'
		    } elsif (exists $annotation_headers{$header} && $annotation_headers{$header} =~ /\w/) {
		    	$value = $annotation_headers{$header};
		    } else {
				#Here the variable name matches
				$value = eval($variable);
			}
			if (!$gene_overlap && exists $annotation_headers{$header}) {
				$value = "No gene overlap";
			} elsif ($value !~ /\w/) {
				if($self->{genome}){
			#		my $chr = $chr_xml_data{'db_variants'}->read("chr");
			#		my $coord = $chr_xml_data{'db_variants'}->read("chr$chr","s");
					print "ERROR $header has no value for $indel_id\n";
				} else {
					my $chr = $indel->chr;
					my $coord = $indel->start_coord;
					print "ERROR $header has no value for $chr $coord\n";
				}
			}
			
			push @entries, $value;
			
			#print "$header $value\n";
		}
		
		my $line_str = join("\t",@entries);
		push @lines,$line_str;
	}
	if ($self->{genome}) {
	    $self->{indel_lookup} = \%indel_lookup;
	}
    
	return @lines;
}

#Get the formatted sv lines
sub _generate_sv_lines {
	my $self = shift;
	my %args = @_;
	my @required_args = (
						'-sv_ids'
			 			);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $sv_ids = $args{-sv_ids};
	my $gene_mapper = $self->{gene_mapper};
	
	my @filter_names;
	if ($self->{genome}) {
		@filter_names = @{$self->{filter_info}};
	} 

	my $debug_chr;
	if (defined $args{-chr}) {
		$debug_chr = $args{'-chr'};
	}


	my @lines = ();
	
	my %chr_xml_data;
	my $current_chr = -1;
	
	foreach my $sv_id (sort {$a<=>$b} keys %{$sv_ids}) {
		
		
		my $sv;
		my $filt_sv;
		my $sv_caller = my $start_coord = my $end_coord = my $sv_type = my $chr;
		
		if ($self->{genome}) {
			($chr,$start_coord,$end_coord,$sv_caller,$sv_type) = split(':',$sv_ids->{$sv_id});
		} 
		
		
		
		next unless $sv_caller eq $self->{sv_caller};
		
		if ($debug_chr) {
			next unless $chr eq $debug_chr;
		}
		
		
		if ($self->{genome}) {
			#Here need to update the db_var xml for the new chromosome
			if ($current_chr ne $chr) {
				
				my $sv_file_name = $self->{confdir} . '/' . $self->{sample_name} . '_' . $self->run->id .'.'.$self->{sv_caller} .'.sv.'.$chr.'.xml';
				if ( !-e $sv_file_name ) {
					modules::Exception->throw("File $sv_file_name doesn't exist");	
				}
				%chr_xml_data = (); #Remove old chr entries
				$chr_xml_data{'sv_info'} = modules::ConfigXML->new("$sv_file_name");
				for my $info_filter (@filter_names) {
					next if $info_filter eq $FILTER_GENE_BP || $info_filter eq $FILTER_SV_EXON_BP;
					my $file_name = $self->{confdir} . '/' . $self->{sample_name} . '_' . $self->run->id . '.'.$info_filter.'.sv.'.$chr.'.xml';
					if ( !-e $file_name ) {
						modules::Exception->warning("Skip file xml $file_name doesn't exist");
					} else {
						$chr_xml_data{$info_filter} = modules::ConfigXML->new($file_name);
					}	
					if ($info_filter eq $FILTER_GENE || $info_filter eq $FILTER_SV_EXON) {
						(my $bp_file = $file_name) =~ s/\.sv/\.bp\.sv/;
						if ( !-e $bp_file ) {
							modules::Exception->warning("Skip file xml $bp_file doesn't exist");
						} else {
							$chr_xml_data{$info_filter.'_breakpoint'} = modules::ConfigXML->new($bp_file);
						}
					}
				}
				$current_chr = $chr;
			}
		}
		
		#chr1,coord1,chr2,coord2,sv_caller,sv_type,sv_id,quality,bases,length,split_reads,paired_reads

		my $chr1 = $chr;
		my $chr2;
		my $coord2;
		
		if ($sv_type eq 'tra') {
			if ($chr_xml_data{'sv_info'}->exists("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'chr2_pair')) {
				$chr2 = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'chr2_pair');
				$coord2 = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'coord2');
			} else {
				next; #With translocations we index via one of two chromosome so skip the one not used
			}
		} else {
			$chr2 = $chr1;
			$coord2 = $end_coord;
		}
		
		if (!$chr_xml_data{'sv_info'}->exists("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'id')) {
			print Dumper $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord");
		}
		
		my $sv_id = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'id');
		my $length = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'len');
		my $split_reads = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'sr');
		my $paired_reads = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'pe');
		my $var_bases = "N/A";
		if ($chr_xml_data{'sv_info'}->exists("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'bases')) {
			$var_bases = $chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'bases');	
		}
		#Quality only reported for delly
		my $quality = $sv_caller eq 'lumpy'?'N/A':$chr_xml_data{'sv_info'}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",'qual');
	
		my $dgv_match;
		my $dgv_freq = "N/A";
		my $dgv_info = "N/A";
		my $breakpoint_gene = 'NONE';
		my $breakpoint_exon = 'NONE';
		my %filter_values;
		$filter_values{'filter_dgv'} = "N/A" if $sv_type eq 'tra' || $sv_type eq 'ins';
		my %bp_genes;
		my %bp_exons;
		my %all_genes;
		
		#Get the variant filter info
		foreach my $filter_name (@filter_names) {
			next unless exists $chr_xml_data{$filter_name};
			if ($self->{genome}) {
				my $lookup_name = 'sv_'.$filter_name.'_string';
				my $attr_string;
				if ($chr_xml_data{$filter_name}->exists("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",$lookup_name)) {
					$attr_string = $chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",$lookup_name);
				}
				if ($filter_name eq $FILTER_DGV) {
					next if $sv_type eq 'tra' || $sv_type eq 'ins';
					($dgv_freq) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'sv',-attr=>$attr_string, -genome=>1);
					$dgv_info = $attr_string;
					if ($attr_string eq 'novel') {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} elsif ($filter_name eq $FILTER_GENE) {
					if ($sv_type ne 'tra' && $sv_type ne 'inv') {
						%all_genes  = map {$_ => 1} split(",",$chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type",$lookup_name)); 
					}
					
				} elsif ($filter_name eq $FILTER_GENE_BP) {
					my @end_genes = my @start_genes = ();
					
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$start_coord","e$start_coord","sv$sv_caller","sv$sv_type",'sv_filter_gene_string')) {
						@start_genes = split(",",$chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$start_coord","sv$sv_caller","sv$sv_type",'sv_filter_gene_string'));
					}
				
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$end_coord","e$end_coord","sv$sv_caller","sv$sv_type",'sv_filter_gene_string')) {
						@end_genes = split(",",$chr_xml_data{$filter_name}->read("chr$chr","s$end_coord","e$end_coord","sv$sv_caller","sv$sv_type",'sv_filter_gene_string'));
					}
					
					for my $gene ( @start_genes ) {
					    $gene =~ s/-DUP//;
					    $bp_genes{$gene}++;
					}
					for my $gene ( @end_genes ) {
					    $gene =~ s/-DUP//;
					    $bp_genes{$gene}++;
					}
					
				} elsif ($filter_name eq $FILTER_SV_EXON_BP) {
					my @end_genes = my @start_genes = ();
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$start_coord","e$start_coord","sv$sv_caller","sv$sv_type",'sv_filter_sv_exon_string')) {
						@start_genes = split(",",$chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$start_coord","sv$sv_caller","sv$sv_type",'sv_filter_sv_exon_string'));
					}
				
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$end_coord","e$end_coord","sv$sv_caller","sv$sv_type",'sv_filter_sv_exon_string')) {
						@end_genes = split(",",$chr_xml_data{$filter_name}->read("chr$chr","s$end_coord","e$end_coord","sv$sv_caller","sv$sv_type",'sv_filter_sv_exon_string'));
					}
					
					for my $gene ( @start_genes ) {
					    $gene =~ s/-DUP//;
					    $bp_exons{$gene}++;
					}
					for my $gene ( @end_genes ) {
					    $gene =~ s/-DUP//;
					    $bp_exons{$gene}++;
					}
					
				} else {
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type","pass") && $chr_xml_data{$filter_name}->read("chr$chr","s$start_coord","e$end_coord","sv$sv_caller","sv$sv_type","pass") == 1) {					
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				}	
			}
		}
		
		my @genes_to_report = ();
		my $gene_count = keys %all_genes;
		
		my %breakpoint_genes;
		if (keys %bp_genes >= 1) {
			@genes_to_report = keys %bp_genes;
			for my $gene (@genes_to_report) {
				my $gene_info = $gene_mapper->get_gene_from_name(-gene_name => $gene);
				for my $gene2 (keys %{$gene_info->{'gene(hgnc)'}}) {
					$breakpoint_genes{$gene2}++;
				}
			}
			$breakpoint_gene = join(",",keys %breakpoint_genes);
		} elsif ($gene_count == 0) {
			#Shouldn't happen
			push @genes_to_report, "NO_GENES";
		} elsif ($gene_count >= 10) {
			push @genes_to_report, "TOO_MANY_GENES ($gene_count)";
		} else {
			@genes_to_report = keys %all_genes;
		}
		my %final_genes;
		
		for my $gene (@genes_to_report) {
			$gene =~ s/-DUP//;
			$final_genes{$gene}++;
		}
		$gene_count =  keys %final_genes;
		#print Dumper \%final_genes;
		
		my $sv_exon_bp;
		if (!keys %bp_exons) {
			$sv_exon_bp = "NO_EXON_BP";
		} else {
			$sv_exon_bp = join(",",keys %bp_exons);
		}
		
		
		
		
		my $first_line = 1;
		
		
		
		
		for my $gene (keys %final_genes) {
			#print "Gene $gene $breakpoint_gene\n";
			my %annotation_headers = ();
			my $gene_info;
			my $gene_overlap = 1;
			
			if ($gene =~ /_GENES/) {
				#Store the annotation_header values for reporting
				for my $annotation_header (@{$self->{annotation_headers}}) {
					$annotation_headers{$annotation_header} = $gene;
				}
				$annotation_headers{ensembl} = $gene;
				$annotation_headers{'gene(hgnc)'} = $gene;
				if ($gene =~ /NO_GENES/) {
					$gene_overlap = 0;
				}
			} else {
				
				$gene_info = $gene_mapper->get_gene_from_name(-gene_name => $gene);
				#If the start coord doesn't overlap use the end coord
				if (!exists $gene_info->{'ensembl'}) {
					modules::Exception->throw("Can't find information for gene");
				}			
			
				#Store the annotation_header values for reporting
			
				$annotation_headers{ensembl} = join(",",keys %{$gene_info->{ensembl}});
				
				if ($self->organism eq 'mouse') {
					$annotation_headers{'gene(mgi)'} = join(",",keys %{$gene_info->{'gene(mgi)'}});
				} else {
					$annotation_headers{'gene(hgnc)'} = join(",",keys %{$gene_info->{'gene(hgnc)'}}) . ' ('.$gene_count .')';
				}
				for my $annotation_header (@{$self->{annotation_headers}}) {
					$annotation_headers{$annotation_header} = join(",",keys %{$gene_info->{$annotation_header}});
				}
			
			}
			
			
			my @entries = ();
			
			my $protein_domains = "N/A";
			my $pubmed = "N/A";
			
			#Iterate over the headers; the variables will be the same name except in the special cases; this allows us to have custom reports based on the headers
			my $ens_header_found = 0;
			my $sv_id_found = 0;
			for my $header (@{$self->{headers}}) {
				if ($header eq 'ensembl') {
					$ens_header_found = 1;
				} 
				if ($header eq 'sv_id') {
					$sv_id_found = 1;
				} 
				
				my $variable = '$'.$header;
				my $value;
				
				if (exists $filter_values{$header}) {
					$value = $filter_values{$header};
				} elsif ($header eq 'chr_1') {
					$value = $chr1;
				} elsif ($header eq 'chr_2') {
					$value = $chr2;
				} elsif ($header eq 'coord1') {
					$value = $start_coord;
				} elsif ($header eq 'coord2') {
					$value = $coord2;
				}  elsif ($header eq 'final_status') {
			    	#Don't know this yet so use a default of pass
			    	$value = 'OVERALL_PASS'
			    } elsif ($header eq 'breakpoint_exon') {
			    	$value = $sv_exon_bp;
			    } elsif (exists $annotation_headers{$header} && $annotation_headers{$header} =~ /\w/) {
			    	$value = $annotation_headers{$header};
			    } else {
					#Here the variable name matches
					$value = eval($variable);
				}
				if (!$gene_overlap && exists $annotation_headers{$header}) {
					$value = "No gene overlap";
				} elsif ($value !~ /\w/) {
					print "ERROR $header has no value for $sv_id\n";
				}
				#Here we don't repeat the same info over and over; just change gene info
				if (!$first_line && !$ens_header_found && $sv_id_found) {
					#$value = '""';
				} 
				push @entries, $value;
				
				#print "$header $value\n";
			}
			
			my $line_str = join("\t",@entries);
			push @lines,$line_str;
			$first_line = 0;
			#print "$line_str\n\n";
		}
		
	}
	return @lines;
}



#Get the formatted summary lines
sub _generate_snv_lines {
	my $self = shift;
	my %args = @_;
	my @required_args = (
						'-snv_ids'
			 			);


    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $debug_chr;
	if (defined $args{-chr}) {
		$debug_chr = $args{'-chr'};
	}

	my %polyphen;
    if (defined $self->{polyphen}) {
    	%polyphen = %{$self->{polyphen}};
    }

    my %polyphen_info;
    if (defined $self->{polyphen_info}) {
    	%polyphen_info = %{$self->{polyphen_info}};
    }
	

	my $snv_ids = $args{-snv_ids};
	
	my $gene_mapper = $self->{gene_mapper};

	my @filter_names;
	if ($self->{genome}) {
		@filter_names = @{$self->{filter_info}};
	} else {
	    foreach my $filter (@{$self->{filter_info}}){
			push @filter_names, $filter->name;
	    }
	}
	#Lookup for creating snv_rows laters (need snv_id from db)
	my %snv_lookup = ();
	
	my @lines = ();
	
	
	my %chr_xml_data;
	my $current_chr = -1;
	
	
	foreach my $snv_id (sort {$a<=>$b} keys %{$snv_ids}) {
		my $snv;
		my $filt_snv;
		my $var_base = my $coord = my $chr;
		
		if ($self->{genome}) {
			($chr,$coord,undef,$var_base) = split(':',$snv_ids->{$snv_id});
		} else {
			($snv) = modules::Adaptors::SNV->search(id => $snv_id);
			$filt_snv = modules::FilteredSNV->new('snv' => $snv,
									    		  'filter_list_order' => $self->{filter_info});
			$var_base = $snv->var_base;
			$coord = $snv->coord;
			$chr = $snv->chr;								
		}

		if ($debug_chr) {
			next unless $chr eq $debug_chr;
		}

		if ($self->{genome}) {
			#Here need to update the db_var xml for the new chromosome
			if ($current_chr ne $chr) {
				my $dbvar_file_name = $self->{confdir} . '/' . $self->{sample_name} . '_' . $self->run->id . '.db_variants.snv.'.$chr.'.xml';
				if ( !-e $dbvar_file_name ) {
					modules::Exception->throw("File $dbvar_file_name doesn't exist");	
				}
				%chr_xml_data = (); #Remove old chr entries
				$chr_xml_data{'db_variants'} = modules::ConfigXML->new("$dbvar_file_name");
				for my $info_filter (@filter_names) {
					my $file_name = $self->{confdir} . '/' . $self->{sample_name} . '_' . $self->run->id . '.'.$info_filter.'.snv.'.$chr.'.xml';
					if ( !-e $file_name ) {
						modules::Exception->warning("Skip file xml $file_name doesn't exist");
					} else {
						$chr_xml_data{$info_filter} = modules::ConfigXML->new($file_name);
					}	
				}
				$current_chr = $chr;
			}
		}
		my $ref_allele_count = 0;
		my $var_allele_count = 0;
		my $tumour_other_alleles = "N/A";
		my $other_alleles = "N/A";
		my $normal_alleles = "N/A";
		my $ref_base;
	
		#Set ref_base for pileup lookup
	    if ($self->{genome}) {
			$ref_base = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","ref_bases");
		} else {
			$ref_base = $snv->ref_base;
		}	
	
	
		#Different handling for paired variant calling....
		if ($self->{paired_variants}) {

			my $tumour_base_string;
			my $normal_base_string;
			if ($self->{genome}) {
				$tumour_base_string = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","tumour_base_string");
				$normal_base_string = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","normal_base_string");

			} else {
				$tumour_base_string = $snv->tumour_base_string;
				$normal_base_string = $snv->normal_base_string;
			}
			
			my $pl_tumour = modules::PileupLine->new();
			$pl_tumour->ref_base($ref_base);
			$pl_tumour->base_string($tumour_base_string);
		
			#my $tumour_other_alleles = join("\:", @{$pl_tumour->base_frequencies->bases}) . '_' . join("\:", @{$pl_tumour->base_frequencies->counts});
			
			my @tumour_bases = @{$pl_tumour->base_frequencies->bases};
			my @tumour_counts = @{$pl_tumour->base_frequencies->counts};
			my @tumour_remaining_alleles = ();
			
			for ( my $count = 0 ; $count < @tumour_bases ; $count++ ) {
				if ($tumour_bases[$count] eq $ref_base) {
					$ref_allele_count = $tumour_counts[$count];
				} elsif ($tumour_bases[$count] eq $var_base) {
					$var_allele_count = $tumour_counts[$count];
				} else {
					push @tumour_remaining_alleles,"$tumour_bases[$count]:$tumour_counts[$count]";
				}
			}
		
			if (@tumour_remaining_alleles) {
				$tumour_other_alleles = join('_', @tumour_remaining_alleles);
			}
		
			my $pl_normal = modules::PileupLine->new();
			$pl_normal->ref_base($ref_base);
			$pl_normal->base_string($normal_base_string);
		
			my @normal_bases = @{$pl_normal->base_frequencies->bases};
			my @normal_counts = @{$pl_normal->base_frequencies->counts};
			my @normal_alleles = ();
		
			for ( my $count = 0 ; $count < @normal_bases ; $count++ ) {
				push @normal_alleles, "$normal_bases[$count]:$normal_counts[$count]"
			}
		
			$normal_alleles = join("_", @normal_alleles);
	
		} else {
			my $normal_base_string;
			if ($self->{genome}) {
				$normal_base_string = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","normal_base_string");
			} else {
				$normal_base_string = $snv->normal_base_string;
			}
			my $pl_normal = modules::PileupLine->new();
			$pl_normal->ref_base($ref_base);
			$pl_normal->base_string($normal_base_string);			
			
			my @normal_bases = @{$pl_normal->base_frequencies->bases};
			my @normal_counts = @{$pl_normal->base_frequencies->counts};
			my @normal_remaining_alleles = ();
			
			for ( my $count = 0 ; $count < @normal_bases ; $count++ ) {
				if ($normal_bases[$count] eq $ref_base) {
					$ref_allele_count = $normal_counts[$count];
				} elsif ($normal_bases[$count] eq $var_base) {
					$var_allele_count = $normal_counts[$count];
				} else {
					push @normal_remaining_alleles,"$normal_bases[$count]:$normal_counts[$count]";
				}
			}
		
			if (@normal_remaining_alleles) {
				$other_alleles = join('_', @normal_remaining_alleles);
			}			
		}
	
		# try and parse the aa change, is possible
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
     
		#Get the allele freq if relevant
		my $dbsnp_var_allele_freq = "N/A";
		my $dbsnp_match = "N/A";
		my $cosmic_coord = 'N/A';
		my $known_variation = "N/A";
		my $read_allele_freq;
		my $filter_exac_snv = "N/A";
		my $filter_gnomad_snv = "N/A";
		my $filter_clinvar_snv = "N/A";
		my $filter_regulatory_custom = "N/A";
		my $filter_mirna = "N/A";
	
		my $gene_info = $gene_mapper->get_gene_from_coord(-chr => $chr, -coord => $coord);;
		if ($self->{genome}) {
			my $ref_freq = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","ref_base_freq");
			my $var_freq = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","var_base_freq");
			$read_allele_freq = "reads:R".$ref_freq.':V'.$var_freq;
		} else {
			$read_allele_freq = "reads:R".$snv->ref_base_freq.':V'.$snv->var_base_freq;
		}
	
		#Get mouse polyphen info
		my $nt_change = $ref_base .'->'. $var_base;
		my $snv_exon_type = 'SYN';
		my $splice_exon_type = 'N/A';
		my %filter_values;
		my $aa_position = "N/A";
		my $aa_length = "N/A";
		
		#Get the relevant snv filter info from each step listed in <filter_info> in report.xml 
		foreach my $filter_name (@filter_names) {
			next unless exists $chr_xml_data{$filter_name};
			if ($self->{genome}) {
				my $lookup_name = 'snv_'.$filter_name.'_string';
				my $attr_string;
				if ($chr_xml_data{$filter_name}->exists("chr$chr","s$coord","e$coord","v$var_base",$lookup_name)) {
					$attr_string = $chr_xml_data{$filter_name}->read("chr$chr","s$coord","e$coord","v$var_base",$lookup_name);
				}

				if ($filter_name eq $EXON_NS_FILTER) { 
					($snv_exon_type,$aa_change,$polyphen_score,$polyphen_prediction,$sift_score,$sift_prediction,$cadd_phred) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string,-genome=>1);		
				} elsif ($filter_name eq $EXON_FILTER) {
					($aa_position,$aa_length) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string,-genome=>1);
				} elsif ($filter_name eq $SPLICE_FILTER) { 
					($splice_exon_type) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string,-genome=>1);
				} elsif ($filter_name eq $VEP_FILTER) {
					($pubmed,$exon_intron_count,$gmaf_1000_genomes,$clinical_significance,$protein_domains,$known_variation) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string,-genome=>1);
					if ($gmaf_1000_genomes =~ /\(([ACTG])\)/) {
						if ($1 eq $ref_base) {
							$gmaf_1000_genomes =~ s/\)/:REF\)/;
						}
					}
				} elsif ($filter_name eq $DBSNP_FILTER) {
					($dbsnp_match,$dbsnp_var_allele_freq) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string,-genome=>1);
					if ($attr_string eq 'novel') {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} elsif ($filter_name eq $EXAC_FILTER) {
					($filter_exac_snv) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $GNOMAD_FILTER) {
					($filter_gnomad_snv) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $CLINVAR_FILTER) {
					($filter_clinvar_snv) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $REGULATORY_CUSTOM_FILTER) {
					($filter_regulatory_custom) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $FILTER_MIRNA) {
					($filter_mirna) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string, -genome=>1);
				} elsif ($filter_name eq $COSMIC_FILTER) {
					($cosmic_coord) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-attr=>$attr_string,-genome=>1);
				} else {
					if ($chr_xml_data{$filter_name}->exists("chr$chr","s$coord","e$coord","v$var_base","pass") && $chr_xml_data{$filter_name}->read("chr$chr","s$coord","e$coord","v$var_base","pass") == 1) {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				}
				
			} else {
				if ($filter_name eq $EXON_NS_FILTER) { 
					($snv_exon_type,$aa_change,$polyphen_score,$polyphen_prediction,$sift_score,$sift_prediction,$cadd_phred) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-filter_obj=>$filt_snv);		
				} elsif ($filter_name eq $EXON_FILTER) {
					($aa_position,$aa_length) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-filter_obj=>$filt_snv);
				} elsif ($filter_name eq $SPLICE_FILTER) { 
					($splice_exon_type) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-filter_obj=>$filt_snv);
				} elsif ($filter_name eq $VEP_FILTER) {
					($pubmed,$exon_intron_count,$gmaf_1000_genomes,$clinical_significance,$protein_domains,$known_variation) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-filter_obj=>$filt_snv);
					if ($gmaf_1000_genomes =~ /\(([ACTG])\)/) {
						if ($1 eq $ref_base) {
							$gmaf_1000_genomes =~ s/\)/:REF\)/;
						}
					}
				} elsif ($filter_name eq $DBSNP_FILTER) {
					($dbsnp_match,$dbsnp_var_allele_freq) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-filter_obj=>$filt_snv);
					if ($filt_snv->snv_filter_hash->{$filter_name} ne 'nodata' && $filt_snv->snv_filter_hash->{$filter_name}->filterpass == 1) {
						$filter_values{$filter_name} = 'PASS';
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				} elsif ($filter_name eq $COSMIC_FILTER) {
					($cosmic_coord) = modules::ReportFilter::get_filter_info(-filter_name=>$filter_name,-var_type=>'snv',-filter_obj=>$filt_snv);
				} else {
					if ($filt_snv->snv_filter_hash->{$filter_name} ne 'nodata') {
						if ($filt_snv->snv_filter_hash->{$filter_name}->filterpass == 1) {
							$filter_values{$filter_name} = 'PASS';
						} else {
							$filter_values{$filter_name} = 'FAIL';
						}
					} else {
						$filter_values{$filter_name} = 'FAIL';
					}
				}
				
			}
			
			
			
			
			
		}
		
		my $polyphen_info = 'N/A';
		if ($self->organism eq 'mouse') {
			if (exists $polyphen{$chr} && exists $polyphen{$chr}{$coord} && exists $polyphen{$chr}{$coord}{$nt_change}) {
				($polyphen_prediction,$polyphen_score) = split(",",$polyphen{$chr}{$coord}{$nt_change});
			}
	
			if (exists $polyphen_info{$chr} && exists $polyphen_info{$chr}{$coord}) {
				$polyphen_info = $polyphen_info{$chr}{$coord};
			}
		}
		
		#Overwrite SYN cases that are actually splice site variants
		if ($splice_exon_type ne 'N/A' && $snv_exon_type eq 'SYN') {
			$snv_exon_type = $splice_exon_type;
		}
		
		#Overwrite SYN cases that are actually intronic
		if ($exon_intron_count =~ /INTRON/ && $snv_exon_type eq 'SYN') {
			$snv_exon_type = 'N/A';
		}

		#Store the annotation_header values for reporting
		my %annotation_headers = ();
	
		#Manually add the annotation fields that aren't listed as annotations in report.xml
		$annotation_headers{ensembl} = join(",",keys %{$gene_info->{ensembl}});
	
		if ($self->organism eq 'mouse') {
			$annotation_headers{'gene(mgi)'} = join(",",keys %{$gene_info->{'gene(mgi)'}});
		} else {
			$annotation_headers{'gene(hgnc)'} = join(",",keys %{$gene_info->{'gene(hgnc)'}});
		}
	
		for my $annotation_header (@{$self->{annotation_headers}}) {
			$annotation_headers{$annotation_header} = join(",",keys %{$gene_info->{$annotation_header}});
		}
		
		#Flag the special cases with no gene info to avoid perl warnings and to add 'no gene overlap'; cosmic overlaps may
		my $gene_overlap = 1;
		if (!keys %{$gene_info->{'ensembl'}}) {
			$gene_overlap = 0;
		}
		
		#Overwrite SYN cases that are actually non-coding site variants
		if (!$gene_overlap && $snv_exon_type eq 'SYN') {
			$snv_exon_type = "N/A";
		}
	
		#Need this lookup for creating snv_rows later
		if (!$self->{genome}) {
			$snv_lookup{$chr}{$coord}{$var_base} = $snv_id;
		}
		my @entries = ();
		
		#Iterate over the headers; the variables will be the same name except in the special cases; this allows us to have custom reports based on the headers
		for my $header (@{$self->{headers}}) {
			my $variable = '$'.$header;
			my $value;
			
				
			if (exists $filter_values{$header}) {
				$value = $filter_values{$header};
			} elsif ($header eq 'chr') {
				$value = $chr;
			} elsif ($header eq 'coord') {
				$value = $coord;
			} elsif ($header eq 'ref_allele') {
				$value = $ref_base;
			} elsif ($header eq 'var_allele') {
				$value = $var_base;
			} elsif ($header eq 'snv_score') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","snv_score");
				} else {
					$value = $snv->snv_score;											
				}
			}  elsif ($header eq 'median_quality_score') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","median_quality_score");
				} else {
					$value = $snv->median_quality_score;
				}
			} elsif ($header eq 'read_depth') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","read_depth");
				} else {
					$value = $snv->read_depth;						
				}
			} elsif ($header eq 'clr_score') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","clr_score");
				} else {				
					$value = $snv->clr_score;
				}
			} elsif ($header eq 'snv_class') {
				if ($self->{genome}) {
					$value = $chr_xml_data{'db_variants'}->read("chr$chr","s$coord","e$coord","v$var_base","snv_class");
				} else {				
					$value = $snv->snv_class;
				}
			} elsif ($header eq 'final_status') {
		    	#Don't know this yet so use a placeholder
		    	$value = 'OVERALL_PASS'
		    } elsif (exists $annotation_headers{$header} && $annotation_headers{$header} =~ /\w/) {
		    	#Most annotation headers are set 
		    	$value = $annotation_headers{$header};
		    } else {
				#Here the variable name matches
				$value = eval($variable);
			}
			

			#Cosmic won't have gene information
			if (!$gene_overlap && exists $annotation_headers{$header}) {
				$value = "No gene overlap";
			} elsif ($value !~ /\w/) {
				if ($self->{genome}){
					my ($chr, $coord,undef,undef) = split(":",$snv_ids->{$snv_id});
					print "ERROR $header has no value for $chr $coord\n";
				} else{
					my $chr = $snv->chr;
					my $coord = $snv->coord;
					print "ERROR $header has no value for $chr $coord\n";
				}
			}
			
			push @entries, $value;
		}
	
		my $line_str = join("\t",@entries);
		push @lines,$line_str;
		
		
	}
	
	if (!$self->{genome}) {
	    $self->{snv_lookup} = \%snv_lookup;
	}
    
	return @lines;
}

#Getter for passed rows
sub get_fail_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};		
	
	for my $field_array ( @{$self->{fail_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} elsif ($variant_type eq 'indel') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr_1}]}{$field_array->[$header_col_lookup{chr_2}]}{$field_array->[$header_col_lookup{coord1}]}{$field_array->[$header_col_lookup{coord2}]}{$field_array->[$header_col_lookup{sv_caller}]} = $line;
		}
	}
	
	return \%sorted_rows;
}

#Getter for passed rows
sub get_pass_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{pass_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} elsif ($variant_type eq 'indel') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr_1}]}{$field_array->[$header_col_lookup{chr_2}]}{$field_array->[$header_col_lookup{coord1}]}{$field_array->[$header_col_lookup{coord2}]}{$field_array->[$header_col_lookup{sv_caller}]} = $line;
		}
	}
	
	return \%sorted_rows;
}

#Getter for allele rows
sub get_allele_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{allele_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} elsif ($variant_type eq 'indel') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		}  else {
			$sorted_rows{$field_array->[$header_col_lookup{chr_1}]}{$field_array->[$header_col_lookup{chr_2}]}{$field_array->[$header_col_lookup{coord1}]}{$field_array->[$header_col_lookup{coord2}]}{$field_array->[$header_col_lookup{sv_caller}]} = $line;
		}
	}
	
	return \%sorted_rows;
}

#Getter for no_freq rows
sub get_no_freq_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{no_freq_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		}
	}
	
	return \%sorted_rows;
}

#Getter for all the rows
sub get_all_rows {
	my $self = shift;
	my %sorted_rows;	
	my %header_col_lookup = %{$self->{header_lookup}};
	my $variant_type = $self->{variant_type};
	
	for my $field_array ( @{$self->{all_rows}}) {
		my $line = join("\t",@{$field_array});
		if ($variant_type eq 'snv') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{coord}]} = $line;
		} elsif ($variant_type eq 'indel') {
			$sorted_rows{$field_array->[$header_col_lookup{chr}]}{$field_array->[$header_col_lookup{start_coord}]}{$field_array->[$header_col_lookup{end_coord}]}{$field_array->[$header_col_lookup{var_type}]}{$field_array->[$header_col_lookup{var_allele}]} = $line;
		} else {
			$sorted_rows{$field_array->[$header_col_lookup{chr1}]}{$field_array->[$header_col_lookup{coord1}]}{$field_array->[$header_col_lookup{coord2}]}{$field_array->[$header_col_lookup{sv_caller}]}{$field_array->[$header_col_lookup{sv_type}]} = $line;
		}
	}

	return \%sorted_rows;
}


#Get the snvs to report on from the primary filter field
sub _get_variants {
	my $self = shift;
	my ($chr,$start,$end) = @_;
    my $var_id_iterator;  
    my  %var_ids = ();
	my $variant_type = $self->{variant_type};
    
    if ($self->{genome}) {
    	my @xml_files;    	
    	my $sample_name = $self->{sample_name};

    	for my $primary_filter (@{$self->{primary_filters}}){
    		my $file_name = $self->{confdir} . '/' . $sample_name . '_' . $self->run->id . '.'.$primary_filter.'.'.$variant_type.'.xml';
    		if ( !-e $file_name ) {
    			modules::Exception->throw("File $file_name doesn't exist");	
    		}
    		push @xml_files,$file_name;
    		if ($primary_filter eq $FILTER_GENE) {
    			#Here we have extra bp xml file
    			(my $file_name_bp = $file_name) =~ s/\.sv/\.bp\.sv/;
    			if ( !-e $file_name_bp ) {
    				modules::Exception->throw("File $file_name_bp doesn't exist");	
    			}
    			push @xml_files,$file_name_bp;
    		}
    	}
    
    	my $var_count = 0;
    	#Now get all the coordinates for lookup
    	for my $xml_file (@xml_files) {
    		my $xml_obj = modules::ConfigXML->new($xml_file);
    		my $xml_struct = $xml_obj->{'xml_ref'};
    		for my $chr (sort keys %{$xml_struct}) {
    			for my $start (sort {my ($a_coord) = $a =~ /(\d+)/;  my ($b_coord) = $b =~ /(\d+)/;$a_coord<=>$b_coord} keys %{$xml_struct->{$chr}}) {
    				for my $end (sort {my ($a_coord) = $a =~ /(\d+)/;  my ($b_coord) = $b =~ /(\d+)/;$a_coord<=>$b_coord} keys %{$xml_struct->{$chr}{$start}}) {
    					if ($variant_type eq 'sv') {
    						for my $sv_caller (keys %{$xml_struct->{$chr}{$start}{$end}}) {
    							(my $sv_caller_tmp = $sv_caller) =~ s/sv//;
    							next unless $sv_caller_tmp eq $self->{sv_caller}; #Make sure sv_caller matches
    							for my $sv_type (keys %{$xml_struct->{$chr}{$start}{$end}{$sv_caller}}) { 
    								(my $start_tmp = $start) =~ s/s//;
	    							(my $end_tmp  = $end) =~ s/e//;
	    							(my $chr_tmp = $chr) =~ s/chr//;
	    							(my $sv_type_tmp = $sv_type) =~ s/sv//;
	    							my $var_base_tmp;
	    							if ($xml_file =~ /\.bp\.sv/) {
										#Here we only load translocations
										$var_ids{$var_count} = "$chr_tmp:$start_tmp:$end_tmp:$sv_caller_tmp:$sv_type_tmp" if $sv_type_tmp eq 'tra';
	    							} else {
		    							$var_ids{$var_count} = "$chr_tmp:$start_tmp:$end_tmp:$sv_caller_tmp:$sv_type_tmp"; #here we store info for lookup
	    							}
							    	$var_count++;
    							}
    						}
    					} else {	
	    					for my $var_base (keys %{$xml_struct->{$chr}{$start}{$end}}) {
	    						(my $start_tmp = $start) =~ s/s//;
	    						(my $end_tmp  = $end) =~ s/e//;
	    						(my $chr_tmp = $chr) =~ s/chr//;
	    						my $var_base_tmp;
	    						
	    						if ($variant_type eq 'snv') {
	    							($var_base_tmp = $var_base) =~ s/v//;
	    						} elsif ($var_base =~ /^vd/) {
	    							($var_base_tmp = $var_base) =~ s/vd/-/;
	    						} elsif ($var_base =~ /^vi/) {
	    							($var_base_tmp = $var_base) =~ s/vi/+/;
	    						}
							    $var_ids{$var_count} = "$chr_tmp:$start_tmp:$end_tmp:$var_base_tmp"; #here we store info for lookup
							    $var_count++;
	    					}
    					}
    				}
    			}
    		}	
    	}
    } else {
	    my @primary_filter_ids;
	    
	    foreach my $primary_filter (@{$self->{primary_filters}}){
			push @primary_filter_ids, $primary_filter->id;
	    }
	
		my $var_count = 0;
	
	
	    foreach my $primary_filter_id (@primary_filter_ids) {
	
			if (defined $chr && defined $start && defined $end) {
				if ($variant_type eq 'snv') {
			    	$var_id_iterator = modules::Adaptors::SNV->search_region_by_experiment_and_passed_filter($self->run->id, 
													     $primary_filter_id,
													     $chr, 
													     $start, 
													     $end);
				} else {
					$var_id_iterator = modules::Adaptors::Variant->search_region_by_experiment_and_passed_filter($self->run->id, 
													     $primary_filter_id,
													     $chr, 
													     $start, 
													     $end);
				}
		
			} elsif (defined $chr) {
				if ($variant_type eq 'snv') {
			    	$var_id_iterator = modules::Adaptors::SNV->search_chr_by_experiment_and_passed_filter($self->run->id, 
													     $primary_filter_id,
													     $chr);
				} else {
					$var_id_iterator = modules::Adaptors::Variant->search_chr_by_experiment_and_passed_filter($self->run->id, 
													     $primary_filter_id,
													     $chr);
				}
				
			} else {
				if ($variant_type eq 'snv') {
			    	$var_id_iterator = modules::Adaptors::SNV->search_by_experiment_and_passed_filter($self->run->id, 
													     $primary_filter_id);
				} else {
					$var_id_iterator = modules::Adaptors::Variant->search_by_experiment_and_passed_filter($self->run->id, 
													     $primary_filter_id);
				}
			}
		
			while (my $var_id = $var_id_iterator->next){
				$var_count++;
			    $var_ids{$var_id}++;
			}
	    }
    }
    return %var_ids;
}

#Print the summary stats to a file; used in the pipeline for passed snvs
sub summarize {
	my $self = shift;
	my %args = @_;
	
	my @required_args = ('-summary_file');

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $variant_type = $self->{variant_type};
	
	my $summary_file = $args{-summary_file};
	my %filter_set = %{$self->{filter_set}};
	my %report_filters = %{$self->{report_set}};
	my @primary_filters =  @{$self->{primary_filters}};
	
	my %filter_count = %{$self->{filter_count}};
	
	my %filter_files = %{$self->{filter_files}};

	
	
	open(my $SUMM, ">$summary_file")
	    or modules::Exception->throw("Unable to open file [$summary_file] for writing");
	
	  # record details of filters applied
	
	print $SUMM "Filters applied:\n";
	
	foreach my $filter (keys %filter_set){
		if ($filter_set{$filter} =~ /allele_freq=([0-9\.]+)/) {
			#Special dbsnp for human case
			print $SUMM "\t$filter (PASS DEF: no_match OR allele_freq <= $1)\n";
		} elsif ($filter_set{$filter} !~ /^[0-9]+$/) {
			#Non numerical filter
			print $SUMM "\t$filter (PASS DEF: $filter_set{$filter})\n";
		} elsif ($filter_set{$filter} == 1) {
			print $SUMM "\t$filter (PASS DEF: no_match)\n";	
		} else {
	    	print $SUMM "\t$filter (PASS DEF: >= $filter_set{$filter})\n";
		}
	}
	
	for my $file_filter ( keys %filter_files ) {
	    print $SUMM "\t$file_filter (PASS DEF: no_match)\n";
	}
	
	my @required_filters_names;
	
	if ($self->{genome}) {
		@required_filters_names = @primary_filters;
	} else {
		for my $primary_filter (@primary_filters) {
			push @required_filters_names, $primary_filter->name;
		}
	}
	
	my $required_filter = join (" OR ",@required_filters_names);
	print $SUMM "\tMatch EITHER $required_filter\n";
	
	
	my $variant_count;
	if ($variant_type eq 'snv') {
		$variant_count = $filter_count{filter_snv}{pass};
	} elsif ($variant_type eq 'indel') {
		$variant_count = $filter_count{filter_indel}{pass};
	} else {
		$variant_count = $filter_count{filter_sv}{pass};
	}
	print $SUMM "\nTotal $required_filter $variant_type called: $variant_count\n";
	
	
	#Breakdown the snv_type_exon fields
	my @all_rows = @{$self->{all_rows}};
		
	my %header_col_lookup = %{$self->{header_lookup}};
	
	
	my %var_classes = ();
	my %var_type_exons = ();
	
	for my $row ( @all_rows ) {
		if ($variant_type eq 'snv') {
			$var_type_exons{$row->[$header_col_lookup{snv_exon_type}]}++;
		} elsif ($variant_type eq 'indel') {
			$var_classes{$row->[$header_col_lookup{var_type}]}++;
			$var_type_exons{$row->[$header_col_lookup{exon_overlap}]}++;
			
		} else {
			$var_classes{$row->[$header_col_lookup{sv_type}]}++;
		}
	}
	
	for my $var_class_key (keys %var_classes) {
		my $pass_percent = sprintf("%.2f",$var_classes{$var_class_key} / $variant_count * 100 );
		my $lc_entry = lc($var_class_key);
		print $SUMM "\t$lc_entry MATCH: $var_classes{$var_class_key} ($pass_percent%)\n";
	}
	
	my %splice_entries = ();
	my $splice_total = 0;
	my $exonic_all = 0;	
	for my $var_exon_key (keys %var_type_exons) {
		if ($var_exon_key !~ /SPLICE/) {
			$exonic_all += $var_type_exons{$var_exon_key};
		} else {
			$splice_entries{$var_exon_key} = $var_type_exons{$var_exon_key};
			$splice_total += $var_type_exons{$var_exon_key};
			next;
		}				
		my $pass_percent = sprintf("%.2f",$var_type_exons{$var_exon_key} / $variant_count * 100 );
		my $lc_entry = lc ($var_exon_key);
		print $SUMM "\t$lc_entry MATCH: $var_type_exons{$var_exon_key} ($pass_percent%)\n";	
		
	}
	if ($splice_total) {
		my $pass_percent = sprintf("%.2f",$splice_total / $variant_count * 100 );
	
		print $SUMM "\ttotal splice MATCH: $splice_total ($pass_percent%) ; ";
	
		my $splice_line = '';
		for my $splice_entry (sort keys %splice_entries) {
			$splice_line .= "$splice_entry->$splice_entries{$splice_entry} , ";
		}
		$splice_line =~ s/ , $//;
		print $SUMM "$splice_line\n"; 
	} else {
		print $SUMM "\tfilter_splice MATCH: 0 (0%)\n";
	}
	
	if ($exonic_all > 0) {
		my $pass_percent = sprintf("%.2f",$exonic_all / $variant_count * 100 );
		print $SUMM "\tfilter_exon MATCH: $exonic_all ($pass_percent%)\n";
	} else {
		print $SUMM "\tfilter_exon MATCH: 0 (0%)\n";
	}
		
	#Print the report filters	
	for my $report_filter ( sort keys %report_filters ) {
		if (defined $filter_count{$report_filter}{pass} ) {
	   		my $pass_percent = sprintf("%.2f",$filter_count{$report_filter}{pass} / $variant_count * 100 );
	  		print $SUMM "\t$report_filter PASS: $filter_count{$report_filter}{pass} ($pass_percent%)\n";
		} else {
			print $SUMM "\t$report_filter PASS: 0 (0%)\n";
		}
	}
	 
	my $gene_filtered_count = 0;
	for my $gene ( keys %{$self->{filtered_genes}} ) {
	    $gene_filtered_count += $self->{filtered_genes}{$gene};
	}
	if ($gene_filtered_count > 0) {
		my $pass_percent = sprintf("%.2f",$gene_filtered_count / $variant_count * 100 );
		print $SUMM "\tproblem genes FAIL: $gene_filtered_count ($pass_percent%) (", join(",", keys %{$self->{filtered_genes}}), ")\n";
	} else {
		print $SUMM "\tproblem genes FAIL: 0 (0%)\n";		
	}
	
	my %sorted_pass_rows = %{$self->{sorted_pass_rows}};
	
	# tally a few items of information between rows
	my $ss_passed = 0;
	my $exonic_passed = 0;
	my $heter_count_passed = 0;
	my $homo_count_passed = 0;
	my $other_count_passed = 0;
	my $insertions_passed = 0;
	my $deletions_passed = 0;
	my $trans_passed = 0;
	my $inv_passed = 0;
	my $dup_passed = 0;
	
	my @depth_values;
	my @quality_values;
	my %genes;

	my @pass_rows = @{$self->{pass_rows}};
	
	#Get some meta info from the passed rows
	foreach my $pass_row (@pass_rows) {
		
		if ($variant_type eq 'indel') {
			if ($pass_row->[$header_col_lookup{var_type}] eq 'DEL') {
				$deletions_passed++;
			} elsif ($pass_row->[$header_col_lookup{var_type}] eq 'INS') {
				$insertions_passed++;
			} else {
				modules::Exception->throw("ERROR: INDEL is none of accepted types");
			}
		}
		
		if ($variant_type eq 'sv') {
			if ($pass_row->[$header_col_lookup{sv_type}] eq 'del') {
				$deletions_passed++;
			} elsif ($pass_row->[$header_col_lookup{sv_type}] eq 'ins') {
				$insertions_passed++;
			} elsif ($pass_row->[$header_col_lookup{sv_type}] eq 'inv') {
				$inv_passed++;
			} elsif ($pass_row->[$header_col_lookup{sv_type}] eq 'dup') {
				$dup_passed++;
			} elsif ($pass_row->[$header_col_lookup{sv_type}] eq 'tra') {
				$trans_passed++;
			} else {
				modules::Exception->throw("ERROR: INDEL is none of accepted types");
			}
		}
		
		
		if ($variant_type eq 'snv') {
			if ($pass_row->[$header_col_lookup{snv_exon_type}] =~ /SPLICE/) {
				$ss_passed++;
			} else {
				$exonic_passed++;
			}
		} elsif ($variant_type eq 'indel') {
			if ($pass_row->[$header_col_lookup{exon_overlap}] =~ /SPLICE/) {
				$ss_passed++;
			} else {
				$exonic_passed++;				
			}
		}
		
		
		my $mutant_depth;
		if ($variant_type eq 'snv') {
			$mutant_depth = $pass_row->[$header_col_lookup{read_depth}];
		
		    if ($mutant_depth =~ /[0-9]/) {
				push @depth_values, $mutant_depth;
		    } else {
				modules::Exception->warning("Non-numeric characters in read depth column");
		    }
		    if ($pass_row->[$header_col_lookup{median_quality_score}] =~ /[0-9\.]/) {
				push @quality_values, $pass_row->[$header_col_lookup{median_quality_score}];
		    } else {
				modules::Exception->warning("Non-numeric characters in median base quality column");
		    }
		}
		
		$genes{$pass_row->[$header_col_lookup{$self->{gene_col_name}}]}++;
	}
	
	# print summary output
	
	print $SUMM "\n";
	print $SUMM "Total passed variants: " . scalar @pass_rows . "\n";
	print $SUMM "\tTotal homozygous: $homo_count_passed\n" if $homo_count_passed>0;
	print $SUMM "\tTotal heterozygous: $heter_count_passed\n" if $heter_count_passed>0;
	print $SUMM "\tTotal other: $other_count_passed\n" if $other_count_passed>0;
	print $SUMM "\tTotal deletions: $deletions_passed\n" if $deletions_passed>0;
	print $SUMM "\tTotal insertions: $insertions_passed\n" if $insertions_passed>0;
	
	print $SUMM "\nBreakdown of passed variants:\n\tDistinct genes: " .  scalar (keys %genes) . "\n";
	print $SUMM "\tExonic: $exonic_passed\n";
	print $SUMM "\tSplice-site: $ss_passed\n";
	print $SUMM "\tMedian Depth: " . modules::Utils->median(\@depth_values) . "\n" if @depth_values;
	print $SUMM "\tMedian Quality: " . modules::Utils->median(\@quality_values) . "\n" if @quality_values;
	

	
	print $SUMM "\tGenes with more than one variant: ";
	my $gene_str = 'NONE';
	foreach my $gene (keys %genes){
	    if ($genes{$gene} > 1){
	    	$gene_str .= ",$gene";
	    }
	}
	$gene_str =~ s/NONE,//;
	print $SUMM "($gene_str)\n";
	
	close($SUMM);
			
}

sub insert_SNV_rows {
	my $self = shift;
	my %args = @_;
	my @snv_rows = ();
	my @snv_annotations = ();
	
	my %header_col_lookup = %{$self->{header_lookup}};
	my @headers = @{$self->{headers}};
	my %snv_lookup = %{$self->{snv_lookup}};
	
	#Too many rows to insert (db timeout) so just insert passed and rare_alleles 
	my @rows = ();
	if (@{$self->{all_rows}} > 1000) {
		if ($self->{pass_rows}) {
			push @rows, @{$self->{pass_rows}};
		}
		if ($self->{allele_rows}) {
			push @rows, @{$self->{allele_rows}};
		}
		if ($self->{no_freq_rows}) {
			push @rows, @{$self->{no_freq_rows}};
		}
	} else {
		@rows = @{$self->{all_rows}}
	}
	
	
		
	for my $line ( @rows) {
		#my @cols = split("\t",$line);
		my $chr = $line->[$header_col_lookup{chr}];
		my $coord = $line->[$header_col_lookup{coord}];
		my $var_base = $line->[$header_col_lookup{var_allele}];
		
		if (!exists $snv_lookup{$chr}{$coord}{$var_base}) {
			modules::Exception->throw("ERROR: Can't get snv id for $chr:$coord $var_base");
		}
		
		my $snv_id = $snv_lookup{$chr}{$coord}{$var_base};
		
		my %snv_row = (
						snv_id => $snv_id,
						run_id => $self->run->id
						);
						
#		for my $db_column (@db_columns) {
#			my $column_name = $db_column->{column_name};
#			next if $column_name eq 'snv_id' || $column_name eq 'run_id' || $column_name eq 'id';
#			$snv_row{$column_name} = $line->[$header_col_lookup{$column_name}];
#		}

		push @snv_rows, \%snv_row;
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'snv_rows', -data=>\@snv_rows) if @snv_rows;
	
	my $last_col_count = keys %header_col_lookup;
	$last_col_count--;
	
	my %snv_row_lookup = ();
	
	my $snv_row_iterator = modules::Adaptors::SNV_Row->search(run_id => $self->run->id);	
	
	#Get all the ids for the inserted records
	while (my $snvrow_db = $snv_row_iterator->next){
		my ($snv_obj) = modules::Adaptors::SNV->search(id => $snvrow_db->snv_id);
		$snv_row_lookup{$snv_obj->chr}{$snv_obj->coord}{$snv_obj->var_base} = $snvrow_db->id; 
	}

	
	
	#Now we need to add the snv_annotations
	for my $line ( @rows) {
		my $chr = $line->[$header_col_lookup{chr}];
		my $coord = $line->[$header_col_lookup{coord}];
		my $var_base = $line->[$header_col_lookup{var_allele}];
		
		if (!exists $snv_row_lookup{$chr}{$coord}{$var_base}) {
			modules::Exception->throw("ERROR: Can't find snv_row entry for $chr:$coord:$var_base");
		}
		
		#Get the snv_row_id using the snv_lookup 
		my $snv_row_id = $snv_row_lookup{$chr}{$coord}{$var_base};

		#This is range of columns not in snv_row; these are added as annotations
		for my $col_count (0..$last_col_count) {
			my $column_value = $line->[$col_count];
			$column_value =~ s/'/"/g if defined $column_value;
			if ($snv_row_id !~ /\d/) {
				print "ERROR: for $line\n";
			}
			my %snv_annotation = (
									column_name => $headers[$col_count],
									column_value => $column_value,
									column_number => $col_count,
									snv_row_id => $snv_row_id
								  );
			push @snv_annotations, \%snv_annotation;
		}
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'snvrow_values', -data=>\@snv_annotations) if @snv_annotations;
	
	
}

sub insert_Indel_rows {
	my $self = shift;
	my %args = @_;
	my @indel_rows = ();
	my @indel_annotations = ();
	
	my %header_col_lookup = %{$self->{header_lookup}};
	my @headers = @{$self->{headers}};
	my %indel_lookup = %{$self->{indel_lookup}};
	
	#Too many rows to insert (db timeout) so just insert passed and rare_alleles 
	my @rows = ();
	if (@{$self->{all_rows}} > 1000) {
		if ($self->{pass_rows}) {
			push @rows, @{$self->{pass_rows}};
		}
		if ($self->{allele_rows}) {
			push @rows, @{$self->{allele_rows}};
		}
		if ($self->{no_freq_rows}) {
			push @rows, @{$self->{no_freq_rows}};
		}		
	} else {
		@rows = @{$self->{all_rows}}
	}
	
	
	for my $line ( @rows) {
		#my @cols = split("\t",$line);
		my @cols = @{$line};
		my $chr = $cols[$header_col_lookup{chr}];
		my $start_coord = $cols[$header_col_lookup{start_coord}];
		my $end_coord = $cols[$header_col_lookup{end_coord}];
		my $var_base = $cols[$header_col_lookup{var_allele}];
		my $var_type = $cols[$header_col_lookup{var_type}];
		
		if (!exists $indel_lookup{$chr}{$start_coord}{$end_coord}{$var_type}{$var_base}) {
			modules::Exception->throw("ERROR: Can't get indel id for $chr:$start_coord $var_base");
		}
		
		my $indel_id = $indel_lookup{$chr}{$start_coord}{$end_coord}{$var_type}{$var_base};
		
		my %indel_row = (
						variant_id => $indel_id,
						run_id => $self->run->id
						);
		
#		for my $db_column (@db_columns) {
#			my $column_name = $db_column->{column_name};
#			next if $column_name eq 'variant_id' || $column_name eq 'run_id' || $column_name eq 'id';
#			$indel_row{$column_name} = $cols[$header_col_lookup{$column_name}];
#		}
		
		push @indel_rows, \%indel_row;
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'variant_rows', -data=>\@indel_rows) if @indel_rows;
	
	#Get the first column for annotations; snv_exon_type is the last in the snv_row table
	my $first_col_count = 0;
	my $last_col_count = keys %header_col_lookup;
	$last_col_count--;
	
	my %indel_row_lookup = ();
	
	my $indel_row_iterator = modules::Adaptors::Variant_Row->search(run_id => $self->run->id);	
	
	#Get all the ids for the inserted records
	while (my $indelrow_db = $indel_row_iterator->next){
		my ($var_db) = modules::Adaptors::Variant->search(id => $indelrow_db->variant_id);
		my $bases;
		if ($var_db->var_type eq 'DEL') {
			$bases = '-'. $var_db->affected_bases;
		} elsif ($var_db->var_type eq 'INS') {
			$bases = '+'. $var_db->affected_bases;
		} else {
			modules::Exception->throw("ERROR: Unknown variant type");
		}
		$indel_row_lookup{$var_db->chr}{$var_db->start_coord}{$var_db->end_coord}{$var_db->var_type}{$bases} = $indelrow_db->id; 
	}
	
	#Now we need to add the indel_annotations
	for my $line ( @rows) {
		my $chr = $line->[$header_col_lookup{chr}];
		my $start_coord = $line->[$header_col_lookup{start_coord}];
		my $end_coord = $line->[$header_col_lookup{end_coord}];
		my $var_base = $line->[$header_col_lookup{var_allele}];
		my $indel_type = $line->[$header_col_lookup{var_type}];
		
		if (!exists $indel_row_lookup{$chr}{$start_coord}{$end_coord}{$indel_type}{$var_base}) {
			modules::Exception->throw("ERROR: Can't get indel row id for $chr:$start_coord $var_base");
		}
		
		#Get the snv_row_id using the snv_lookup 
		my $indel_row_id = $indel_row_lookup{$chr}{$start_coord}{$end_coord}{$indel_type}{$var_base};

		#This is range of columns not in snv_row; these are added as annotations
		for my $col_count ($first_col_count..$last_col_count) {
			my $column_value = $line->[$col_count];
			$column_value =~ s/'/"/g if defined $column_value;
			my %indel_annotation = (
									column_name => $headers[$col_count],
									column_value => $column_value,
									column_number => $col_count,
									variant_row_id => $indel_row_id
								  );
			push @indel_annotations, \%indel_annotation;
		}
	}
	modules::Adaptors::BulkInsert->insert(-table_name=>'variantrow_values', -data=>\@indel_annotations) if @indel_annotations;
	
	
	
}


#Create the final xml file from the calculated lines
sub final_xml {
	my $self = shift;
	my %args = @_;
    
    
    my @required_args = (
    					'-xml_file',
			 			);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    my $xml_file = $args{-xml_file};
    
    my %header_col_lookup = %{$self->{header_lookup}};	
    my @headers = @{$self->{headers}};
    
    my $first_col_count = 0;
	my $last_col_count = keys %header_col_lookup;
	$last_col_count--;
	
    my %xml_data = ();
    
	if ($self->{variant_type} eq 'snv') {
		for my $line (@{$self->{all_rows}}) {
			my $chr = $line->[$header_col_lookup{chr}];
			my $coord = $line->[$header_col_lookup{coord}];
			my $var_base = $line->[$header_col_lookup{var_allele}];
			my $xml_key = $chr.":".$coord;
			for my $col_count ($first_col_count..$last_col_count) {
				my $column_value = $line->[$col_count];
				$column_value =~ s/'/"/g if defined $column_value;
				my $column_name = $headers[$col_count];
				$column_name =~ s/\(/_/g; #xml tags can't have brackets
				$column_name =~ s/\)//g;
				$xml_data{$xml_key}{$coord}{$var_base}{$column_name} = $column_value;
			}
		}
	} elsif ($self->{variant_type} eq 'indel') {
		for my $line (@{$self->{all_rows}}) {
			my $chr = $line->[$header_col_lookup{chr}];
			my $start_coord = $line->[$header_col_lookup{start_coord}];
			my $end_coord = $line->[$header_col_lookup{end_coord}];
			my $var_base = $line->[$header_col_lookup{var_allele}];
			my $xml_key = $chr.":".$start_coord;
			for my $col_count ($first_col_count..$last_col_count) {
				my $column_value = $line->[$col_count];
				$column_value =~ s/'/"/g if defined $column_value;
				my $column_name = $headers[$col_count];
				$column_name =~ s/\(/_/g; #xml tags can't have brackets
				$column_name =~ s/\)//g;
				$xml_data{$xml_key}{$end_coord}{$var_base}{$column_name} = $column_value;
			}
		}
	} elsif ($self->{variant_type} eq 'sv') {
		my $sv_caller = $self->{sv_caller};
		for my $line (@{$self->{all_rows}}) {
			my $chr1 = $line->[$header_col_lookup{chr_1}];
			my $chr2 = $line->[$header_col_lookup{chr_2}];
			my $start_coord = $line->[$header_col_lookup{coord1}];
			my $sv_type = $line->[$header_col_lookup{sv_type}];
			my $end_coord; 
			if ($sv_type eq 'tra') {
				$end_coord = $start_coord;
			} else {
				$end_coord = $line->[$header_col_lookup{coord2}];
				
			}
			
			my $xml_key = $chr1.":".$start_coord;
			for my $col_count ($first_col_count..$last_col_count) {
				my $column_value = $line->[$col_count];
				$column_value =~ s/'/"/g if defined $column_value;
				my $column_name = $headers[$col_count];
				$column_name =~ s/\(/_/g; #xml tags can't have brackets
				$column_name =~ s/\)//g;
				$xml_data{$xml_key}{$end_coord}{$sv_caller}{$sv_type}{$column_name} = $column_value;
			}
		}
	} else {
		modules::Exception->throw("ERROR: variant type must be snv or indel");
	}
	
	#Now create the xml file
	my $var_xml = modules::VariantXML->new(-outdir=>$self->{confdir});
	if ($self->{variant_type} eq 'sv') {
		$var_xml->create_var_xml(-file_name=>$xml_file,-data=>\%xml_data, -chr=>'all',-sv=>1);
	} else {
		$var_xml->create_var_xml(-file_name=>$xml_file,-data=>\%xml_data, -chr=>'all');
	}
	my $full_xml_file = $self->{confdir} . '/'.$xml_file;
	$var_xml->split_xml(-file_name=>$full_xml_file);
}


return 1;
