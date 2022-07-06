package modules::VEP;

use strict;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;
#./variant_effect_predictor/variant_effect_predictor.pl -i vep.in --poly b --sift b --canonical --cache --dir ~/Desktop/ --offline -o STDOUT --coding_only

sub new {
    my ($class, @args) = @_;

	my @required_args = (
			             -exon,
			             -input_file,   
						 -executable_path,
						 -db_dir,        
						 -working_dir,     
						 -source_type
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

    my $self = bless \%args, $class;
    #exon flag
    $self->{exon} = $args{-exon};
    $self->input_file($args{-input_file});
    $self->executable_path($args{-executable_path});
    $self->db_dir($args{-db_dir});
    $self->working_dir($args{-working_dir});
    $self->source_type($args{-source_type});
    
    if (defined $args{-output_file}) {
    	$self->output_file($args{-output_file});
    } else {
    	$self->output_file('vep.out');
    }
    
    
    return $self;
}



sub source_type {
	my ($self, $source_type) = @_;

    if (defined $source_type){
		$self->{'source_type'} = $source_type;
    } elsif (! defined $self->{'source_type'}) {
		modules::Exception->throw("variant_predictor source_type not set");
    }

    return $self->{'source_type'};
}

sub executable_path {
    my ($self, $executable_path) = @_;
    if (defined $executable_path 
	&& -e $executable_path){
		$self->{'executable_path'} = $executable_path;
    } elsif (! defined $self->{'executable_path'}) {
		$self->{'executable_path'} = 'variant_effect_predictor.pl';
    }

    return $self->{'executable_path'};
}

sub db_dir {
    my ($self, $db_dir) = @_;

    if (defined $db_dir){
		$self->{'db_dir'} = $db_dir;
    } elsif (! defined $self->{'db_dir'}) {
		modules::Exception->throw("variant_predictor database directory not set");
    }

    return $self->{'db_dir'};
}

sub working_dir {
    my ($self, $working_dir) = @_;

    if (defined $working_dir){
		$self->{'working_dir'} = $working_dir;
    } elsif (! defined $self->{'working_dir'}) {
		modules::Exception->throw("variant_predictor working directory not set");
    }

    return $self->{'working_dir'};
}

sub output_file {
    my ($self,$output_file) = @_;
	if (defined $output_file){
		$self->{'output_file'} = $output_file;
    } elsif (! defined $self->{'output_file'}) {
	    $self->{'output_file'} = $self->{'working_dir'} . '/vep.out';
    }
    return $self->{'output_file'};
}

sub input_file {
    my ($self, $input_file) = @_;

    if (defined $input_file){
		$self->{'input_file'} = $input_file;
    } elsif (! defined $self->{'input_file'}) {
		modules::Exception->throw("No input file set");
    }

    return $self->{'input_file'};
}

sub run {
    my ($self) = @_;
	my $pipe_config = modules::Pipeline::get_pipe_conf();
	my $args;
	
	#vep args set in pipe.xml for each source type and whether using exon args or not
	if ($self->{exon}) {
		$args = $pipe_config->read($self->source_type,'binaries','variant_predictor','args_exon');
	} else {
		$args = $pipe_config->read($self->source_type,'binaries','variant_predictor','args_all');
	}


    my $command 
	= $self->executable_path 
	#Get sift and polyphen scores, canonical state, use the cache file(no db), and only find coding exon entries
	. ' ' . $args
	. ' --dir ' . $self->db_dir
	. ' --o ' . $self->output_file
	. ' --i ' . $self->input_file;

    print STDERR "Running command: $command\n";
	my $sys_call = modules::SystemCall->new();
    
    my $return_value = $sys_call->run($command);
    sleep(60);
    
    if ($return_value) {
		return 1;    	
    } else {
	    return 0;
    }
    
#    if (system($command)){
#		modules::Exception->throw("Command: $command\nError: Non-zero exit status.");
#		return 0;
#    }

}


sub parse_result {
    my ($self) = @_;
    my $output_file = $self->output_file;
    open(my $OUTPUT, $output_file) || modules::Exception->throw("ERROR: Can't open output file $output_file");

	my %grouping_data = ();

	while (<$OUTPUT>){
		chomp;
		#Skip headers
		next if $_ =~ /^#/;
		
		my ($identifier, $coord_str, $var_base, $ens_gene, $ens_transcript, undef, $aa_type, undef, undef, $aa_pos, $aa_change, $codon_change, $rs, $attribute_str ) = split /\t/;
		
		if ($rs eq '-') {
			$rs = "N/A";
		}
		
		my $var_type;
		if ($identifier =~ /\/[ACTG]([ACTG]+)/) {
			$var_type = 'INS';	
			$var_base = '+' . $1;
		}  elsif ($identifier =~ /\-$/) {
			$var_type = 'DEL';
			($var_base) = $identifier =~ /_(\D+)\/\-$/;
			$var_base = '-' . $var_base;
		} elsif ($identifier =~ /_[TCGA]{1,3}\/[TCGA]{1,3}$/) {
			$var_type = 'SNV'; #up to a codon long
		} else {
			modules::Exception->throw("ERROR: Identifier $identifier doesn't match any var type\n");
		}
		
		my $chr = my $start = my $end;
		if ($coord_str =~ /([0-9XYMT]+):(\d+)\-(\d+)/) {
			$chr = $1;
			$start = $2;
			$end = $3;
			if ($identifier =~ /\/[ACTG][ACTG]/) {
				#insertions need to change insert start coord
				$start++;
			}
		} else {
			($chr,$start) = $coord_str =~ /([0-9XYMT]+):(\d+)/;
			$end = $start;
		}
		
		
		if ($self->{exon}) {
			#Handling for exon specific fields	

			#Only record results from canonical transcripts and skip synonomous mutants
			next unless $attribute_str =~ /CANONICAL/; 

			my $aa_ref = my $aa_var;
			
			if ($aa_change =~ /([A-Z\*])\/([A-Z\*])/){
				$aa_ref = $1;
				$aa_var = $2;
			} else {
				#synonomous variants
			    next;
			}
			
			my $polyphen_pred = "N/A";
			my $polyphen_score = "N/A";
			my $sift_pred = "N/A";	
			my $sift_score = "N/A";
			my $cadd_phred = "N/A";
		
			my @attribute_pairs = split /;/, $attribute_str;
			#Sample attribute line: 
			#PolyPhen=possibly_damaging(0.593);CANONICAL=YES;SIFT=deleterious(0);EXON=5/11
			for my $attribute_pair ( @attribute_pairs ) {			    
			    if ($attribute_pair =~ /PolyPhen/) {
			    	($polyphen_pred, $polyphen_score) = $attribute_pair =~ /PolyPhen=(.+)\(([0-9\.]+)\)/;
			    } elsif ($attribute_pair =~ /SIFT/) {
			    	($sift_pred, $sift_score) = $attribute_pair =~ /SIFT=(.+)\(([0-9\.]+)\)/;
			    } elsif ($attribute_pair =~ /CADD_PHRED=([^\;]+)/) {
			    	$cadd_phred = $1;
			    }
			}
		
			
		
		
			$aa_ref = 'Stop' if ($aa_ref eq '*');
			$aa_var = 'Stop' if ($aa_var eq '*');
		
			#Just keep one copy of redundant entries
			$grouping_data{$chr}{$start}{$end}{"$aa_ref->$aa_var"} = [ $chr, $start, $end, $var_type, $aa_type, $var_base, $ens_gene, $ens_transcript, $aa_ref, $aa_var, $aa_pos, $polyphen_pred, $polyphen_score, $sift_pred, $sift_score, $cadd_phred];
		
			#push @parsed_result, [$aa_type, $chr, $start, $ens_gene, $ens_transcript, $aa_ref, $aa_var, $polyphen_pred, $polyphen_score, $sift_pred, $sift_score];
	    } else {
	    	#Handling for running on 'all' bases
	    	
	    	my $gmaf = "N/A"; #human only
	    	my $domain = "N/A"; 
	    	my $pubmed = "N/A";
	    	my $clinical = "N/A"; #possibly human only
	    	my $exon_str = "N/A";
	    	my $cadd_phred = "N/A";
	    	my $canonical = $attribute_str =~ /CANONICAL/ ? 1:0;

			my @attribute_pairs = split /;/, $attribute_str;
			#Sample attribute line: 
			#DOMAINS=Pfam_domain:PF01108,Superfamily_domains:SSF49265;CCDS=CCDS33544.1;PUBMED=16690980,16885196;CLIN_SIG=non-pathogenic;GMAF=G:0.0325
			my $exon = my $intron = 0;
			
			for my $attribute_pair ( @attribute_pairs ) {
			    
			    if ($attribute_pair =~ /DOMAINS=(.*)/ && $canonical) {
			    	$domain=$1;
			    } elsif ($attribute_pair =~ /AF=([0-9\.]+)/) {
			    	#my ($base,$freq) = split(':',$1);
			    	#$gmaf = $freq .'('. $base.')';
			    	$gmaf=$1;
			    } elsif ($attribute_pair =~ /CLIN_SIG=(.*)/) {
			    	$clinical = $1;
			    } elsif ($attribute_pair =~ /PUBMED=(.*)/) {
			    	$pubmed = $1;
			    } elsif ($attribute_pair =~ /EXON/ && $canonical){
			    	($exon_str) = $attribute_pair =~ /(EXON=\d+\/\d+)/;
			    	$exon_str =~ s/=/->/;
			    	$exon = 1;
			    } elsif ($attribute_pair =~ /INTRON/ && $canonical){
			    	($exon_str) = $attribute_pair =~ /(INTRON=\d+\/\d+)/;
			    	$exon_str =~ s/=/->/;
			    	$intron = 1;
			    } elsif ($attribute_pair =~ /CADD_PHRED=([^\;]+)/) {
			    	$cadd_phred = $1;
			    } 
			}
			
			if ($intron) {
		    	$grouping_data{$chr}{$start}{$end}{intron} = [$chr, $start, $end, $var_type, $aa_type, $var_base, $rs, $gmaf, $domain, $pubmed, $clinical, $exon_str, $ens_gene, $ens_transcript,$aa_type,$cadd_phred];
			} elsif ($exon) {
				$grouping_data{$chr}{$start}{$end}{exon} = [$chr, $start, $end, $var_type, $aa_type, $var_base, $rs, $gmaf, $domain, $pubmed, $clinical, $exon_str, $ens_gene, $ens_transcript,$aa_type,$cadd_phred];				
			} else {
				$grouping_data{$chr}{$start}{$end}{neither} = [$chr, $start, $end, $var_type, $aa_type, $var_base, $rs, $gmaf, $domain, $pubmed, $clinical, $exon_str, $ens_gene, $ens_transcript,$aa_type,$cadd_phred];
			}
	    	
	    }
	    
	}
    close($OUTPUT);

	my @parsed_result;

	for my $chr ( sort keys %grouping_data ) {
	    for my $start_coord ( sort {$a<=>$b} keys %{$grouping_data{$chr}} ) {
	    	for my $end_coord (keys %{$grouping_data{$chr}{$start_coord}}) {
	    		if ($self->{exon}) {
	    			for my $aa_change (keys %{$grouping_data{$chr}{$start_coord}{$end_coord}}) {
	    				push @parsed_result, $grouping_data{$chr}{$start_coord}{$end_coord}{$aa_change};
	    			}
	    		} else {
	    			#Preferably use exon overlap
	    			if (exists $grouping_data{$chr}{$start_coord}{$end_coord}{exon}) {
	    				push @parsed_result, $grouping_data{$chr}{$start_coord}{$end_coord}{exon};
	    			} elsif (exists $grouping_data{$chr}{$start_coord}{$end_coord}{intron}) {
	    				push @parsed_result, $grouping_data{$chr}{$start_coord}{$end_coord}{intron};
	    			} else {
	    				push @parsed_result, $grouping_data{$chr}{$start_coord}{$end_coord}{neither};
	    			}
	    		}
	    	}
		}
	}
    return \@parsed_result;
}

return 1;
