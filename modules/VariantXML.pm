package modules::VariantXML;
use modules::ConfigXML;
use modules::Exception;
use Data::Dumper;
use vars qw(%SIG);
use strict;

sub new { 
	my ($class, @args) = @_;

    my $self = bless {}, $class;

    my %args = @args;

	my @required_args = (
			             -outdir
						 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

	my $outdir = $args{'-outdir'};

	if (!-d $outdir) {
		modules::Exception->throw("ERROR: outdir $outdir doesn't exist");
	}
	
    $self->outdir($outdir);
   	
    return $self;
}

#Add the outdir
sub outdir {
	my ($self, $outdir) = @_;

    if (defined $outdir) {
		$self->{'outdir'} = $outdir;
    } elsif (! defined $self->{'outdir'}) {
		modules::Exception->throw("outdir not set");
    }

    return $self->{'outdir'};
}

#subroutine requires a file_name, and hash of unique key->value pairs, and optional append flag and chr flag
#hash must have keys of the form {chr:start_coord:end_coord}{var_base1} etc
sub create_var_xml { 
	my ($self, @args) = @_;
	my %args = @args;
	my @required_args = (
			             -file_name,
			             -data,
			             '-chr'
					 );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument $required_arg not set");
		}
    }
    
	my $chr_flag = $args{'-chr'};
	
	#Set first and last to use for opening tags
    my $first = my $last = 0;
    if ($chr_flag eq '1' || $args{-chr} eq 'all') {
    	$first = 1;
    } 
    if ($chr_flag eq 'Y' || $args{-chr} eq 'all') {
    	$last = 1;
    }
    
    my %data = %{$args{'-data'}};
   	
    
    my $out_file = $self->outdir() . '/' . $args{'-file_name'};

    if ($first) {
    	#Open the file first time
    	open(OUT,">$out_file") || modules::Exception->throw("Can't open file to write $out_file\n");
		print OUT "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n<vars>\n\n";
    } else {
    	#Otherwise append
    	open(OUT,">>$out_file") || modules::Exception->throw("Can't open file to write $out_file\n");
    }

	#No overlaps; make sure to close vars tags
	if (!keys %data) {
		if ($last) {
			print OUT "\n</vars>\n";
			close OUT;
		}
		return;
	}

   	my ($first_key) = keys %data;
    if (my $bad_key   = $self->_check_key_format($first_key)) {
    	modules::Exception->throw("ERROR: Bad key for $bad_key; Must be of the form chr:start_coord");
    }

    
	my $last_chr = -1; #flag for when we change chromosome
	
	for my $chr_str (sort {my ($a_chr,$a_coord) = $a =~ /(.+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /(.+):(\d+)/; $a_chr cmp $b_chr ||  $a_coord <=> $b_coord} keys %data) {
		my ($chr,$start);
		if ($chr_str =~ /(.+):(\d+)/) {
			($chr,$start) = ($1,$2);
		} else {
			modules::Exception->throw("ERROR: chr_str $chr_str is the wrong format");
		}

		if ($last_chr eq '-1') {
			#First entry for file
			print OUT "\t<chr$chr>\n";
		} elsif ($chr ne $last_chr) {
			print OUT "\t</chr$last_chr>\n\n\t<chr$chr>\n";
		}
		
		print OUT "\t\t<s$start>\n";
		
		for my $end (sort {$a<=>$b} keys %{$data{$chr_str}}) {
			
			if ($end < $start) {
				modules::Exception->throw("ERROR: End coord ($end) is less than start ($start)\n");
			}
			
			print OUT "\t\t\t<e$end>\n";
		
			if (defined $args{'-sv'}) {
				#Different structure for SVs
				for my $sv_caller (keys %{$data{$chr_str}{$end}}) {
					print OUT "\t\t\t\t<sv$sv_caller>\n";
					for my $sv_type (keys %{$data{$chr_str}{$end}{$sv_caller}}) {
						print OUT "\t\t\t\t\t<sv$sv_type>\n";
						for my $key (sort keys %{$data{$chr_str}{$end}{$sv_caller}{$sv_type}}) {
							my $filter_key = $self->_delimit($key);
							my $filter_value = $self->_delimit($data{$chr_str}{$end}{$sv_caller}{$sv_type}{$key});
							print OUT "\t\t\t\t\t\t<$filter_key>$filter_value</$filter_key>\n" 
						}
						print OUT "\t\t\t\t\t</sv$sv_type>\n";
					}	
					print OUT "\t\t\t\t</sv$sv_caller>\n";
					
				}
				
			} else {
				for my $var_base (keys %{$data{$chr_str}{$end}}) {
					my $var_base_tmp = $var_base;
					$var_base_tmp =~ s/^\+/i/;
					$var_base_tmp =~ s/^\-/d/;
					print OUT "\t\t\t\t<v$var_base_tmp>\n";
					
					for my $key (sort keys %{$data{$chr_str}{$end}{$var_base}}) {
						my $filter_key = $self->_delimit($key);
						my $filter_value = $self->_delimit($data{$chr_str}{$end}{$var_base}{$key});
						print OUT "\t\t\t\t\t<$filter_key>$filter_value</$filter_key>\n" 
					}
					
					print OUT "\t\t\t\t</v$var_base_tmp>\n";
				}
				
			}
		
			print OUT "\t\t\t</e$end>\n";
		}
		
		print OUT "\t\t</s$start>\n";
		$last_chr = $chr;
	}
		
	print OUT "\t</chr$last_chr>\n";
	
	#If it's last chr or we've done all of them we need to close the vars tag
	if ($last) {
		print OUT "\n</vars>\n";
	}     
	close OUT;    
	
	#Store the file if we need to interogate later
	push @{$self->{files}}, $out_file;
	    
}



#Delimit all relevant characters to make sure xml is legit
sub _delimit {
	my ($self,$string) = @_;
	$string =~ s/&/&amp;/g;
	$string =~ s/</&lt;/g;
	$string =~ s/>/&gt;/g;
	return $string;
}

#Check keys are the correct format
sub _check_key_format {
	my ($self,$key) = @_;
	
	if (!defined $key) {
		return 'undefined key';
	}
	
	if ($key !~ /[0-9XYM]+:\d+$/) {
		return $key;
	}
	
	return 0;
	
}

#split by chr into individual reports
sub split_xml {
	my ($self, @args) = @_;
	my %args = @args;
	my @required_args = (
			             -file_name
					 	);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument $required_arg not set");
		}
    }
    
    my $file = $args{-file_name};
    if ( !-e $file ) {
    	modules::Exception->throw("File $file doesn't exist");	
    }
    open(FILE,"$file") || modules::Exception->throw("Can't open file $file\n");
    
    my $chr = -1;
    my $fh;
    while (<FILE>) {
    	if (/<chr([0-9XY]+)>/) {
    		close $fh unless $chr eq '-1';
    		$chr = $1;
    		(my $chr_file = $file) =~ s/.xml/.$chr.xml/;
    		open($fh,">$chr_file") || die "Can't open file $chr_file for writing"; 
    		print $fh "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n\n<vars>\n\n$_";	
    	} elsif (/<\/chr[0-9XY]+>/) {
    		print $fh "$_</vars>\n";
    	} elsif (/\/vars/) {
			last;
    	} elsif ($chr ne '-1') {
    		print $fh $_;
    	}
    	
    }
    
    
}

#load db_var xml into object
sub load_xml_keys {
	my ($self, @args) = @_;
	my %args = @args;
	my @required_args = (
			             -file_name
					   );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument $required_arg not set");
		}
    }
    my $xml_file = $self->outdir() . '/' . $args{'-file_name'};

    
    if (!-e $xml_file) {
    	modules::Exception->throw("ERROR: Couldn't open file for parsing $xml_file");
    }
    
    if (exists $self->{xml_loaded} && exists $self->{xml_loaded}{$xml_file}) {
    	return;
    } else {
    	$self->{xml_loaded}{$xml_file}++;
    }
    
    my $var_type = $args{'-file_name'} =~ /indel/?'indel':'snv'; 
    $var_type = 'sv' if defined $args{'-sv'};
    
    my $xml_obj = modules::ConfigXML->new($xml_file);
    my $xml_struct = $xml_obj->{'xml_ref'};
    
	#Check we haven't already loaded the file
#	if (exists $self->{dbvar} && exists $self->{dbvar}{$var_type}) {
#		return;
#	}
	    
    
 	for my $chr (sort keys %{$xml_struct}) {
		my @array_of_keys = ();
 		(my $chr_tmp = $chr) =~ s/chr//;
    	for my $start (sort {my ($a_coord) = $a =~ /(\d+)/;  my ($b_coord) = $b =~ /(\d+)/;$a_coord<=>$b_coord} keys %{$xml_struct->{$chr}}) {
    		for my $end (sort {my ($a_coord) = $a =~ /(\d+)/;  my ($b_coord) = $b =~ /(\d+)/;$a_coord<=>$b_coord} keys %{$xml_struct->{$chr}{$start}}) {
    			if ($var_type eq 'sv') {
    				for my $sv_caller (keys %{$xml_struct->{$chr}{$start}{$end}}) {
	    				for my $sv_type (keys %{$xml_struct->{$chr}{$start}{$end}{$sv_caller}}) {
	    					my %key_data = ();
		    				(my $start_tmp = $start) =~ s/s//;
		    				(my $end_tmp  = $end) =~ s/e//;
		    				$key_data{chr} = $chr_tmp;
							$key_data{start_coord} = $start_tmp;
							$key_data{end_coord} = $end_tmp;
							$key_data{sv_type} = $sv_type;
							$key_data{sv_caller} = $sv_caller;
							push @array_of_keys, \%key_data;
	    				}
    				}
    			} else {
	    			for my $var_base (keys %{$xml_struct->{$chr}{$start}{$end}}) {
	    				my %key_data = ();
	    				(my $start_tmp = $start) =~ s/s//;
	    				(my $end_tmp  = $end) =~ s/e//;
	    				my $var_base_tmp;
	    						
	    				if ($var_type eq 'snv') {
	    					($var_base_tmp = $var_base) =~ s/v//;
	    				} elsif ($var_base =~ /^vd/) {
	    					($var_base_tmp = $var_base) =~ s/vd/-/;
	    				} elsif ($var_base =~ /^vi/) {
	    					($var_base_tmp = $var_base) =~ s/vi/+/;
	    				}
						#Generate struct like db obj to keep lookup consistent
						$key_data{chr} = $chr_tmp;
						$key_data{start_coord} = $start_tmp;
						$key_data{end_coord} = $end_tmp;
						$key_data{var_bases} = $var_base_tmp;
						push @array_of_keys, \%key_data;
	    			}
    			}
    		}	
    	}
    	$self->{dbvar}{$var_type}{$chr_tmp} = \@array_of_keys;
    }
}

#search xml keys by chr/coord and optionally var_bases to see if entry exists
sub search_xml_by_coord {
	my ($self, @args) = @_;
	my %args = @args;
	my @required_args = (
			             -var_type,
			             -chrom,
			             -start_coord,
			             -end_coord
					   );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument $required_arg not set");
		}
    }
    
    my $var_type = $args{'-var_type'};
    my $chr = $args{'-chrom'};
    my $start_coord = $args{'-start_coord'};
    my $end_coord = $args{'-end_coord'};
    
	if (!exists $self->{dbvar} || !exists $self->{dbvar}{$var_type} || !exists $self->{dbvar}{$var_type}{$chr}) {
		return 0;
	}
	
	if (defined $args{'-var_bases'}) {
		my $var_bases = $args{'-var_bases'};
		
		for my $var_obj (@{$self->{dbvar}{$var_type}{$chr}}) {
			if ($var_obj->{start_coord} == $start_coord && $var_obj->{end_coord} == $end_coord && $var_obj->{var_bases} eq $var_bases) {
				return 1;
			}
		}
		return 0;
	} elsif ($var_type eq 'sv') {
		if (! defined $args{'-sv_caller'}) {
			modules::Exception->throw("Required argument -sv_caller not set");
		}
		for my $var_obj (@{$self->{dbvar}{$var_type}{$chr}}) {
			if ($var_obj->{start_coord} == $start_coord && $var_obj->{end_coord} == $end_coord && $var_obj->{sv_caller} eq $args{'sv_caller'}) {
				return 1;
			}
		}
	} else {
		for my $var_obj (@{$self->{dbvar}{$var_type}{$chr}}) {
			if ($var_obj->{start_coord} == $start_coord && $var_obj->{end_coord} == $end_coord) {
				return 1;
			}
		}
		return 0;
	}
	
	
	
}

#search xml keys from dbvar
sub search_xml_keys {
	my ($self, @args) = @_;
	my %args = @args;
	my @required_args = (
			             -var_type
					   );

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument $required_arg not set");
		}
    }
    
    my $var_type = $args{'-var_type'};
    
    #Check we haven't already loaded the file
	if (!exists $self->{dbvar} || !exists $self->{dbvar}{$var_type}) {
		if($var_type eq 'sv'){
			modules::Exception->warning("WARNING: no lookup data for $var_type");
		}else{
			modules::Exception->throw("ERROR: No lookup data for $var_type");
		}
	}
	
	my @return_objs = ();
	#Either return the a single chromosome or all the chromosomes
	if (defined $args{'-chr'}) {
		my $chr = $args{'-chr'};
		if (!exists $self->{dbvar}{$var_type}{$chr}) {
			modules::Exception->throw("ERROR: No lookup data for $var_type for chr $chr");
		}
		@return_objs = @{$self->{dbvar}{$var_type}{$chr}};
	} else {
		for my $chr_data (sort keys %{$self->{dbvar}{$var_type}}) {
			push @return_objs, @{$self->{dbvar}{$var_type}{$chr_data}};
		}
	}
	return \@return_objs;
}



1;
