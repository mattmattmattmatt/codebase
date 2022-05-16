package modules::Lumpy;

use strict;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;

sub new {
    my ($class, @args) = @_;

	my @required_args = (
						 -executable_path,
						 -input_bam,
						 -output_vcf,
						 -tmpdir
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

    my $self = bless \%args, $class;
    $self->output_vcf($args{-output_vcf});
    $self->executable_path($args{-executable_path});
    $self->input_bam($args{-input_bam});
    $self->tmpdir($args{-tmpdir});
    
    if (defined $args{-lumpy_exclude}) {
		$self->lumpy_exclude($args{-lumpy_exclude});    	
    }
    return $self;
}

sub output_vcf {
	my ($self, $output_vcf) = @_;

    if (defined $output_vcf){
		$self->{'output_vcf'} = $output_vcf;
    } elsif (! defined $self->{'output_vcf'}) {
		modules::Exception->throw("output_vcf not set");
    }

    return $self->{'output_vcf'};
}

sub executable_path {
    my ($self, $executable_path) = @_;
    if (defined $executable_path 
	&& -e $executable_path){
		$self->{'executable_path'} = $executable_path;
    } elsif (! defined $self->{'executable_path'}) {
		$self->{'executable_path'} = 'lumpyexpress';
    }

    return $self->{'executable_path'};
}

sub lumpy_exclude {
    my ($self, $lumpy_exclude) = @_;

    if (defined $lumpy_exclude){
		$self->{'lumpy_exclude'} = $lumpy_exclude;
    } elsif (! defined $self->{'lumpy_exclude'}) {
		return 0;
    }

    return $self->{'lumpy_exclude'};
}

sub input_bam {
    my ($self,$input_bam) = @_;
	if (defined $input_bam){
		$self->{'input_bam'} = $input_bam;
    } elsif (! defined $self->{'input_bam'}) {
	    modules::Exception->throw("input_bam not set");
    }
    return $self->{'input_bam'};
}

sub tmpdir {
    my ($self,$tmpdir) = @_;
        if (defined $tmpdir){
                $self->{'tmpdir'} = $tmpdir;
    } elsif (! defined $self->{'tmpdir'}) {
            modules::Exception->throw("tmpdir not set");
    }
    return $self->{'tmpdir'};
}


sub run {
    my ($self) = @_;

    my $command 
	= $self->executable_path 
	. ' -T ' . $self->tmpdir
	. ' -o ' . $self->output_vcf
	. ' -B ' . $self->input_bam . ' --verbose';
	
	if ($self->lumpy_exclude) {
		$command .= ' -x '.$self->lumpy_exclude;
	} 

    print STDERR "Running command: $command\n";
	my $sys_call = modules::SystemCall->new();
    
    my $return_value = $sys_call->run($command);
    sleep(60);
    
    if ($return_value) {
		return 1;    	
    } else {
	    return 0;
    }
 

}


sub parse_result {
    my ($self) = @_;
    my $output_vcf = $self->output_vcf;
    open(my $OUTPUT, $output_vcf) || modules::Exception->throw("ERROR: Can't open output file $output_vcf");
	my %sv_data = ();
	my %lumpy_ids = ();
	my %trans_seen = (); #Lumpy reports trans twice so only record once -> handled automatically for other events because we sort by chromosome
	
	while (<$OUTPUT>) {
		chomp;
		next if /^#/;
		my ($chr1,$start1,$id,$ref_base,$var_base,undef,undef,$gff_str) = split("\t");
		next if $chr1 =~ /GL0/;
		next if /SECONDARY/;  #partner event 
		$chr1 =~ s/chr//;
		(my $short_id = $id) =~ s/_\d+$//;
		if (exists $lumpy_ids{$short_id}) {
			next;
		}
		my $event_type;
		my $chr2 = my $start2;
		my $length = 0;
		if ($gff_str =~ /SVTYPE=BND/) {
			
			($chr2,$start2) = $var_base =~ /c?h?r?([0-9XYM]+):(\d+)/;
			next if $var_base =~ /GL0/;
			#Determine the type
			if ($chr1 ne $chr2) {
				$event_type = 'tra';
				if (exists $trans_seen{"$chr1:$start1"}) {
					next;
				}
				$trans_seen{"$chr2:$start2"}++;
			} elsif ($gff_str =~ /\+\+/ || $gff_str =~ /\-\-/) {
				$event_type = 'inv';
			} else {
				modules::Exception->throw("ERROR: Cannot classify SV event $_");
			}
			if ($chr1 eq $chr2) {
				$length = abs($start2-$start1);
			}
		} else {
			($event_type) = $gff_str =~ /SVTYPE=([A-Z]+);/;
			$event_type = lc($event_type);
			$chr2 = $chr1;
			($start2) = $gff_str =~ /END=(\d+);/; 
			($length) = $gff_str =~ /SVLEN=\-?(\d+)/;
		}
		if ($start1 > $start2 && $event_type ne 'tra') {
			my $start2_tmp = $start2;
			$start2 = $start1;
			$start1 = $start2_tmp;
		}
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{length} = $length;
		my ($pe_count) = $gff_str =~ /PE=(\d+)/;
		my ($sr_count) = $gff_str =~ /SR=(\d+)/;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{pe} = $pe_count;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{sr} = $sr_count;
		
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{gff} = $gff_str;	
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{id} = $short_id;
		
		
	}

	
    return \%sv_data;
}

return 1;
