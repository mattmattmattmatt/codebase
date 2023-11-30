package modules::Delly;

use strict;
use modules::Exception;
use modules::SystemCall;
use Data::Dumper;

sub new {
    my ($class, @args) = @_;

	my @required_args = (
						 -executable_path,
						 -bcftools,
						 -input_bam,
						 -output_bcf,
						 -call_args,
						 -sv_type,
						 -ref_fasta,
						 -bcftools
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }

    my $self = bless \%args, $class;
    $self->output_bcf($args{-output_bcf});
    $self->executable_path($args{-executable_path});
    $self->input_bam($args{-input_bam});
    $self->call_args($args{-call_args});
    $self->sv_type($args{-sv_type});
    $self->ref_fasta($args{-ref_fasta});
    $self->bcftools($args{-bcftools});
    
    if (defined $args{-delly_exclude}) {
		$self->delly_exclude($args{-delly_exclude});    	
    }
    if (defined $args{-fast}) {
	$self->{fast} = 1;
    }
    return $self;
}

sub output_bcf {
	my ($self, $output_bcf) = @_;

    if (defined $output_bcf){
		$self->{'output_bcf'} = $output_bcf;
    } elsif (! defined $self->{'output_bcf'}) {
		modules::Exception->throw("output_bcf not set");
    }

    return $self->{'output_bcf'};
}

sub bcftools {
	my ($self, $bcftools) = @_;

    if (defined $bcftools){
		$self->{'bcftools'} = $bcftools;
    } elsif (! defined $self->{'bcftools'}) {
		modules::Exception->throw("bcftools not set");
    }

    return $self->{'bcftools'};
}

sub executable_path {
    my ($self, $executable_path) = @_;
    if (defined $executable_path 
	&& -e $executable_path){
		$self->{'executable_path'} = $executable_path;
    } elsif (! defined $self->{'executable_path'}) {
		$self->{'executable_path'} = 'delly';
    }

    return $self->{'executable_path'};
}

sub call_args {
    my ($self, $call_args) = @_;

    if (defined $call_args){
		$self->{'call_args'} = $call_args;
    } elsif (! defined $self->{'call_args'}) {
		modules::Exception->throw("call_args not set");
    }

    return $self->{'call_args'};
}

sub delly_exclude {
    my ($self, $delly_exclude) = @_;

    if (defined $delly_exclude){
		$self->{'delly_exclude'} = $delly_exclude;
    } elsif (! defined $self->{'delly_exclude'}) {
		return 0;
    }

    return $self->{'delly_exclude'};
}

sub sv_type {
    my ($self, $sv_type) = @_;

    if (defined $sv_type){
		$self->{'sv_type'} = $sv_type;
    } elsif (! defined $self->{'sv_type'}) {
		modules::Exception->throw("sv_type not set");
    }

    return $self->{'sv_type'};
}

sub ref_fasta {
    my ($self, $ref_fasta) = @_;

    if (defined $ref_fasta){
		$self->{'ref_fasta'} = $ref_fasta;
    } elsif (! defined $self->{'ref_fasta'}) {
		modules::Exception->throw("ref_fasta not set");
    }

    return $self->{'ref_fasta'};
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

sub _convert_to_vcf {
	my ($self, $bcf_file) = @_;
    (my $outvcf = $bcf_file) =~ s/.bcf/.vcf/;
	my $cmd = $self->bcftools . ' view '. $bcf_file . ' > '.$outvcf;
	my $syscall = modules::SystemCall->new();
	$syscall->run($cmd);
	sleep 60;
	return($outvcf);	
}

sub run {
    my ($self) = @_;

    my $command 
	= $self->executable_path 
	. ' ' . $self->call_args
	. ' -o ' . $self->output_bcf
	. ' -g ' . $self->ref_fasta;
	
	if ($self->delly_exclude) {
		$command .= ' -x '.$self->delly_exclude;
	}
	if ($self->{fast}) {
		$command .= " -q 20 -s 15";
	}  
	$command .= ' '.$self->input_bam;

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
    my ($self,$in_file) = @_;
    
    my $output_vcf;
    
    if ($in_file =~ /.bcf$/) {
    	$output_vcf = $self->_convert_to_vcf($in_file);
    } else {
    	$output_vcf = $in_file;
    }
    
    open(my $OUTPUT, $output_vcf) || modules::Exception->throw("ERROR: Can't open output file $output_vcf");
	my %sv_data = ();
	
	while (<$OUTPUT>) {
		chomp;
		next if /\tLowQual/;
		next if /^#/;
		my ($chr1,$start1,$id,$ref_base,undef,undef,undef,$gff_str) = split("\t");
		$chr1 =~ s/chr//;
		my $chr2 = $chr1;
		if ($gff_str =~ /CHR2=c?h?r?([0-9XY]+)/) {
			$chr2 = $1;
		}
		my ($start2) = $gff_str =~ /END=(\d+);/; 
		my $length = my $sr_count = my $pe_count = 0;
		if (/SR=/) {
			($sr_count) = $gff_str =~ /SR=(\d+)/;
		}
		($pe_count) = $gff_str =~ /PE=(\d+)/;
		my ($event_type) = $gff_str =~ /SVTYPE=([A-Z]+)/;
		#my $event_type = lc $event_type_uc;
		my ($qual) = $gff_str =~ /MAPQ=(\d+)/;
		$length  = abs($start2 - $start1) unless $event_type eq 'tra';
		
		if ($event_type eq 'INS') {
			$start2 = $start1;
		}
		
		
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{pe} = $pe_count;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{sr} = $sr_count;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{length} = $length;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{gff} = $gff_str;	
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{id} = $id;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{qual} = $qual;
		
		if ($event_type eq 'INS') {
			my ($seq) = $gff_str =~ /CONSENSUS=([A-Z]+)/;
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr1:$start2"}{length} = length($seq);
			$seq = 'TOO_LONG' if length($seq) > 1000;
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr1:$start2"}{seq} = '+'.$seq;
		}
		
	}

	
    return \%sv_data;
}

return 1;
