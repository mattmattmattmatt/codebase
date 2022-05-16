package modules::PED;

use strict;
use modules::Exception;
use Data::Dumper;
use modules::PipelineSample;

sub new {
    my ($class) = @_;

    my $self = bless {}, $class;

    return $self;
}

#Create a ped file a source input line; needed for production pipeline
sub create_ped {
	my ($self, @args) = @_;

    my @required_args = (
			             -source_line,
			             -output
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    chomp;
	if (my $error = modules::PipelineSample->validate_source_line(-line=>$_)) {
		modules::Exception->throw("ERROR $error with line $_");
	}
	my @fields = split(",");
    
    
}

#Parse a ped file
sub parse_ped {
    my ($self, @args) = @_;

    my @required_args = (
			             -ped_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-ped_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    open(PED,$args{-ped_file}) || modules::Exception->throw("Can't open file $args{-ped_file}\n");
    
    my %ped_data;
    
    while (<PED>) {
    	chomp;
    	my @fields = split;
    	if (@fields < 6) {
    		modules::Exception->throw("ERROR: Expecting Ped file with at least six fields (family_id, id, father, mother, sex, and affected");
    	}
    	my ($family,$id,$father,$mother,$sex,$affected) = split @fields;
    	
    	if ($affected !~ /[12]/) {
    		modules::Exception->throw("ERROR: Affected status must be 1 or 2");
    	}
    	if ($sex !~ /[12]/) {
    		modules::Exception->throw("ERROR: Sex must be 1 or 2");
    	}
    	
    	$ped_data{$family}{$id}{father} = $father;
    	$ped_data{$family}{$id}{mother} = $mother;
    	
    	if ($affected == 2) {
    		$ped_data{$family}{$id}{affected} = 1;
    	} else {
    		$ped_data{$family}{$id}{affected} = 0;
    	}
    	
    	if ($sex == 1) {
    		$ped_data{$family}{$id}{sex} = 'male';
    	} else {
    		$ped_data{$family}{$id}{sex} = 'female';
    	}
    	
    	
    }
    
    $self->{data}{$args{-ped_file}} = \%ped_data;
}

#Get ped data from a file
sub get_ped {
	my ($self, @args) = @_;

    my @required_args = (
			             -ped_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-ped_file} ) {
    	modules::Exception->throw("File $self->{data}{$args{-ped_file}} doesn't exist");	
    }

   return $self->{data}{$args{-ped_file}};

}


1;