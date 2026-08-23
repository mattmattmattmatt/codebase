#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "report=s",
	   "outfile=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{report});



=pod

=head1 SYNOPSIS

script_name [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

script_name.pl -> One line description

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut



my $report = $OPT{report};
my @samples = ();
my $sample_flag =  0;
my %headers = ();
my %header_names = ();

my %key_count = ();


if ( !-e $report) {
	modules::Exception->throw("File $report doesn't exist");
}


my $out;
if (defined $OPT{outfile}) {
	$out = $OPT{outfile};
} else {
	($out = $report) =~ s/.tsv/_fixsamples.tsv/;
}

open(OUT,">$out") || modules::Exception->throw("Can't open file $out\n");

open(FILE,$report) || modules::Exception->throw("Can't open file $report\n");

while (<FILE>) {
	chomp;
	
	my @fields = split("\t");
	my $count = 0;
	if (/^chr/) {
		
		for my $field (@fields) {
			if ($sample_flag == 1) {
				push @samples, $field;
			}
			
			if ($field eq 'mouse_phenotype') {
				$sample_flag = 1;
			}
			$header_names{$field} = $count;
			$headers{$count} = $field;
			$count++;
			
			
			
		}
		print OUT join("\t",@fields) ."\n\n";
		next;
	}
	if (!@samples) {
		modules::Exception->throw("ERROR: Need 'mouse_phenotype' column to set sample flag $sample_flag\n");
	}
	next if /^$/;

	
	for my $field (@fields) {
		if ($headers{$count} eq 'variant_read_count') {
			my $no_data = 0;
			my $var = 0;
			my $het = 0;
			my $hom = 0;
			my $ref = 0;
			my @read_counts = split(",",$field);
			my $read_count = 0;
			my @variant_samples = ();
			for my $sample ( @samples ) {
			    my ($var_base,$total_base) = split('/',$read_counts[$read_count]);
			    if ($var_base == 0 && $total_base == 0) {
			    	$no_data++;
			    	$fields[$header_names{$sample}] = 'no_data';
			    } elsif ($var_base == 0) {
			    	$ref++;
			    	$fields[$header_names{$sample}] = 'ref';
			    } elsif ($var_base > 0) {
			    	$var++;
			    	push @variant_samples, $sample;
			    	my $af = sprintf("%.4f", $var_base/$total_base);
			    	if ($af < 0.7) {
			    		$het++;
			    		$fields[$header_names{$sample}] = 'het';
			    	} else {
			    		$hom++;
			    		$fields[$header_names{$sample}] = 'hom';
			    	}
			    }
			    $read_count++;
			}
			$fields[$header_names{'variant_count (het/hom)'}] = $var .' ('.$het.'/'.$hom.')';
			$fields[$header_names{'ref_count'}] = $ref;
			$fields[$header_names{'no_data_count'}] = $no_data;
			$fields[$header_names{'variant_samples'}] = join(",",@variant_samples);
			
		}
		$count++;
	}
	print OUT join("\t", @fields) ."\n";

}
	



