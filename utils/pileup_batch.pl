#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use modules::Utils;
use vars qw(%OPT);
use modules::PileupLine;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"pileup_file=s",
	   	"min_freq=s",
	   	"min_reads=i",
	   	);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{pileup_file});



=pod

=head1 SYNOPSIS

pileup_batch.pl -help -man [options] -pileup pileup_file_from_samtools -min_freq min_var_freq_to_include -min_reads min_reads_for_variant_to_include

Required flags: -pileup_file pileup_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

FILE.pl 

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./FILE.pl    

=cut

my $pileup_file = $OPT{pileup_file};

if ( !-e $pileup_file ) {
	modules::Exception->throw("File $pileup_file doesn't exist");
}


my $min_freq = defined $OPT{min_freq}?$OPT{min_freq}:'0';
my $min_reads = defined $OPT{min_reads}?$OPT{min_reads}:'0';

open(PILEUP,"$pileup_file") || modules::Exception->throw("Can't open file $pileup_file\n");

while (<PILEUP>) {
    chomp;
    my @fields = split(" ");
    my $pl = modules::PileupLine->new();
    $pl->ref_base($fields[2]);
	$pl->base_string($fields[4]);
	
	if ($fields[3] == 0) {
		#print join ("\t",$fields[0],$fields[1],$fields[3],"NO_READS")."\n";
		next; 
	}
	
	my $ref_base = uc($fields[2]);
	my $zyg = 'het';
	my %freq = ();
	if (defined %{$pl->base_frequencies->lookup}{$ref_base}) {
		$freq{$ref_base."(ref)"} = %{$pl->base_frequencies->lookup}{$ref_base};
	} else {
		$freq{$ref_base."(ref)"} = '0.00000';
		$zyg = 'hom';
	}
	
	my $max_alt = 0;
	my $max_reads = 1;
		
	for my $base (keys %{$pl->base_frequencies->lookup}) {
		next if $base eq $ref_base;
		$freq{$base} = %{$pl->base_frequencies->lookup}{$base};
		if ($freq{$base} > $max_alt && $base ne 'REFSKIP') {
			$max_alt = $freq{$base};
			$max_reads = int($fields[3]*$max_alt);
		}
			
	}

	next if ($max_alt < $min_freq);
	next if ($max_reads < $min_reads);
	#print Dumper $pl;
	

	my $base_str;
 	for my $base ( sort {$freq{$b}<=>$freq{$a}} keys %freq ) {
 	    $base_str .= $base.":".sprintf("%.5f",$freq{$base}).';';
 	}
 	$base_str =~ s/;$//;
	print join ("\t",$fields[0],$fields[1],$fields[3],$fields[4], $base_str, $zyg) . "\n";
	
	#my @counts = @{$pl->base_frequencies->counts};
	#my @freqs = @{$pl->base_frequencies->counts};
	
	
}
