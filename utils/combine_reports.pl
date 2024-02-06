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
	   "report1=s",
	   "report2=s",
	   "combine_file=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{report1} || !$OPT{report2});



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


my %data = ();

my $report1 = $OPT{report1};
my $report2 = $OPT{report2};



my $combine = defined $OPT{combine_file}?$OPT{combine_file}:"combined.tsv"; 

if ( !-e $report1 || !-e $report2) {
	modules::Exception->throw("File $report1 or $report2 doesn't exist");
}

open(COMBINE,">$combine") || modules::Exception->throw("Can't open file $combine\n");

for my $file ($report1 $report2) {
	open(FILE,"$file") || modules::Exception->throw("Can't open file $file\n");
	
	while (<FILE>) {
		chomp;
		next if /^chr/;
		next if /^$/;
		my @fields = split("\t");
		my $key = "$fields[0]:$fields[1]:$fields[2]:$fields[13]";
		$data{$key}{$file} = "$fields[7]:$fields[8]:$fields[9]";
	}
	
}

while (<COMBINE>) {
	chomp;
	my @fields = split("\t");
	next if /^$/;
	if (/^chr/) {
		print join("\t",
					@fields[0..6],
					'Day24_sc_nodata',
					'Day24_sc_ref',
					'Day24_sc_variants (singlecell_variant %)',
					'Day548_sc_nodata',
					'Day548_sc_ref',
					'Day548_sc_variants (singlecell_variant %)',
					'Product_sc_nodata',
					'Product_sc_ref',
					'Product_sc_variants (singlecell_variant %)',
					@fields[7..40]
					) ."\n\n";
	} else {
		my $key = "$fields[0]:$fields[1]:$fields[2]:$fields[10]";
		
		my @sc_info = ();
		
		for my $sample (@samples) {
			if (exists $data{$sample}{$key}) {
				my @sc_fields = split(':',$data{$sample}{$key});
				push @sc_info, $sc_fields[0],$sc_fields[1],$sc_fields[2];
			} else {
				push @sc_info, "N/A","N/A","N/A";
			}
			
		}
		print join("\t",
					@fields[0..6],
					@sc_info,
					@fields[7..40]
					) . "\n";
		
		
		
		
		
		
	}
	
	
	
}







