#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use modules::Utils;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
			"master_file=s",
			"add_file=s",	   		
			"out=s",
	   		"debug=s",
	   		);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{master_file} || !$OPT{add_file});



=pod

=head1 SYNOPSIS

combine_files.pl -master_file master_input_file -add_file file_to_add_to_novel_entries_to_master -out outfile(default=master.combined) -debug [options]

Required flags: -master_file -add_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

combine_polyphen.pl -> combine files of genomic coordinates with a master and add file.  Master entry is maintained when both files have same coordinate

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./combine_files.pl -master_file /drive2/variantdb/trunk/conf/human/GRCh37/dbsnp/138/GRCh37.dbsnp.overlap.snv.1 -add_file /drive2/variantdb/trunk/conf/human/GRCh37/dbsnp/137/GRCh37.dbsnp.overlap.snv.1 -out /drive2/variantdb/trunk/conf/human/GRCh37/dbsnp/137_138/GRCh37.dbsnp.overlap.snv.1    

=cut

my $master_infile = $OPT{master_file};
my $add_infile = $OPT{add_file};
my $outfile = defined $OPT{out}?$OPT{out}:$master_infile.'.combine';




my %data;

open (ADD, $add_infile) || modules::Exception->throw("Cannot open file $add_infile");

while (<ADD>) {
	chomp;
	my @fields = split(" ");
	push @{$data{$fields[0]}{$fields[1]}{$fields[2]}{add}}, $fields[3]; 
}
close ADD;

open (MASTER, $master_infile) || modules::Exception->throw("Cannot open file $master_infile");

while (<MASTER>) {
	chomp;
	my @fields = split(" ");
	push @{$data{$fields[0]}{$fields[1]}{$fields[2]}{master}}, $fields[3];
}
close MASTER;

open(OUT,">$outfile") || die "Can't open out $outfile for writing\n";

for my $chr (sort keys %data) {
	for my $start_coord (sort {$a<=>$b} keys %{$data{$chr}}) {
		for my $end_coord (sort {$a<=>$b} keys %{$data{$chr}{$start_coord}}) {
			if (exists $data{$chr}{$start_coord}{$end_coord}{master}) {
				for my $master_data ( @{$data{$chr}{$start_coord}{$end_coord}{master}} ) {
					print OUT join(" ",
										$chr,
										$start_coord,
										$end_coord,
										$master_data
										). "\n";
				    
				}
			} else {
				for my $add_data ( @{$data{$chr}{$start_coord}{$end_coord}{add}} ) {
					print OUT join(" ",
										$chr,
										$start_coord,
										$end_coord,
										$add_data
										). "\n";
				    
				}
			}
		}
	}
}

