#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);
#use Bio::SearchIO;
#use Bio::SeqIO;
use modules::Exception;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "dir=s",
	   "tsv_in=s",
	   "anno_file=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || (!$OPT{dir} && !$OPT{tsv_in}) || !$OPT{anno_file});



=pod

=head1 SYNOPSIS

annotate_consensusDE.pl -tsv_in input_tsv -dir dir_of_tsvs -anno_file annotation_file_starting_with_ENSID_in_first_column [options]

Required flags: (-dir || -tsv_in) -anno_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

annotate_consensusDE.pl -> Annotate consensusDE output of DE genes

=head1 DESCRIPTION

Jan 10th, 2019

Adds all annotations in supplied file to tsv files

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

annotate_consensusDE.pl -anno_file /drive2/vardb/trunk/conf/gene_names/all_gene_info/110119/mm10_anno.gene -dir /drive3/work/paragen/output_group_filter
annotate_consensusDE.pl -anno_file /drive2/vardb/trunk/conf/gene_names/all_gene_info/110119/mm10_anno.gene -tsv_in /drive3/work/paragen/output_group_filter/DE.tsv

=cut

my @files_to_process;

if (defined $OPT{tsv_in}) {
	my $tsv_in = $OPT{tsv_in};
	
	if (!-e $tsv_in) {
		modules::Exception->throw("Can't open file $tsv_in\n");
	}
	
	push @files_to_process, $tsv_in;
}

if (defined $OPT{dir}) {
	my $dir = $OPT{dir};
	
	if (!-d $dir) {
		die("Can't find dir $dir");
	}
	
	opendir(DIR,$dir) || modules::Exception->throw("Can't open dir $dir\n");
	my @tsv_files = grep {/.tsv$/} readdir DIR;
	for my $tsv_file (@tsv_files) {
		next if $tsv_file =~ /^normalised_/;
		next if $tsv_file =~ /^raw_/;
		next if $tsv_file =~ /^sample_/;
		next if $tsv_file =~ /_annotate.tsv$/;
		push @files_to_process, "$dir/$tsv_file";
	}
}

#Lookup table for converting old header to new header names         
my %data = ();   

my $anno_file = $OPT{anno_file};

if (!-e $anno_file) {
		modules::Exception->throw("Can't open file $anno_file\n");
}
                           
open(ANNO,$anno_file) || modules::Exception->throw("Can't open $anno_file");

my @anno_headers = ();
my $first_line = 1;
my $anno_count;

while (<ANNO>) {
    chomp;
    if ($first_line) {
    	#Check there is header info (different from gene entry)
    	if ($_ !~ /ENS[A-Z]+\d+/) {
    		@anno_headers = split("\t") 
    	}
    	$first_line = 0;
    } else {
	    my @fields = split("\t");
	    $anno_count = @fields;
	    my ($ens_gene)  = $fields[0] =~ /(ENS[A-Z]+\d+)/;
	    $data{$ens_gene} = \@fields;
    }
}

my @no_anno_line = ();
push @no_anno_line, 'NO_ANNO' for (1..$anno_count); 

for my $tsv_in_file (@files_to_process) {
	open(TSV_IN,$tsv_in_file) || modules::Exception->throw("Can't open file $tsv_in_file\n");
	(my $tsv_out = $tsv_in_file) =~ s/\.tsv$/_annotate.tsv/;
	open(TSV_OUT,">$tsv_out") || modules::Exception->throw("Can't write file $tsv_out\n");
	my $header_flag = 1;
	
	while (<TSV_IN>) {
		chomp;
		my @fields = split("\t");
		
		if ($header_flag) {
			$header_flag = 0;
			my @no_quote_headers = ();
			for my $header (@fields) {
				$header =~ s/"//g;
				push @no_quote_headers, $header;
			}
			print TSV_OUT join ("\t",
						@no_quote_headers,
						@anno_headers
						) ."\n\n";
			
			
		} else {
			
			if ($fields[0] !~ /ENS[A-Z]+\d+/) {
				modules::Exception->throw("Error with format of line $_\nMust start with ENSEMBL gene name");
			}
			$fields[0] =~ s/"//g;
			$fields[0] =~ s/\.[0-9]+$//;
			my @local_anno_line = defined $data{$fields[0]}?@{$data{$fields[0]}}:@no_anno_line;
			
			print TSV_OUT join("\t",
								@fields,
								@local_anno_line
								) ."\n";
			 
		}
	}
	
}

