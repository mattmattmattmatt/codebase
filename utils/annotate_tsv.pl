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
	   "anno_file=s",
	   "anno_lookup_colnum=i",
	   "file_lookup_colnum=i"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || (!$OPT{dir} && !$OPT{tsv_in}));



=pod

=head1 SYNOPSIS

annotate_tsv.pl -tsv_in input_tsv -dir dir_of_tsvs -anno_file annotation_file -anno_lookup_colnum default=0 -file_lookup_colnum default=0 [options]

Required flags: (-dir || -tsv_in) 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

annotate_tsv.pl -> Annotate tsv files to incorporate gene files

=head1 DESCRIPTION

Jan 10th, 2019

Adds all annotations in supplied file to tsv files

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

annotate_tsv.pl -anno_file /drive2/vardb/trunk/conf/gene_names/all_gene_info/110119/mm10_anno.gene -dir /drive3/work/paragen/output_group_filter
annotate_tsv.pl -tsv_in pre_post_topTable.csv -anno_file /drive2/codebase/utils/GRCh38.tsv -anno_lookup_colnum 2 -file_lookup_colnum 0

=cut

my @files_to_process;

my $anno_colnum = defined $OPT{anno_lookup_colnum}?$OPT{anno_lookup_colnum}:0;
my $file_colnum = defined $OPT{file_lookup_colnum}?$OPT{file_lookup_colnum}:0;

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
	my @tsv_files = grep {/.[ct]sv$/} readdir DIR;
	for my $tsv_file (@tsv_files) {
		next if $tsv_file =~ /_annotate.[ct]sv$/;
		push @files_to_process, "$dir/$tsv_file";
	}
}

#Lookup table for converting old header to new header names         
my %data = ();   
my $svndir = $ENV{'SVNDIR'}; 

my $anno_file = defined $OPT{anno_file}?$OPT{anno_file}:"${svndir}/utils/GRCh38.tsv";

if (!-e $anno_file) {
		modules::Exception->throw("Can't open file $anno_file\n");
}
                           
open(ANNO,$anno_file) || modules::Exception->throw("Can't open $anno_file");

my @anno_headers = ();
my $first_line = 1;
my $anno_count;
my $gene_header;

while (<ANNO>) {
    chomp;
    if ($first_line) {
    	#Check there is header info (different from gene entry)
    	if ($_ =~ /Gene/i) {
    		@anno_headers = split("\t") 
    	}
#    	for ( my $header_count = 0 ; $header_count < @anno_headers ; $header_count++ ) {
#    	 	if ($anno_headers[$header_count] =~ /Gene$/) {
#    	 		$gene_header = $header_count;
#    	 		
#    	 	}
#    	 	   
#    	}
    	$first_line = 0;
    } else {
	    my @fields = split("\t");
	    
	    $anno_count = @fields;
	    my $gene;
	    if ($anno_headers[$anno_colnum] =~ /URL/) {
			($gene) = $fields[$anno_colnum] =~ /(ENS[A-Z0-9]+)/;	    	
	    } else {
		    ($gene)  = $fields[$anno_colnum] =~ /(\S+)/;
	    }
	    $data{uc($gene)} = \@fields;
    }
}




my @no_anno_line = ();
push @no_anno_line, 'NO_ANNO' for (1..$anno_count); 

my $found = my $not_found = 0;

for my $tsv_in_file (@files_to_process) {
	open(TSV_IN,$tsv_in_file) || modules::Exception->throw("Can't open file $tsv_in_file\n");
	my ($base,$suffix) = split('\.',$tsv_in_file);
	my $tsv_out = $base.'_annotate.'.$suffix;
	open(TSV_OUT,">$tsv_out") || modules::Exception->throw("Can't write file $tsv_out\n");
	my $header_flag = 1;
	
	while (<TSV_IN>) {
		chomp;
		my @fields;
		if ($suffix eq 'csv') {
			@fields = split(",");
		} else {
			@fields = split("\t");
		}
		if (@fields == 1) {
	    	modules::Exception->throw("ERROR: File $tsv_in_file doesn't split on tab delim\n");
	    }
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
			
			
			my $lookup_gene = $fields[$file_colnum];
			$lookup_gene =~ s/"//g;
			
			#print "LOOK $fields[$file_colnum]\n";
			
			if (defined $data{uc($lookup_gene)}) {
				$found++;
				print TSV_OUT join ("\t",
									@fields,
									@{$data{uc($lookup_gene)}}
										) ."\n"; 
			} else {
				$not_found++;
				print TSV_OUT join ("\t",
									@fields,
									@no_anno_line
									) ."\n";
			}
				
			
			 
		}
	}
	close TSV_IN;
	print "Matched: $found lines\nNot matched: $not_found lines\n";	
}


