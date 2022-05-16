#! /usr/bin/perl -w

use strict;
use modules::Report;
use modules::Adaptors::Filter;
use modules::Adaptors::BulkInsert;
use modules::VariantXML;
use modules::Annotation;
use modules::Pipeline;
use modules::ConfigXML;
use Getopt::Long;
use modules::Exception;
use modules::Utils;
use File::Basename;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

GetOptions(\%OPT, 
		    "help|h",
		    "man|m",
		 	"runid=i",
		    "chr=s",
		 	"start=i",
		 	"end=i",
		 	"ref=s",
		 	"filter_summary_file=s",
		 	"annotation_file=s",
		 	"exon_coord_file=s",
		    "tsv_file=s",
		    "pass_summary_file=s",
		    "pass_file=s",
		    "genefilter_files=s",
		    "polyphen_file=s",
		    "pass_filter=s",
		    "writeDB=i",
			"report_xml=s",
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{tsv_file} || !$OPT{runid} || !$OPT{filter_summary_file} || !$OPT{annotation_file} || !$OPT{exon_coord_file} || !$OPT{pass_summary_file} || !$OPT{pass_file});

	   
=pod

=head1 SYNOPSIS

generate_snv_eport.pl -runid runid -polyphen_file polyphen_score_file -pass_filter filter_report_name(default=filter_pass_snv) -ref ref_name(default=NCBIM37) -annotation_file annotation_file -exon_coord_file exon_coord_file -debug print_stderr -filter_summary_file filter_summary_output_file -tsv_file output_tsv_file -chr specific_chr(default=all) -start chr_start_coord(default=all) -end chr_end_coord(default=all) -pass_summary_file pass_rows_summary [options]

Required flags: -runid -tsv_file -filter_summary_file -annotation_file -exon_coord_file -pass_summary_file -pass_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

generate_snv_report.pl -> Script to generate a snv report

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Dan Andrews

=head1 EXAMPLE

generate_snv_report.pl

=cut

my $pipe_conf = modules::Pipeline::get_pipe_conf();

#Lots of variables to load...
my $output_file = $OPT{output_file};
my $summary_file = $OPT{pass_summary_file};
my $pass_file = $OPT{pass_file};
my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:1;
my $pass_filter = defined $OPT{pass_filter}?$OPT{pass_filter}:'filter_pass_snv';

my $runid = $OPT{runid};

my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);






my @filter_files = ();
if (defined $OPT{genefilter_files}) {
	@filter_files = split(",",$OPT{genefilter_files});
}

#Flag tells us whether to report polyphen values
my $polyphen_score = my $polyphen_info = 0;
my $polyphen_score_file;
my $polyphen_info_file;
if (defined $OPT{polyphen_file}) {
	$polyphen_score_file = $OPT{polyphen_file};
	($polyphen_info_file = $polyphen_score_file) =~ s/match/nomatch/;
	#One of the two must exist...
	if (!-e $polyphen_score_file && !-e $polyphen_info_file) {
		modules::Exception->throw("ERROR: Polyphen files $polyphen_score_file and $polyphen_info_file don't exist");
	} 
	
	if (-e $polyphen_score_file) {
		$polyphen_score = 1;
	}
	
	if (-e $polyphen_info_file) {
		$polyphen_info = 1;
	}
	
} 

my $tsv_file = $OPT{tsv_file};
my $chr = defined $OPT{chr}?$OPT{chr}:'';
my $start = defined $OPT{start}?$OPT{start}:'';
my $end = defined $OPT{end}?$OPT{end}:'';
my $annotation_file = $OPT{annotation_file};
my $ref = defined $OPT{ref}?$OPT{ref}:'GRCh37';
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);
my $sample_type = modules::Pipeline::get_sample_type(-run_id=>$runid);
my $organism = $pipe_conf->read($source_type,'organism');

my $splice_size = $pipe_conf->read($source_type,'cutoffs','splice_length_cutoff');
my $exon_coord_file = $OPT{exon_coord_file};

if ( !-e $annotation_file ) {
	modules::Exception->throw("File $annotation_file doesn't exist");	
}
	

if ( !-e $exon_coord_file ) {
	modules::Exception->throw("File $exon_coord_file doesn't exist");	
}

my ($run_obj) = modules::Adaptors::Run->search('id' => $runid);
unless (defined $run_obj) {
    modules::Exception->throw("Did not get all necessary things from database [run_id:$runid]");
}

#Get the gene column name
my $config = defined($OPT{report_xml}) ? modules::Pipeline->get_report_conf(-report_xml=>$OPT{report_xml}) : modules::Pipeline->get_report_conf();
if (!$config->exists('snv','sample_types',$sample_type)) {
   	modules::Exception->throw("ERROR: Cannot get annotation columns for snv $sample_type");
}
 
my $gene_col_name;
if ($organism eq 'human') {
	($gene_col_name) = split(",",$config->read('snv','common','human_file_annotations'));
} else {
	($gene_col_name) = split(",",$config->read('snv','common','mouse_file_annotations'));
}

    
#itenerate the gene mapper 
my $annotation = modules::Annotation->new(-annotation_file=>$annotation_file,-exon_coord_file=>$exon_coord_file,-splice_size=>$splice_size,-sample_type=>$sample_type,-mutant_type=>'snv',-organism => $organism, -config => $config);

#Create the snv report params

my %snv_report_param = ( -run => $run_obj,
	                     -gene_mapper => $annotation,
	                     -sample_type => $sample_type,
	                     -source_type => $source_type,
	                     -gene_col_name => $gene_col_name
	                     );


if (@filter_files) {
	$snv_report_param{-filter_files} = \@filter_files
}

if ($polyphen_score) {
	$snv_report_param{-polyphen_file} = $polyphen_score_file;
}

if ($polyphen_info) {
	$snv_report_param{-polyphen_info_file} = $polyphen_info_file;
}

my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);

$snv_report_param{-confdir} = $run_dir . '/conf';

my $snv_report = modules::Report->new(%snv_report_param);


#Load the filters, this is specific to each sample group
if ($OPT{report_xml}){
	$snv_report->load(-report_xml=>$OPT{report_xml});
} else {
	$snv_report->load();
}

if ($chr && $start) {
	$snv_report->generate_pass_fail(-chr=>$chr, -start=>$start, -end=>$end);
} elsif ($chr) {
	$snv_report->generate_pass_fail(-chr=>$chr);
} else {
	$snv_report->generate_pass_fail();
}

$snv_report->print_to_files(-tsv_file=>$tsv_file,-pass_file=>$pass_file);
$snv_report->summarize(-summary_file=>$summary_file);

my $overlap = $OPT{filter_summary_file};
my $overlapdir = dirname($overlap);
$overlapdir =~ s/summary//;
opendir(DIR,$overlapdir) || die "Can't open directory $overlapdir\n";
my @local_filter_files = grep {/filter_/} readdir DIR;
closedir (DIR);

my %all_snvs = ();
my %bp_lookup = ();

for my $filter_file ( @local_filter_files ) {
	next if $filter_file =~ /indel/;
	next if $filter_file =~ /^vep/;
	#This filter is covered by filter_basic
    open(FILTER,"$overlapdir/$filter_file") || modules::Exception->throw("Can't open file $overlapdir/$filter_file\n");
    while (<FILTER>) {
        my ($chr,$start,$end,$rest) = split;
        if ($rest =~ /([AGTCN]->[ATGC])/) {
        	$bp_lookup{$chr}{$start} = $1;
        }
        my ($filter) = $filter_file =~ /filter_(.*)$/;
        $filter =~ s/\.snv//;
        push @{$all_snvs{$chr}{$start}},$filter;
    }
    close FILTER;
}

my @overlap_lines = ();
open(OVERLAP,">$overlap") || die "Can't open file to write $overlap\n";
for my $chr (sort keys %all_snvs) {
	for my $coord ( sort {$a<=>$b} keys %{$all_snvs{$chr}} ) {
		my $filter_str = join(",",sort @{$all_snvs{$chr}{$coord}});
		if (exists $bp_lookup{$chr}{$coord}) {
			my $bp_change = $bp_lookup{$chr}{$coord};		
		    print OVERLAP "$chr $coord $coord $bp_change $filter_str\n";
		}
	}
}

close OVERLAP;

#Get the snv filter
my @snv_filter_inserts;


#First the passed rows
my $pass_row_data = $snv_report->get_pass_rows();
my $rare_allele_rows = $snv_report->get_allele_rows();
my $no_freq_rows = $snv_report->get_no_freq_rows();
my $fail_row_data = $snv_report->get_fail_rows();


my %rows = (
			'pass'=>$pass_row_data,
			'allele'=>$rare_allele_rows,
			'no_freq'=>$no_freq_rows,
			'fail'=>$fail_row_data
			);

my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');


for my $row_type (keys %rows) {
	my $row_count = 0;
	for my $chr ( sort keys %{$rows{$row_type}} ) {
			my $file_name = $sample_name . '_' . $runid . '.db_variants.snv.'.$chr.'.xml';
			$var_xml->load_xml_keys(-file_name=>$file_name);
		    for my $coord (sort {$a<=>$b} keys %{$rows{$row_type}->{$chr}}) {
		    	
	    		if (!$var_xml->search_xml_by_coord(-var_type=>'snv',-chrom=>$chr,-start_coord=>$coord,-end_coord=>$coord)) {
	    			modules::Exception->throw("ERROR: Can't retrieve snv with above search parameters");
	    		}
		    	
				$row_count++;
		    }
	}
	print STDERR "Reporting $row_count $row_type rows\n";
}

my $final_xml = $sample_name . '_' . $runid . '.report_snv.snv.xml';
$snv_report->final_xml(-xml_file=>$final_xml);


