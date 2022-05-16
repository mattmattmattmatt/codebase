#! /usr/bin/perl -w

use strict;
use modules::Polyphen;
use modules::Adaptors::BulkInsert;
use modules::Adaptors::Run;
use modules::Adaptors::SNV;
use modules::Adaptors::Sample;
use modules::Adaptors::Filter;
use modules::VariantXML;
use modules::Report;
use modules::Overlap;
use modules::Utils;
use modules::Pipeline;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "overlap_outfile=s",
	   "tmpdir=s",
	   "runid=i",
	   "writeDB=i",
	   "annotation_file=s",
	   "exon_coord_file=s",
	   "polyphen_filter=s",
	   "ref=s"
    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{overlap_outfile} || !$OPT{tmpdir} || !$OPT{annotation_file} || !$OPT{exon_coord_file});


=pod

=head1 SYNOPSIS

filter_polyphen.pl -polyphen_filter polyphen_filter(default=filter_polyphen) -overlap_outfile <overlap_outfile> -writeDB 1|0 -tmpdir scratch_dir -annotation_file annotation_file -exon_coord_file exon_coord_file [options]

Required flags: -overlap_outfile -tmpdir -annotation_file -exon_coord_file -runid

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_polyphen.pl -> Script to drive polyphen and obtain the polyphen score

=head1 DESCRIPTION

Mar 30, 2011

Single processor polyphen wrapper

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_polyphen.pl -overlap_outfile test_polyphen/filter_polyphen.match -tmpdir test_polyphen/ -exon_coord_file /home/matt/work/pipeline/conf/mouse/NCBIM37/ccds/040411/NCBIM37.ccdsGene.txt -annotation_file /home/matt/work/pipeline/conf/mouse/NCBIM37/gene/160611/NCBIM37.BiomartMapping.txt -runid 15

=cut

# Put command line options into the right places

my $outfile = $OPT{overlap_outfile};
(my $polyphen_info_file = $outfile) =~ s/match/nomatch/;



# Put command line options into the right places
my $runid = $OPT{runid};
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $clus_conf = modules::Pipeline::get_cluster_conf();
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);
my $sample_type = modules::Pipeline::get_sample_type(-run_id=>$runid);

if ($source_type !~ /mouse/) {
	modules::Exception->throw("ERROR: Only run polyphen for mouse; for human we use variant_effect_predictor");
}

my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');


#Get the run_obj for the report
my ($run_obj) = modules::Adaptors::Run->search(id=>$runid);

my $polyphen_executable = $pipe_conf->read($source_type,'binaries','polyphen','binary');
my $polyphen_conf = $pipe_conf->read($source_type,'binaries','polyphen','polyphen_conf');
my $conf_base = $clus_conf->read($source_type,'svn','conf_dir') . '/polyphen/' . $pipe_conf->read($source_type,'binaries','polyphen','version');

#Scratch directory for the program; create if needed, otherwise reuse alignments already there
if (!-d $OPT{tmpdir}) {
	mkdir($OPT{tmpdir});
}
my $working_dir = `readlink -f $OPT{tmpdir}`;
chomp $working_dir;

if (!-d $working_dir) {
	modules::Exception->throw("Directory $working_dir doesn't exist");
} else {
    #Remove locks from previous run if they exist
    my $lockdir = $working_dir . '/lock';
    if (-d $lockdir) {
        my $lock_files = $lockdir .'/*lock';
        my $lock_flag = `ls $lock_files 2>/dev/null`;
        if ($lock_flag) {
            print STDERR "Removing locks...\n";
            system("rm $lockdir/*");
        }
    }
}

if (!-d $polyphen_conf) {
	modules::Exception->throw("Polyphen conf $polyphen_conf doesn't exist");
}

if (!-x $polyphen_executable) {
	modules::Exception->throw("Polyphen executable $polyphen_executable problems");
}

my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;
my $annotation_file = $OPT{annotation_file};
my $ref = defined $OPT{ref}?$OPT{ref}:'mm10';

my $exon_coord_file = $OPT{exon_coord_file};

my $splice_length = $pipe_conf->read($source_type,'cutoffs','splice_length_cutoff');

my $polyphen_filter_name = defined $OPT{polyphen_filter}?$OPT{polyphen_filter}:'filter_polyphen';

if ( !-e $annotation_file ) {
	modules::Exception->throw("File $annotation_file doesn't exist");	
}


if ( !-e $exon_coord_file ) {
	modules::Exception->throw("File $exon_coord_file doesn't exist");	
}






#Get the gene column name
my $config = modules::Pipeline->get_report_conf();
if (!$config->exists('snv','sample_types',$sample_type)) {
   	modules::Exception->throw("ERROR: Cannot get annotation columns for snv $sample_type");
}
my $gene_col_name;
if ($config->exists('snv','sample_types',$sample_type,'annotations')) {
	($gene_col_name) = split(",",$config->read('snv','sample_types',$sample_type,'annotations'));
} else {
	($gene_col_name) = split(",",$config->read('snv','common','annotations'));
} 

#generate the gene mapper 
my $annotation = modules::Annotation->new(-annotation_file=>$annotation_file,-exon_coord_file=>$exon_coord_file,-splice_size=>$splice_length,-sample_type=>$sample_type,-mutant_type=>'snv',-organism=>'mouse');

#my $annotation = modules::Annotation->new(-annotation_file=>$annotation_file,-exon_coord_file=>$exon_coord_file,-splice_size=>$splice_length);
#Always report for exons and splice
my $snv_report = modules::Report->new(
											 -run        => $run_obj,
                                             -gene_mapper => $annotation,
                     						 -sample_type => $sample_type,
                     						 -gene_col_name => $gene_col_name,
                     						 -source_type => $source_type
                                             );


#Load the filters, this is specific to each sample group
$snv_report->load();

$snv_report->generate_pass_fail();

#Now generate the polyphen input files
my $polyphen_input_mapped = $working_dir . '/polyphen.mapped.input.txt';
my $polyphen_input_unmapped = $working_dir. '/polyphen.unmapped.input.txt';
my $fasta = $working_dir . '/polyphen.fa';

open(my $ph_mapped, ">$polyphen_input_mapped")
    or modules::Exception->throw("Unable to open file for writing polyphen input format [$polyphen_input_mapped]");

open(my $ph_unmapped, ">$polyphen_input_unmapped")
    or modules::Exception->throw("Unable to open file for writing polyphen input format [$polyphen_input_unmapped]");

open(my $fh_fasta, ">$fasta")
    or modules::Exception->throw("Unable to open file for writing polyphen input format [$fasta]");


# Create the polyphen object for both mapped and unmapped cases
my $polyphen
    = modules::Polyphen->new(
    							-input_mapped_file    => $polyphen_input_mapped,
    							-input_unmapped_file    => $polyphen_input_unmapped,
								-executable_path => $polyphen_executable,
								-config_path => $polyphen_conf,
								-working_dir     => $working_dir,
								-fasta => $fasta					
								);
								
#Get the passed cases from the snvreport object
my ($polyphen_data,$polyphen_count) = $snv_report->generate_polyphen_input();
print STDERR "Analysing $polyphen_count entries for polyphen...\n";

#Create the overlap structure for generating the mappings
my $overlap = modules::Overlap->new();
my $mapped_found = my $unmapped_found = 0;

for my $chr ( sort keys %{$polyphen_data} ) {
	#Open an tmp file for overlapping
	my $chr_ref_file = $polyphen_input_mapped.".$chr";
	open(CHR,">$chr_ref_file") || modules::Exception->throw("Can't open file $chr_ref_file\n");
	
	for my $coord ( sort keys %{$polyphen_data->{$chr}} ) {
		#1	3062302	3062302	C G QU43432 V->L
		#1	3083416	3083416	A T QZ85904 D->S
		
		my $refbase = $polyphen_data->{$chr}{$coord}{ref}; 
		my $varbase = $polyphen_data->{$chr}{$coord}{var};
		my $uniprot = $polyphen_data->{$chr}{$coord}{uniprot};
		my $aa_change = $polyphen_data->{$chr}{$coord}{aachange};
			
		my $snv_name = "$chr:$coord:$refbase->$varbase:$aa_change:$uniprot";
		print CHR "$chr\t$coord\t$coord\t$snv_name\n";
   	}
   	close CHR;

	my $chr_mapped_file = $conf_base. '/' . $ref . '.uniprot_mapping.' . $chr;
	my $chr_unmapped_file = $conf_base. '/' . $ref . '.uniprot_nomapping.' . $chr;
	
	#First check if it's mapped
	my ($overlaps_mapped,$mapped_flag) =  $overlap->overlap(-ref=>$chr_ref_file,-coord=>$chr_mapped_file,-silent=>1);
   	if ($mapped_flag) {
   		$mapped_found  = 1;
		#Get it from the uniprot mapping file
		$polyphen->generate_input(-overlap=>$overlaps_mapped,-mapped=>1,-fh=>$ph_mapped);
   	} 

   	my ($overlaps_unmapped,$unmapped_flag) =  $overlap->overlap(-ref=>$chr_ref_file,-coord=>$chr_unmapped_file, -silent=>1);
   	if ($unmapped_flag) {
   		#Get it from the ccds file
   		$unmapped_found = 1;
   		$polyphen->generate_input(-overlap=>$overlaps_unmapped,-mapped=>0,-fh=>$ph_unmapped,-fh_fasta=>$fh_fasta);
   	}
	system("rm $chr_ref_file");
}

close($fasta);
close($ph_mapped);
close($ph_unmapped);

#Now generate the polyphen_info file; this occurs for all summary entries allowing users to run polyphen if required
open(INFO,">$polyphen_info_file") || modules::Exception->throw("Can't open file $polyphen_info_file\n");
my $row_data = $snv_report->get_all_rows();

for my $chr ( sort keys %{$row_data} ) {
	#print STDERR "CHR $chr\n";
	my $chr_info_file = $polyphen_input_mapped.".info.".$chr;
	open(CHR,">$chr_info_file") || modules::Exception->throw("Can't open file $chr_info_file\n");
    for my $coord (sort {$a<=>$b} keys %{$row_data->{$chr}}) {
    	print CHR "$chr $coord $coord\n";	
    }
     
    close CHR;
    
    my $chr_mapped_file = $conf_base. '/' . $ref . '.uniprot_mapping.' . $chr;
	my $chr_unmapped_file = $conf_base. '/' . $ref . '.uniprot_nomapping.' . $chr;
	
	my ($overlaps_mapped,$mapped_flag) =  $overlap->overlap(-coord=>$chr_info_file,-ref=>$chr_mapped_file,-silent=>1);
   	if ($mapped_flag) {
   		for my $chr (keys %{$overlaps_mapped->{PASS}}) {
			for my $overlap_coord (keys %{$overlaps_mapped->{PASS}{$chr}}) {
				for my $snv_str (keys %{$overlaps_mapped->{PASS}{$chr}{$overlap_coord}{$overlap_coord}}) {
					print INFO join("\t",
									$chr,
									$overlap_coord,
									$overlap_coord,
									$snv_str
									) . "\n";
				}
			}
   		}
   	} 

   	my ($overlaps_unmapped,$unmapped_flag) =  $overlap->overlap(-coord=>$chr_info_file,-ref=>$chr_unmapped_file, -silent=>1);
   	if ($unmapped_flag) {
   		for my $chr (keys %{$overlaps_unmapped->{PASS}}) {
			for my $overlap_coord (keys %{$overlaps_unmapped->{PASS}{$chr}}) {
				for my $snv_str (keys %{$overlaps_unmapped->{PASS}{$chr}{$overlap_coord}{$overlap_coord}}) {
					print INFO join("\t",
									$chr,
									$overlap_coord,
									$overlap_coord,
									$snv_str
									) . "\n";
				}
			}
   		}
   	}
    
    system("rm $chr_info_file");
    
}

#The two output files
my $polyphen_output_mapped = $working_dir . '/polyphen.mapped.output.txt';
my $polyphen_output_unmapped = $working_dir . '/polyphen.unmapped.output.txt';

if ($mapped_found) {
	$polyphen->run(-output=>$polyphen_output_mapped,-mapped=>1);
	$polyphen->parse_result(-output=>$polyphen_output_mapped);
}

if ($unmapped_found) {
	$polyphen->run(-output=>$polyphen_output_unmapped,-mapped=>0);
	$polyphen->parse_result(-output=>$polyphen_output_unmapped);
}

open(FILE,">$outfile") || modules::Exception->throw("Can't open file to write $outfile\n");

my $combined_results = $polyphen->get_results();

my @polyphen_inserts = ();

my %mutant_snv_data = ();

for my $chr (sort keys %{$combined_results}) {
	
	my %chr_snvs;

	for my $coord (sort {$a<=>$b} keys %{$combined_results->{$chr}}) {
		for my $nt_change (keys %{$combined_results->{$chr}{$coord}}) {
			my @fields = split(':',$combined_results->{$chr}{$coord}{$nt_change}{name});
			
			my $attribute = 'poly_pred=' . $combined_results->{$chr}{$coord}{$nt_change}{prediction} . ';poly_score=' . $combined_results->{$chr}{$coord}{$nt_change}{score};
			
		
			my $mutant_snv_key =  $chr.":".$coord;
			my ($var_base) = $nt_change =~ /->(.*)/;
			$mutant_snv_data{$mutant_snv_key}{$coord}{$var_base}{snv_filter_polyphen_string} = $attribute;
			$mutant_snv_data{$mutant_snv_key}{$coord}{$var_base}{pass} = 1;
			
			print FILE join("\t",
								$chr,
								$coord,
								$coord,
								$combined_results->{$chr}{$coord}{$nt_change}{score},
								$combined_results->{$chr}{$coord}{$nt_change}{prediction},
								$combined_results->{$chr}{$coord}{$nt_change}{name}
								). "\n";
			
		}
	}
}

close FILE;

my $file_name = $sample_name . '_' . $runid . '.filter_polyphen.snv.xml';
$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_snv_data, -chr=>'all');
my $full_xml_file = $run_dir . '/conf/'.$file_name;
$var_xml->split_xml(-file_name=>$full_xml_file);










