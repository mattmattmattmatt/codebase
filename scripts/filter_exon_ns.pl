#! /usr/bin/perl -w

use strict;
use modules::VEP;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::BulkInsert;
use modules::VariantXML;
use modules::SystemCall;
use modules::Pipeline;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Env qw($ENSEMBL_REGISTRY);
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "overlap_outfile_snv=s",
		   "ref_file_snv=s",
		   "exonns_filter_name=s",
		   "tmpdir=s",
		   "runid=i",
		   "writeDB=i",
		   "chr=s"
	    );

	   
pod2usage(1) if ($OPT{help} || !$OPT{overlap_outfile_snv} || !$OPT{ref_file_snv} || !$OPT{runid} || !$OPT{tmpdir});


=pod

=head1 SYNOPSIS

filter_exon_ns.pl -ref_file_snv <overlap_infile_snv> -overlap_outfile_snv <overlap_outfile_snv> -chr chr_to_run_on -exonns_filter_name <filter_name for non-synonymous variants, eg. filter_exon_ns> -writeDB 1|0

Required flags: -runid -overlap_outfile_snv -ref_file_snv -tmpdir

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

filter_exon_ns.pl -> Script to drive variant_effect_predictor.pl and parse details of exonic and non-synonymous variants for snvs and indels

=head1 DESCRIPTION

Mar 30, 2011

It's ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

./filter_exon_ns.pl

=cut

# Put command line options into the right places
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my $clus_conf = modules::Pipeline::get_cluster_conf();
my $runid = $OPT{runid};
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);


my $exonns_file = $OPT{overlap_outfile_snv};
#Need exon file to get ensembl codon information
(my $exon_file = $exonns_file) =~ s/_ns//;

my $ref_file = $OPT{ref_file_snv};
my $vep_executable = $pipe_conf->read($source_type,'binaries','variant_predictor','binary');
my $vep_db_dir = $clus_conf->read($source_type,'svn','conf_dir') . '/vep_index';
my $working_dir = $OPT{tmpdir};
my $exonns_filter_name = defined $OPT{exonns_filter_name}?$OPT{exonns_filter_name}:'filter_exon_ns';
my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;
my $organism = $pipe_conf->read($source_type,'organism');

my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);



my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	print "ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n";
	exit;
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $vep_input_file = $working_dir . '/vep.filter_exon_ns.in';
my $vep_output_file = $working_dir . '/vep.filter_exon_ns.out';

open(my $DBVAR_FILE, $ref_file)
    or modules::Exception->throw("Unable to open input file [$ref_file]");

if (!-e $exon_file) {
	print STDERR "No exon overlap entries\n";
	exit;
}

open(my $EXON_FILE, $exon_file)
    or modules::Exception->throw("Unable to open input file [$exon_file]");
        
    
open(my $VEP_INPUT_FILE, ">$vep_input_file")
    or modules::Exception->throw("Unable to open file for writing annovar input format [$vep_input_file]");

#Connect to database

my ($run) = modules::Adaptors::Run->search('id' => $runid);

my %vep_input_info = ();

#Build up the string for input from the db_variants file
while (<$DBVAR_FILE>){
    chomp;

	#1	55776	55776	SNV^^^C->T^^^Q202^^^TUM81^^^D34^^^het^^^16/34^^^MQ40^^^TBS:,TTt,,,.,T,,.,,tT,.T.Ttt.TtTT.tT,,NBS:.,.,.,.,.,..,,,,,....,,,.,.^^^41:40:41:41:40:41:41:41:41:41:34:34:41:41:40:38:34:39:41:41:41:38:32:37:41:41:36:39:39:37:28:26:34:35

    my ($chr, $start, $end, $attributes) = split /\t/;
	my (undef,$base_change) = split('\^\^\^',$attributes);
	my ($ref_base,$var_base) = $base_change =~ /(.*)->(.*)/;


    if (!defined $chr || ! defined $start || ! defined $end || ! defined $attributes){
		modules::Exception->warning("Did not parse input line successfully [$_]");
    }

    unless ($var_base =~ /[ATGC]/) {
		modules::Exception->warning("Did not parse attributes successfully [$attributes]");
    }

	my $base_str = $ref_base . '/' . $var_base;
	$vep_input_info{$chr}{$start} = join("\t", $chr, $start, $end, $base_str, '+' );
    #print $VEP_INPUT_FILE join("\t", $chr, $start, $end, $base_str, '+' ) . "\n";
}

close($DBVAR_FILE);

#Now get the peptide info from ensembl to see if we need to combine any variants
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($ENSEMBL_REGISTRY);
my $slice_adaptor = $registry->get_adaptor($organism, 'core', 'Slice');
# get the adaptors
my $gene_adaptor = $registry->get_adaptor($organism, "Core", "Gene");

my %pep_to_coord = ();
my %coord_to_pep = ();
my %aa_length = ();

#filter_exon file -> map the coordinates to peptide location for combining adjacent variants
while (<$EXON_FILE>) {
	chomp;
	#1	1269831	1269831	SNV^^^G->A^^^SOMATIC^^^Q205^^^REF_NORM99^^^D35^^^het^^^15/35^^^MQ36^^^TBS:,.A^^^NBS:,,...,..,,,,...^^^39:30:34	ENSG00000169962_exon6
	my @fields = split("\t");
	next unless /ENS/;
	my $ens_gene_name;
	if ($fields[4] =~ /(ENSG\d+)/) {
		$ens_gene_name = $1;
	} elsif ($fields[4] =~ /(ENSMUSG\d+)/) {
		$ens_gene_name =$1;
	} else {
		modules::Exception->throw("ERROR: Can't get ensembl gene name from $fields[4]");
	}
	my $gene = $gene_adaptor->fetch_by_stable_id($ens_gene_name);
	next unless $gene;
	my $canonical_transcript = $gene->canonical_transcript();
	my $transcript_length = 0;
	my @exon_objs = @{$canonical_transcript->get_all_Exons()};
    foreach my $exon (@exon_objs) {
    	if ($exon->coding_region_start($canonical_transcript)) {
    		$transcript_length += abs( $exon->coding_region_end($canonical_transcript) - $exon->coding_region_start($canonical_transcript)) + 1;
    	}
    }
	$aa_length{$ens_gene_name} = $transcript_length/3 - 1; # Account for stop codon
	my $strand = $canonical_transcript->strand();
	my $transcriptMapper = $canonical_transcript->get_TranscriptMapper();
	my @pepCoords = $transcriptMapper->genomic2pep($fields[1], $fields[1], $strand);
	foreach my $pepinfo (@pepCoords){
		push @{$pep_to_coord{$ens_gene_name . ':' . $pepinfo->start}}, "$fields[0]:$fields[1]"; #combine the entries
		$coord_to_pep{$fields[0]}{$fields[1]} = $ens_gene_name . ':' . $pepinfo->start; #store the peptide starts
	}
	
}

close($EXON_FILE);
my %reported_coords = ();
#Check all db_variant entries
for my $chr ( sort keys %vep_input_info ) {
    for my $coord ( sort {$a<=>$b} keys %{$vep_input_info{$chr}} ) {
        next if exists $reported_coords{"$chr:$coord"}; 
        
        #Check if the coordinate overlaps an exon
        if (exists $coord_to_pep{$chr} && exists $coord_to_pep{$chr}{$coord}) {
                if (exists $pep_to_coord{$coord_to_pep{$chr}{$coord}} && @{$pep_to_coord{$coord_to_pep{$chr}{$coord}}} > 1) {
                        #Here we combine the entries from the same codon
                        my $full_ref = my $full_var;
                        my $max_coord =  0;
                        my $min_coord = 100000000000;

                        my @fields = split("\t",$vep_input_info{$chr}{$coord});
						my $last_coord = $coord;

                        for my $vep_input ( @{$pep_to_coord{$coord_to_pep{$chr}{$coord}}} ) {
                                my ($vep_chr,$vep_coord) = split(':',$vep_input);
                                if ($vep_coord < $min_coord) {
                                        $min_coord = $vep_coord;
                                }
                                if ($vep_coord > $max_coord) {
                                        $max_coord = $vep_coord;
                                }
			
                                my @vep_fields = split("\t",$vep_input_info{$vep_chr}{$vep_coord});
                                my ($vep_ref,$vep_var) = split('/',$vep_fields[3]);
            					if (abs($last_coord-$vep_coord) == 2) {
            						#If it's base 1 and 3 of a peptide
            						my $middle_coord = $vep_coord - 1;
            						my $middle_base = $slice_adaptor->fetch_by_region('chromosome',$chr, $middle_coord, $middle_coord)->seq();
            						$full_ref .= $middle_base;
            						$full_var .= $middle_base;
            					}
                                $full_ref .= $vep_ref;
                                $full_var .= $vep_var;
                                $reported_coords{"$vep_chr:$vep_coord"}++;
                                $last_coord = $vep_coord;
                        }

                        my $full_base = $full_ref . '/' . $full_var;

                        print $VEP_INPUT_FILE join("\t",$chr,$min_coord,$max_coord,$full_base,$fields[4]) . "\n";
                        
                        next;
                }
        }
        #Normal case here with nothing to combine
        print $VEP_INPUT_FILE $vep_input_info{$chr}{$coord} . "\n";
        }
}

close($VEP_INPUT_FILE);


# Run annovar and parse result

my $vep 
    = modules::VEP->new(-input_file      => $vep_input_file,
				-executable_path => $vep_executable,
				-db_dir          => $vep_db_dir,
				-working_dir     => $working_dir,
				-source_type		  => $source_type,
				-exon => 1,
				-output_file => $vep_output_file
				);

if (!$vep->run) {
	modules::Exception->throw("ERROR: VEP failed to run");
}

my $results = $vep->parse_result;

# Get a list of chromosomes included in the results

my %chrs;
foreach my $result (@$results){
    $chrs{$result->[0]}++;
}

# Work through results by chromosome. Write exon_ns matches to database and intermediate files.
open(my $EXONNS, ">$exonns_file") 
    or modules::Exception->throw("Unable to open file for writing");


my %mutant_snv_data;
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');


my @exonns_inserts = ();

foreach my $chr (sort keys %chrs){

    my %chr_snvs;


	my $file_name = $sample_name . '_' . $runid . '.db_variants.snv.'.$chr.'.xml';
	$var_xml->load_xml_keys(-file_name=>$file_name);
	my $snvs = $var_xml->search_xml_keys(-chr=>$chr,-var_type=>'snv');
	
	for my $snv_obj (@{$snvs}) {
		$chr_snvs{$snv_obj->{start_coord}}{$snv_obj->{var_bases}} = $snv_obj->{chr};
	}	    				  
	

    foreach my $result (@$results){
		next unless $result->[0] eq $chr;
	
		my $no_middle = 0;
	
		my ($snv_chr,
		    $snv_start,
		    $snv_end,
		    $var_type,
			$se_type,
		    $var_base_str,
		    $snv_gene,
		    $snv_exon,
		    $snv_ref_aa,
		    $snv_var_aa,
		    $poly_predict,
		    $poly_score,
		    $sift_predict,
		    $sift_score,
		    $cadd_phred) = @$result;
	
		my %coords = ();
		my @var_bases = split("",$var_base_str);
		if (@var_bases == 3 && $snv_end-$snv_start > 2) {
			#Determine whether the middle base borders the start or the end coord when spanning introns
			my $start_plus = $snv_start+1;
			my $end_minus = $snv_end-1;
			
			if (exists $chr_snvs{$start_plus}{$var_bases[1]}) {
				$coords{$start_plus}++;
			} elsif (exists $chr_snvs{$end_minus}{$var_bases[1]}) {
				$coords{$end_minus}++;
			} else {
        		modules::Exception->throw("ERROR: Can't find snv for three base codon mutations spanning introns");
			}
		} elsif (@var_bases == 3) {
			my $middle_coord = $snv_start + 1;
			#Check it's not a middle base we padded
			my $middle_base = $slice_adaptor->fetch_by_region('chromosome',$chr, $middle_coord, $middle_coord)->seq();
			if ($middle_base ne $var_bases[1]) {
				$coords{$middle_coord}++;
			} else {
				$no_middle = 1;
			}
		}
        
        my $info_field = join('^^^',$snv_gene, $snv_exon, $poly_predict, $poly_score, $sift_predict, $sift_score);
		
		my $combined = $snv_start == $snv_end?0:1;
		$coords{$snv_start}++;
		$coords{$snv_end}++;
		if ($se_type =~ /missense_variant/ || $se_type =~ /stop_gained/){
		
			
			my $var_base_count = 0;
			
			for my $snv_count (sort {$a<=>$b} keys %coords) {
			    my $attribute;
			    my $var_base = $var_bases[$var_base_count];
			    my $snv_id = $chr_snvs{$snv_count}{$var_base};
			    
			    if ($snv_id !~ /\w/) {
					print Dumper \@var_bases;
				    modules::Exception->throw("ERROR: Can't retrieve snv for $snv_chr $snv_count $var_base $var_base_count\n");
			    }
			
				if ($combined) {
					$attribute = 'gene=' . $snv_gene 
						                    . ';exon=' . $snv_exon 
						                    . ';aa_change=' . $snv_ref_aa . '->' . $snv_var_aa
						                    . ';poly_pred='. $poly_predict
						                    . ';poly_score='. $poly_score
						                    . ';sift_pred='. $sift_predict
						                    . ';sift_score='.$sift_score
						                    . ';cadd_phred='.$cadd_phred
						                    . ';combined='.$snv_start.'-'.$snv_end;
					print $EXONNS 
								join("\t", 
				     			$snv_chr, 
				     			$snv_count, 
				     			$snv_count,
				     			'') 
								. join('^^^', 
				       					$info_field,
				       					$snv_ref_aa . '->' . $snv_var_aa, 'combined') . "\n";
				} else {
					$attribute = 'gene=' . $snv_gene 
						                    . ';exon=' . $snv_exon 
						                    . ';aa_change=' . $snv_ref_aa . '->' . $snv_var_aa
						                    . ';poly_pred='. $poly_predict
						                    . ';poly_score='. $poly_score
						                    . ';sift_pred='. $sift_predict
						                    . ';sift_score='.$sift_score
						                    . ';cadd_phred='.$cadd_phred;
						                    
						                    
					print $EXONNS 
								join("\t", 
				     			$snv_chr, 
				     			$snv_count, 
				     			$snv_count,
				     			'') 
								. join('^^^', 
				       					$info_field,
				       					$snv_ref_aa . '->' . $snv_var_aa) . "\n";
				}
			    
		 		my $mutant_snv_key =  $snv_chr.":".$snv_count;
				$mutant_snv_data{$mutant_snv_key}{$snv_count}{$var_base}{snv_filter_exon_ns_string} = $attribute;
				$mutant_snv_data{$mutant_snv_key}{$snv_count}{$var_base}{pass} = 1;
			 	
			 	if ($no_middle) {
					$var_base_count += 2; #Skip the middle base in this rare case		 		
			 	} else {
			    	$var_base_count++;
			 	}
			}
		}
    }
}

my $file_name = $sample_name . '_' . $runid . '.filter_exon_ns.snv.xml';
$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_snv_data, -chr=>'all');
my $full_xml_file = $run_dir . '/conf/'.$file_name;
$var_xml->split_xml(-file_name=>$full_xml_file);

close($EXONNS);


