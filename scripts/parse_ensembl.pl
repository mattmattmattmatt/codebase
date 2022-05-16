#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;
use modules::Utils;
use modules::Pipeline;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Env qw($ENSEMBL_REGISTRY_LOAD);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	 	"ref=s",
	 	"splice_length=i",
	 	"phenotype_file=s",
		"homolog_file=s",
		"immgen_file=s",
		"gnf_file=s",
		"omim_file=s",
		"cosmic_file=s",
		"vogel_file=s",
	   	"debug=s"
	   	);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=poddb

=head1 SYNOPSIS

parse_ensembl.pl -ref reference_genome(default=GRCh37) [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_ensembl.pl -> Script to generate all the gene files required for the pipeline

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./parse_ensemlb.pl

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $pipe_conf = modules::Pipeline::get_pipe_conf();


my %homolog_map = ();
my $ens_link;
my $organism;
my $source_type;
my $ref = defined $OPT{ref}?$OPT{ref}:"GRCh38";

if ($ref =~ /GRCh/) {
	$organism = "human";
	$ens_link = "http://ensembl.org/Homo_sapiens/Gene/Summary?db=core;g="; 
} elsif ($ref =~ /mm/) {
	$organism = "mouse";
	$ens_link = "http://ensembl.org/Mus_musculus/Gene/Summary?db=core;g=";
} else {
	modules::Exception->throw("ERROR: Only set up for mouse and human");
}

my $splice_length = defined $OPT{splice_length}?$OPT{splice_length}:10;

my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_all($ENSEMBL_REGISTRY_LOAD);

# get the adaptors
my $gene_adaptor = $registry->get_adaptor($organism, "Core", "Gene");
my $slice_adaptor = $registry->get_adaptor($organism, 'Core', 'Slice');
my @chr_slices = @{ $slice_adaptor->fetch_all('chromosome') };
my $go_adaptor = $registry->get_adaptor("Multi","Ontology","OntologyTerm");

#print Dumper $go_adaptor;

#Stores the GO term mappings to prevent continually looking up the same db value
my %go_mapping = ();

my $conf_dir = "$svndir/conf/gene_names";
#my $conf_dir = "$svndir/conf/$organism/$ref";

#The output files

#File for filter_exon (per chr files)
my $exon_file = $ref . '.exon.overlap';

#Single exon file
my $exon_all_file = $ref . '.exon.overlap.all';

#File for filter_splicesite (per chr files)
my $splice_file = $ref . '.splice.overlap';

#File for exon coords (used for coverage stats)
my $coord_file = $ref . '.exon.coord';
open(COORD,">$coord_file") || modules::modules->throw("Can't open file to write $coord_file\n");

#gene merging file
my $all_gene_info = $ref. ".gene.all";
open(ALLINFO,">$all_gene_info") || modules::modules->throw("Can't open file to write $all_gene_info\n");

#Store all the data for sorting and writing to files
my %gene_data = ();



#URLs for linking
my $cosmic_url = "http://cancer.sanger.ac.uk/cosmic/gene/overview?ln=";
my $omim_link = 'http://omim.org/entry/';

#We only want to report canonical transcripts so keep track of the ones we see
my %transcripts = ();

#map gene names
my %gene_name_mapper = ();

#Human specific stuff
if ($organism eq 'human') {
	#Use the latest COSMIC
	my $cosmic_file;
	if ($OPT{cosmic_file}) {
		$cosmic_file = $OPT{cosmic_file};
	} else {
		my $cosmic_dir = &GetLatest('cosmic');
		$cosmic_file = $cosmic_dir . '/Cosmic.gene';
	}
	
	#Use the latest Vogelstein
	my $vogel_file;
	if ($OPT{vogel_file}) {
		$vogel_file = $OPT{vogel_file};
	} else {
		my $vogel_dir = &GetLatest('vogelstein');
		$vogel_file = $vogel_dir .'/Vogelstein.gene';
	}

	my $phenotype_file;
	if ($OPT{phenotype_file}) {
		$phenotype_file = $OPT{phenotype_file};
	} else {
		my $phenotype_dir = &GetLatest('phenotype');
		$phenotype_file = $phenotype_dir . '/Phenotype.gene';
	}


	#COSMIC entries are joined by gene name; these non coordinate entries are only reported if there is no coordinate match
	open(COSMIC,"$cosmic_file") || modules::Exception->throw("Can't open file $cosmic_file\n");
	
	while (<COSMIC>) {
		chomp;
		next unless /\w/;
		my ($cosmic_gene) = $_;
		$gene_name_mapper{uc($cosmic_gene)}{cosmic} = $cosmic_url . $cosmic_gene;
	}
	
	
	#VOGELSTEIN entries are joined by gene name
	open(VOGEL,"$vogel_file") || modules::Exception->throw("Can't open file $vogel_file\n");
	
	while (<VOGEL>) {
		chomp;
		next unless /\w/;
		my ($vogel_gene,$info) = split(" ");
		$gene_name_mapper{uc($vogel_gene)}{vogel} = $info;
	}
	
	open(PHEN,"$phenotype_file") || modules::Exception->throw("Can't open file $phenotype_file\n");
	while (<PHEN>) {
		chomp;
		my @fields = split("\t");
		$gene_name_mapper{uc($fields[0])}{phenotype} = $fields[1];
	}
	
} elsif ($organism eq 'mouse') {
	#Mouse specific stuff
	my $omim_file;
	if ($OPT{omim_file}) {
		$omim_file = $OPT{omim_file};
	} else {
		my $omim_dir = &GetLatest('omim');
		$omim_file = $omim_dir . '/Omim.gene';
	}
	
	my $phenotype_file;
	if ($OPT{phenotype_file}) {
		$phenotype_file = $OPT{phenotype_file};
	} else {
		my $phenotype_dir = &GetLatest('phenotype');
		$phenotype_file = $phenotype_dir . '/Phenotype.gene';
	}
	
	my $homolog_file;
	if ($OPT{homolog_file}) {
		$homolog_file = $OPT{homolog_file};
	} else {
		my $homolog_dir = &GetLatest('organism_maps');
		$homolog_file = $homolog_dir . '/Mouse_to_human.gene';
	}
	
	my $immgen_file;
	if ($OPT{immgen_file}) {
		$immgen_file = $OPT{immgen_file};
	} else {
		my $immgen_dir = &GetLatest('expression');
		$immgen_file = $immgen_dir . '/Expression.ImmGen.gene';
	}
	
	my $gnf_file;
	if ($OPT{gnf_file}) {
		$gnf_file = $OPT{gnf_file};
	} else {
		my $gnf_dir = &GetLatest('expression');
		$gnf_file = $gnf_dir . '/Expression.GNF.gene';
	}
	
	
	open(PHEN,"$phenotype_file") || modules::Exception->throw("Can't open file $phenotype_file\n");
	while (<PHEN>) {
		chomp;
		my @fields = split("\t");
		$gene_name_mapper{uc($fields[0])}{phenotype} = $fields[1];
	}
	
	open(HOMOLOG,"$homolog_file") || modules::Exception->throw("Can't open file $homolog_file\n");
	while (<HOMOLOG>) {
		chomp;
		my @fields = split();
		$gene_name_mapper{uc($fields[0])}{homolog} = uc($fields[1]);
	}
	
	#Can't use ensembl for omim details as they're tied to human details
	open(OMIM,"$omim_file") || modules::Exception->throw("Can't open file $omim_file\n");
	while (<OMIM>) {
		chomp;
		my @fields = split('\|',$_);
		my @genes = split(",",$fields[5]);
		for my $gene ( @genes ) {
		    $gene =~ s/\s//g;
			$gene_name_mapper{uc($gene)}{omim} =  $omim_link.$fields[9];
			
			#Map omim homolog
			if (exists $gene_name_mapper{uc($gene)}{homolog}) {
				$gene_name_mapper{uc($gene_name_mapper{uc($gene)}{homolog})}{omim} =  $omim_link.$fields[9];
			}
			
		}
	}
	
	open(IMMGEN,"$immgen_file") || modules::Exception->throw("Can't open file $immgen_file\n");
	while (<IMMGEN>) {
		chomp;
		my @fields = split("\t");
		
		$gene_name_mapper{uc($fields[0])}{immgen} = uc($fields[1]);
	}
	
	open(GNF,"$gnf_file") || modules::Exception->throw("Can't open file $gnf_file\n");
	while (<GNF>) {
		chomp;
		my @fields = split("\t");
		$gene_name_mapper{uc($fields[0])}{gnf} = uc($fields[1]);
	}
} else {
	modules::Exception->throw("ERROR: Only works for mouse and human");
}


# traverse chromosomes
for my $chr_slice (@chr_slices) {
	
	#Get the genes for that chromosome
	my @genes = @{$gene_adaptor->fetch_all_by_Slice($chr_slice)};
	my $chr = $chr_slice->seq_region_name();
	next unless $chr =~ /[0-9XYM]/;
	if ($OPT{chr}) {
		next unless $OPT{chr} eq $chr;
	}
	my $count = 0;
	
	
	print "Processing " . scalar(@genes) . " gene IDs for chr $chr...\n";

	#traverse genes
    for my $gene (@genes) {
	    # let user know count
	    local $| = 1;
	    print "[$count/" . scalar(@genes) . "]\r";
	    $count++;
	
		#ENSG id
		my $gene_id = $gene->stable_id();
		
		if ($OPT{debug}) {
			next unless $OPT{debug} eq $gene_id;
		}

	
	    # get canonical transcript
	    my $canonical_transcript = $gene->canonical_transcript();
	
		my $transcript_id = $canonical_transcript->stable_id();
	
		#only protein coding biotypes	
		my $biotype = $canonical_transcript->biotype();
		#next unless $biotype eq 'protein_coding';
		
		#Stores the sequence info
		#my $sequence;
		my $coding = 0;
	    
		my %uniprot_ids = ();
	
		my @exon_objs = @{$canonical_transcript->get_all_Exons()};
		my $exon_count = 1;
	    foreach my $exon (@exon_objs) {
	    	my $coding_start = $exon->coding_region_start($canonical_transcript);
			next unless defined $coding_start; #Skip non coding exons
			$coding = 1;  #Set the coding flag
			my $coding_end = $exon->coding_region_end($canonical_transcript);
			my $full_ensembl_name = $gene_id . '_exon'.$exon_count;
			push @{$gene_data{$chr}{$coding_start}{$coding_end}},$full_ensembl_name;
			
			
			#get uniprot id from this API call
			my @exon_sf = @{ $exon->get_all_supporting_features() };
			foreach my $sf (@exon_sf) {
				if ($sf->analysis()->display_label() =~ /uniprot/i){
		    		$uniprot_ids{$sf->hseqname()}++;
				}
			}
			
			$exon_count++;
				
	    }
	    
	    #next unless $coding;
		
		
		#gene description
		my $desc = $gene->description();
		if ($desc !~ /\w/) {
			$desc = 'NO_DESC'; #In case there is no description
		}
	
		#external db links
		my @dblinks = @{$gene->get_all_DBLinks()};
	
		if (!@dblinks) {
			modules::Exception->throw("ERROR: Can't get dblinks");
		}
	
		#hugo, go, ccds, and omim values
		my %proper_name = ();
		my %go_terms = ();
		my %ccds = ();
		my %omims = ();
	
		#Get OMIM, CCDS, GO TERMS, and HGNC 
	  	foreach my $dbe (@dblinks){
	  		if ($dbe->dbname eq 'HGNC') {
	  			$proper_name{$dbe->display_id}++;
	  		} elsif ($dbe->dbname eq 'MGI') {
	  			$proper_name{$dbe->display_id}++;
	  		} elsif ($dbe->dbname eq 'CCDS') {
	  			$ccds{$dbe->display_id}++;
	  		} elsif ($dbe->dbname eq 'MIM_GENE') {
	  			my ($omim_id) = $dbe->display_id =~ /\[\*(\d+)/;
	  			$omims{"http://omim.org/entry/$omim_id"}++;
	  		} elsif ($dbe->dbname eq 'GO') {
	  			if (exists $go_mapping{$dbe->display_id}) {
	  				#use the lookup if we've already seen code
	  				$go_terms{$go_mapping{$dbe->display_id}}++;
	  			} else {
					#use the db if never seen code before
					my $tmp = $dbe->display_id;
					
					my $go_obj = $go_adaptor->fetch_by_accession($dbe->display_id);
					if (!defined $go_obj) {
						#print "WARNING NO GO: $gene_id ".$dbe->display_id ."\n";
						next;
					} 
	  				$go_mapping{$dbe->display_id} = $go_obj->name;
	  				$go_terms{$go_obj->name}++;
	  			}
	  		}
	  		
	    	#print "\t".$dbe->dbname."\t".$dbe->display_id."\n";
	  	}
	  	
	  	
	  	#join by comma if multiple values
	  	my $ccds = keys %ccds?join(",",keys %ccds):"NO_CCDS";
		my $go_string = keys %go_terms?join("; ",keys %go_terms):"NO_GO";
	  	my $proper_name;
	  	
	  	#Use the external name if no HGNC name
	  	if (!keys %proper_name) {
	  		$proper_name = $gene->external_name();
	  	} else {
	  		$proper_name = join(",",keys %proper_name);
	  	}
	  	
	  	my $lookup_name = uc($proper_name);
	  	
	  	#Get refseq id from this API call
	  	my @ct_sf = @{ $canonical_transcript->get_all_supporting_features() };


		my %refseq_ids;

		foreach my $ct_sf (@ct_sf) {

	    	if ($ct_sf->analysis()->display_label() =~ /$organism cDNAs/i){
				$refseq_ids{$ct_sf->hseqname()}++
	    	}
		}
	  	
	  	my $refseq = keys %refseq_ids?join(",",keys %refseq_ids):"NO_REFSEQ";
	    my $uniprot_name = keys %uniprot_ids?join(",",keys %uniprot_ids):"NO_UNIPROT";
		my $ens_entry = $ens_link . $gene_id;
		
		#Different reorts for mouse and human
		if ($organism eq 'human') {
			my $omim = keys %omims?join(",",keys %omims):"NO_OMIM";
			my $vogelstein = exists $gene_name_mapper{$lookup_name}{vogel}?$gene_name_mapper{$lookup_name}{vogel}:"NO_VOGEL";
			my $cosmic = exists $gene_name_mapper{$lookup_name}{cosmic}?$gene_name_mapper{$lookup_name}{cosmic}:"NO_COSMIC";
			my $phenotype = exists $gene_name_mapper{$lookup_name}{phenotype}?$gene_name_mapper{$lookup_name}{phenotype}:"NO_PHENOTYPE";
			print ALLINFO join("\t",
	    					$ens_entry,
	    					$transcript_id,
							$proper_name,
							$cosmic,
							$vogelstein,
							$uniprot_name,
							$ccds,
							$refseq,
							$desc,
							$omim,
							$go_string,
							$phenotype
	    					) . "\n";
		
		} elsif ($organism eq 'mouse') {
			my $omim = exists $gene_name_mapper{$lookup_name}{omim}?$gene_name_mapper{$lookup_name}{omim}:"NO_OMIM";
			my $phenotype = exists $gene_name_mapper{$lookup_name}{phenotype}?$gene_name_mapper{$lookup_name}{phenotype}:"NO_PHENOTYPE";
			my $homolog = exists $gene_name_mapper{$lookup_name}{homolog}?$gene_name_mapper{$lookup_name}{homolog}:"NO_HOMOLOG";
			my $immgen = exists $gene_name_mapper{$lookup_name}{immgen}?$gene_name_mapper{$lookup_name}{immgen}:"NO_IMMGEN";
			my $gnf = exists $gene_name_mapper{$lookup_name}{gnf}?$gene_name_mapper{$lookup_name}{gnf}:"NO_GNF";
			print ALLINFO join("\t",
	    					$ens_entry,
							$transcript_id,
							$proper_name,
							$uniprot_name,
							$ccds,
							$refseq,
							$desc,
							$omim,
							$go_string,
							$phenotype,
							$homolog,
							$immgen,
							$gnf
	    					) . "\n";
		}
		
	
		#my $tag = ">$gene_id";
		#print FASTA "$tag\n$sequence\n";
	

	}
}


open(EXONALL,">$exon_all_file") ||modules::Exception->throw("Can't open file to write $exon_all_file\n");

for my $chr (sort keys %gene_data) {
	my %exon_coord;
	my $exon_chr_file = $exon_file . '.' . $chr;
	my $splice_chr_file = $splice_file . '.' . $chr;
	open(EXON,">$exon_chr_file") || modules::Exception->throw("Can't open file to write $exon_chr_file\n");
	open(SPLICE,">$splice_chr_file") || modules::Exception->throw("Can't open file to write $splice_chr_file\n");
	
	for my $start (sort {$a<=>$b} keys %{$gene_data{$chr}}) {		
		
		for my $end (sort {$a<=>$b} keys %{$gene_data{$chr}{$start}}) {
			#Holds the overlap exon info
			
			
			my $exon_name = join(",",@{$gene_data{$chr}{$start}{$end}});
			(my $splice_name1 = $exon_name) =~ s/(exon\d+)/$1_splice1/g;
			(my $splice_name2 = $exon_name) =~ s/(exon\d+)/$1_splice2/g;
				
			
		
			    
		    print EXONALL join(" ",
		    				$chr,
		    				$start,
		    				$end,
		    				$exon_name
		    				). "\n";
		    
		    print EXON join(" ",
		    				$chr,
		    				$start,
		    				$end,
		    				$exon_name
		    				). "\n";
		    
		    my $distance = 10;
		    for my $splice_coord ($start-$splice_length..$start-1) {
			    print SPLICE join(" ",
			    				  $chr,
			    				  $splice_coord,
			    				  $splice_coord,
			    				  $splice_name1.'_'.$distance
			    				) . "\n";
		    	
		    	$distance--;
		    }
		    $distance = 1;
		    for my $splice_coord ($end+1..$end+$splice_length) {
				print SPLICE join(" ",
			    				  $chr,
			    				  $splice_coord,
			    				  $splice_coord,
			    				  $splice_name2.'_'.$distance
			    				) . "\n";				
		    	$distance++;
		    }
		    				
		    				
		    
		    my $start_count = $start;
		    #Avoid duplicate coords; use hash and print at end
		   	while ($start_count <= $end) {
	    		$exon_coord{$start_count}++;
	    		$start_count++;
			}
			    
			
		}
	}
	
	for my $coord (sort {$a<=>$b} keys %exon_coord) {
		print COORD "$chr $coord\n";
	}
	
	close SPLICE;
	close EXON;
	
}


close EXONALL;
close COORD;
close ALLINFO;
#close FASTA;

#Get the latest directory based on formats DDMMYY
sub GetLatest {
	my ($name) = @_;
	opendir(DIR,"$conf_dir/$name/") || modules::Exception->throw("ERROR: Cannot open directory $conf_dir/$name/");
	my @files = grep {/^\d/} readdir DIR;
	closedir DIR;
	my ($dir_tmp) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @files);
	return "$conf_dir/$name/$dir_tmp";
}
