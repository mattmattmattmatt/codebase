#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use modules::Utils;
use Bio::EnsEMBL::Registry;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   		"help|h",
	   		"man|m",
	   		"gene_file=s",
			"uniprot_file=s",
			"mapped_out=s",
			"unmapped_out=s",
	   		"debug=s",
	   		"organism=s"
	   		);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{gene_file} || !$OPT{uniprot_file});



=pod

=head1 SYNOPSIS

create_polyphen.pl -organism organism(default=mouse) -uniprot_out uniprot_output_mapping_file(default=uniprot.out) -ensembl_out ensembl_output_mapping_file(default=ensembl.out) -uniprot_file uniprot_to_ens_mapping_file -gene_file gene_mapping_file(see README.conf for format) -debug [options]

Required flags: -gene_file -uniprot_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

create_polyphen.pl -> Script to generate a polyphen aa mapping chart

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

/home/matt/work/pipeline/conf/parsers/create_polyphen.pl -gene_file ~/work/pipeline/conf/mouse/NCBIM37/gene/160611/NCBIM37.BiomartMapping.txt -uniprot_file ~/software/ngs/suites/polyphen/uniprot/mouse.seq -uniprot_out NCBIM37.uniprot_mapping.txt -ensembl_out NCBIM37.uniprot_nomapping.txt

=cut

my $gene_file = $OPT{gene_file};
if ( ! -e $gene_file ) {
	modules::Exception->throw("File $gene_file doesn't exist");
}

my $uniprot_file = $OPT{uniprot_file};
if ( ! -e $uniprot_file ) {
	modules::Exception->throw("File $uniprot_file doesn't exist");
}

my $organism = defined $OPT{organism}?$OPT{organism}:'mouse';

my $unmapped_out = defined $OPT{unmapped_out}?$OPT{unmapped_out}:$organism."_unmapped";
my $mapped_out = defined $OPT{mapped_out}?$OPT{mapped_out}:$organism . "_mapped";

open(NOMAP,">$unmapped_out") || modules::Exception->throw("Can't open file to write $unmapped_out\n");
open(MAP,">$mapped_out") || modules::Exception->throw("Can't open file to write $mapped_out\n");


my %gene_mapper = ();
my $registry = 'Bio::EnsEMBL::Registry';

print STDERR "Loading ensembl registry (this could be slow)\n";

my $reg_file = $ENV{SVNDIR} . '/conf/ensembl_registry2.conf';

$registry->load_all($reg_file);
my $slice_adaptor = $registry->get_adaptor($organism, 'Core', 'Slice');

my $gene_adaptor = $registry->get_adaptor($organism, 'Core', 'Gene');

my %mapped_ensembl = ();
my %total_ensembl = ();

#Map the two gene sets
my %uniprot_to_ensembl = ();

open(GENE,$gene_file) || modules::Exception->throw("Can't open file $gene_file");
#Get the biomart to ccds mapping info
#while (<GENE>) {
#	chomp;
#	my @fields = split("\t");
#	
#	my ($ensembl) = $fields[0] =~ /(ENS.*)$/;
#	$total_ensembl{$ensembl} = 0;
#	next if /NO_UNIPROT/;
#	
#	my @uniprots = split(",",$fields[2]);
#
#
#	for my $uniprot ( @uniprots ) {
#		$uniprot =~ s/\.\d+//;
#		$uniprot =~ s/-\d+//;
#		$uniprot_to_ensembl{$uniprot}{$ensembl}++; 
#	}	
#}

my %transcript_to_gene = ();
my $tmp_count = 0;

#Get the canonical transcript names and map to uniprot ids
while (<GENE>) {
	my ($ens_gene,$uniprot) = split();
	if ($OPT{debug}) {
		next unless $OPT{debug} eq $ens_gene;
	}
	my $gene = $gene_adaptor->fetch_by_stable_id($ens_gene);
   	my $transcript = $gene->canonical_transcript();    
	my $name = $transcript->stable_id();
	$transcript_to_gene{$name} = $ens_gene;
	$uniprot_to_ensembl{$uniprot}{$name}++;
	$total_ensembl{$name} = 1;
	$tmp_count++;
	if ($tmp_count % 1000 == 0) {
		print "Getting canonical transcripts $tmp_count ($name->$uniprot)`\n";
	}
	if ($OPT{debug}) {
		print "Gene $ens_gene Canonical trans $name\n";
	}
}



#Get the protein sequence from the database polyphen uses
open(UNIPROT,"$uniprot_file") || die "Can't open file $uniprot_file\n";
my $sp_accession;

#Map the entries you can to uniprot
my $uniprot_count = 0;
my $mapped_count = 0;
while (<UNIPROT>) {
	chomp;
	if (/>(\S+)/) {
		my @fields = split('\|',$1);
		$sp_accession = $fields[1];
		$uniprot_count++;
	} else {
		if (exists $uniprot_to_ensembl{$sp_accession}) {
			if ($OPT{debug}) {
				for my $ens_trans (keys %{$uniprot_to_ensembl{$sp_accession}}) {
					next unless $OPT{debug} eq $transcript_to_gene{$ens_trans};
				}
			}
			my ($match_ensembl,$ens_coord,$ens_strand,$ens_chr,$ens_name) = &Check_Ensembl($sp_accession,$_);
	    		    			    	
			if ($match_ensembl) {
    			&Report($ens_coord,$ens_strand,$ens_chr,$_,$ens_name,$sp_accession);
    			$mapped_count++;
    			$mapped_ensembl{$ens_name} = 1;	
			   	print "PASS_ENS: SP:$sp_accession ENS:$ens_name MAPPED:$mapped_count TOTAL:$uniprot_count\n";
			} 
		}
	}
	#last if $uniprot_count == 100;
}

my $total_ensembl = keys %total_ensembl;
my $total_mapped = keys %mapped_ensembl;
my $unmapped_total = 0;

for my $ensembl (sort keys %total_ensembl) {
	next if exists $mapped_ensembl{$ensembl}; #Skip mapped entries
	$unmapped_total++;
	my ($chr,$coord,$strand,$aa) = Get_Ensembl_Info($transcript_to_gene{$ensembl});
	print "NO_MAP: SP:NO_SP ENS:$ensembl UNMAPPED:$unmapped_total\n";
	&Report($coord,$strand,$chr,$aa,$ensembl .'_NO_MAP','NO_SP');
}

print "TOTAL: $total_ensembl MAPPED: $total_mapped\n";

my $total = keys %gene_mapper;
print "ANALYSIS $total entries...\n\n";
my $count = 0;


#Report the matching entries; only report perfect matches
sub Report {
	my ($coord_str,$strand,$chr,$aa_uniprot,$mapped_entry,$sp) = @_;
	my @aa = split("",$aa_uniprot);	
	my $aa_length = @aa;
	my $aa_number = 1;
	my $aa_index_count = 0;
	my @aa_uniprot_coord = ();

	my %numbers_in_set = &Get_Numbers($coord_str);
	my @sorted_numbers = ();
	
	#Product strand affects 
	if ($strand eq '+') {
		@sorted_numbers  = sort {$a<=>$b} keys %numbers_in_set;
	} else {
		@sorted_numbers = reverse sort {$a<=>$b} keys %numbers_in_set;
	}

	#my $numbers = @sorted_numbers;
	#print "$coord_str $numbers\n";
	#print Dumper \@sorted_numbers;

	my $local_count = 0;
	for my $number ( @sorted_numbers ) {
	    if ($local_count % 3 == 0 && $local_count != 0) {
	    	#Always report things on the positive strand b/c snps are on the positive strand
	    	my @sorted_coord = sort {$a<=>$b} @aa_uniprot_coord;
			my $bp_seq;
			#First get the sequence
	    	for my $coord (@sorted_coord) {
	    		$bp_seq .= $slice_adaptor->fetch_by_region('chromosome',$chr, $coord, $coord)->seq();
	    	}
	    	my $position = 0;
			for my $coord (@sorted_coord) {
				my @bases = split("",$bp_seq);
				my $refbase = $bases[$position];
				if ($mapped_entry =~ /NO_MAP/) {
					#Here we weren't able to map the SP entry so we use the ccds aa sequence
					
					#Don't print huge proteins like TTN
					my $aa_print = $aa_uniprot;
					if ($aa_length > 2000) {
						$aa_print = "SEQ_TOO_LONG(".$aa_length." aa)";
					}
						
			
					print NOMAP "$chr $coord $coord $sp:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry:$aa_print\n";
					#print "$chr $coord $coord $sp_accession:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry:$aa_uniprot\n";
				} else {
		    		print MAP "$chr $coord $coord $sp:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry\n";
		    		#print "$chr $coord $coord $sp_accession:$aa_number:$aa[$aa_index_count]:$bp_seq:$refbase:$position:$strand:$mapped_entry\n";
				}
		    	$position++;
			}
	    	@aa_uniprot_coord = ();
	    	$aa_number++;
	    	$aa_index_count++;
	    }
	    push @aa_uniprot_coord, $number;
	    $local_count++;
	}	
	$count++;
}

sub Get_Ensembl_Info {
	my ($ensembl_id) = @_;
	my $match_ensembl = 0;
    my @ensembl_coords = ();
    my $ensembl_coord_str;
    
   
	#This entry crashes the script
	my $gene = $gene_adaptor->fetch_by_stable_id($ensembl_id);
	if (!$gene) {
		print "ERROR: Skip checking ensembl for $ensembl_id\n";
		next;
	}
   	my $transcript = $gene->canonical_transcript();
    
	my $strand = $transcript->strand();
	my $name = $transcript->stable_id();
	my $chr = $transcript->seq_region_name();
	my $bp_ensembl_length = 0;
	my $exon_count = 0;
	
	my $aa_ensembl;
	if ($transcript->translate()) {
		$aa_ensembl = $transcript->translate()->seq();
	} else {
		next;
	}
	
	#If we've found a match get the coordinates
	foreach my $exon (@{$transcript->get_all_Exons}) {
		#my $estring = feature2string($exon);
		my $coding_start = $exon->coding_region_start($transcript);
		my $coding_end = $exon->coding_region_end($transcript);
				
		#Skip if not coding
        if (!defined $coding_start) {
        	next;
        } 
			
		#my $genome_end = $exon->seq_region_end();
		#my $genome_start = $exon->seq_region_start();
       	#my $length = $coding_end - $coding_start + 1;
		#my $genome_length = $genome_end - $genome_start + 1;
		#if ($length != $genome_length) {			
			#Adjust the coordinates for UTRs
#			if ($strand == 1) {
#				if ($exon_count == 0) {
#					$genome_start = $genome_end - $length + 1;
#				} else {
#					$genome_end = $genome_start + $length - 1;
#				}
#			} else {
#				if ($exon_count == 0) {
#					$genome_end = $genome_start + $length - 1;
#				} else {
#					$genome_start = $genome_end - $length + 1;
#				}
#			}
#			$exon_count++;
#		}
		#push @ensembl_coords, $genome_start. '-' . $genome_end;
		push @ensembl_coords, $coding_start. '-' . $coding_end;		
	}
	$ensembl_coord_str = join(",",@ensembl_coords);
	my $return_strand =  $strand == 1?'+':'-';
	return ($chr,$ensembl_coord_str,$return_strand,$aa_ensembl);

			
}

#check whether there is a corresponding ensembl entry when ccds mapping fails
sub Check_Ensembl {
	my ($sp_accession,$aa_uniprot) = @_;
	my $match_ensembl = 0;
    my @ensembl_coords = ();
    my $ensembl_coord_str;
    
    #here we try to use ensembl to resolve the difference in ccds/uniprot transcripts
	for my $ensembl_id ( keys %{$uniprot_to_ensembl{$sp_accession}} ) {
		#This entry crashes the script
		my $gene = $gene_adaptor->fetch_by_stable_id($transcript_to_gene{$ensembl_id});
		if (!$gene) {
			print "ERROR: Skip checking ensembl for $ensembl_id\n";
			next;
		}
	    my $transcript = $gene->canonical_transcript();
	    
	   
		my $strand = $transcript->strand();
		my $name = $transcript->stable_id();
		next if exists $mapped_ensembl{$name}; #Skip transcripts we've already mapped
		my $chr = $transcript->seq_region_name();
		my $bp_ensembl_length = 0;
		my $exon_count = 0;
		
		my $aa_ensembl;
		if ($transcript->translate()) {
			$aa_ensembl = $transcript->translate()->seq();
		} else {
			next;
		}

		#If we've found a match get the coordinates
		if ($aa_ensembl eq $aa_uniprot) {
			foreach my $exon (@{$transcript->get_all_Exons}) {
				#my $estring = feature2string($exon);
				my $coding_start = $exon->coding_region_start($transcript);
				my $coding_end = $exon->coding_region_end($transcript);
				
				#Skip if not coding
            	if (!defined $coding_start) {
            		next;
            	} 
			
				#my $genome_end = $exon->seq_region_end();
				#my $genome_start = $exon->seq_region_start();
#            	my $length = $coding_end - $coding_start + 1;
#				my $genome_length = $genome_end - $genome_start + 1;
#				if ($length != $genome_length) {			
#					#Adjust the coordinates for UTRs
#					if ($strand == 1) {
#						if ($exon_count == 0) {
#							$genome_start = $genome_end - $length + 1;
#						} else {
#							$genome_end = $genome_start + $length - 1;
#						}
#					} else {
#						if ($exon_count == 0) {
#							$genome_end = $genome_start + $length - 1;
#						} else {
#							$genome_start = $genome_end - $length + 1;
#						}
#					}
#					$exon_count++;
#				}
#				push @ensembl_coords, $genome_start. '-' . $genome_end;
				push @ensembl_coords, $coding_start. '-' . $coding_end;		
			}
			$ensembl_coord_str = join(",",@ensembl_coords);
			#print "ENS $ensembl_coord_str\n";
			my $return_strand =  $strand == 1?'+':'-';
			return (1,$ensembl_coord_str,$return_strand,$chr,$name);
		}
		#print "ENS:$aa_ensembl\nUNI:$aa_uniprot\n\n";
	}
	
	return (0,'','','','');
}
	

#Get the numbers in a range of coordinates
sub Get_Numbers {
	my ($coord_str) = @_;
	my @num_groups = split(",",$coord_str);
	my %numbers;
	for my $num_group ( @num_groups ) {
	    my ($start,$end) = $num_group =~ /(\d+)-(\d+)/;
	    while ($start <= $end) {
	    	$numbers{$start}++;
	    	$start++;
	    }
	}
	return %numbers;
	
	
}


