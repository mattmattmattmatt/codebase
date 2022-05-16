#! /usr/bin/perl -w

use strict;
use modules::Adaptors::VariantDB;
use modules::Exception;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);
use modules::Adaptors::Sample;
use modules::Adaptors::SNV;
use modules::Adaptors::Variant;

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "ref=s",
	   "gene_file=s",
	   "sample_name_regex=s",
	   "gene=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || (!$OPT{gene_file} && !$OPT{gene}));

=pod

=head1 SYNOPSIS

check_for_gene.pl -ref default=GRCh37 -gene_list file_with_gene_names -gene search_a_single_gene -sample_name_regex only_search_regex_samples

Required flags: (-gene || -gene_file)

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

check_for_gene.pl -> find variants for a particular gene(s)

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

=cut

my $svndir = $ENV{'SVNDIR'};
if (!-d $svndir) {
	modules::Exception->throw("ERROR:You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory");
}

if (defined $OPT{gene} && defined $OPT{gene_list}) {
	modules::Exception->throw("ERROR: Only define -gene or -gene_list");
}

my %input_genes = ();

if ($OPT{gene}) {
	$input_genes{uc($OPT{gene})}++;
} else {
	open(GENE,"$OPT{gene_file}") || modules::Exception->throw("Can't open file $OPT{gene_file}\n");
	while (<GENE>) {
		chomp;
		my ($gene_name) = $_ =~ /(\S+)/;
		$input_genes{uc($gene_name)}++;
	}
}




my $ref = defined $OPT{ref}?$OPT{ref}:'GRCh37';
my $conf_dir;

if ($ref =~ /GRCh37/) {
	$conf_dir = $svndir."/conf/human/".$ref;
} else {
	$conf_dir = $svndir."/conf/mouse/".$ref;
}

if (!-d $conf_dir) {
	modules::Exception->throw("ERROR: dir $conf_dir doesn't exist");
}

my $gene_dir = &GetLatest('gene');

my $gene_map_file = $gene_dir . '/' . $ref . '.gene.all';

if (!-e $gene_map_file) {
	modules::Exception->throw("ERROR: Gene map file $gene_map_file doesn't exist");
}

my $exon_dir = &GetLatest('exon');

my $exon_map_file = $exon_dir . '/' . $ref . '.exon.overlap.all';

if (!-e $exon_map_file) {
	modules::Exception->throw("ERROR: Exon map file $exon_map_file doesn't exist");
}

my %gene_info = ();

open(GENEMAP,"$gene_map_file") || modules::Exception->throw("Can't open file $gene_map_file\n");

while (<GENEMAP>) {
	chomp;
	my @fields = split("\t",$_);
	if (exists $input_genes{$fields[2]}) {
		my ($ens) = $fields[0] =~ /(ENS\S+)/;
		$gene_info{$ens}{NAME} = $fields[2];
	}
}

close GENEMAP;

open(EXONMAP,"$exon_map_file") || modules::Exception->throw("Can't open file $exon_map_file\n");

while (<EXONMAP>) {
	chomp;
	my @fields = split(" ",$_);
	
	my @exons = split(',',$fields[3]);
	
	my ($ens) = $exons[0] =~ /(ENS\S+)_/;
	if (exists $gene_info{$ens}) {
		my $start = $fields[1]-10;
		my $end = $fields[2]+10;
		$gene_info{$ens}{chr} = $fields[0];
		push @{$gene_info{$ens}{ranges}}, "$start-$end";
	}
}

close EXONMAP;

my %gene_output = ();

my @samples_objs = modules::Adaptors::Sample->search_all();
my %sample_type_types = ();
for my $sample_obj ( @samples_objs ) {
	my $sample_id = $sample_obj->id;
    my $sample_name = $sample_obj->sample_name;
    my $sample_type = $sample_obj->sample_type;
    
    if ($OPT{sample_name_regex}) {
    	next unless $sample_name =~ /$OPT{sample_name_regex}/;
    }
    
    if ($ref =~ /GRCh/) {
    	next if $sample_type =~ /G1/ || $sample_type =~ /ENU/ || $sample_type =~ /mouse/;
    } else {
    	next unless $sample_type =~ /G1/ || $sample_type =~ /ENU/ || $sample_type =~ /mouse/;
    }
    
    print "Search $sample_name\n";
    
	my @run_objs = modules::Adaptors::Run->search(sample_id=>$sample_obj->id);
	if (@run_objs > 1) {
		modules::Exception->throw("ERROR: >1 run for sample $sample_name");
	} elsif (!@run_objs) {
		modules::Exception->warning("No runs for $sample_name");
	} else {
		if ($run_objs[0]->production()) {
			for my $ens (keys %gene_info) {
				my $gene = $gene_info{$ens}{NAME};
				my $chr = $gene_info{$ens}{chr};
				
				for my $range (@{$gene_info{$ens}{ranges}}) {
					my ($start,$end) = $range =~ /(\d+)-(\d+)/;
					my @snv_objs = modules::Adaptors::SNV->search_region($run_objs[0]->id(), 
															     $chr,
															     $start,
															     $end);
					if (@snv_objs) {
						my $snv_coord = $snv_objs[0]->coord;
						my $snv_str = $snv_objs[0]->ref_base.'->'.$snv_objs[0]->var_base;
						
						push @{$gene_output{$gene}{snvs}{$chr}{$snv_coord}{$snv_str}}, $sample_name;		
					}
					
					my @indel_objs = modules::Adaptors::Variant->search_region($run_objs[0]->id(), 
															     $chr,
															     $start,
															     $end);
					if (@indel_objs) {
						my $indel_type = $indel_objs[0]->var_type;
						my $indel_start = $indel_objs[0]->start_coord;
						my $indel_end = $indel_objs[0]->end_coord;
						my $bases = $indel_objs[0]->affected_bases;
						push @{$gene_output{$gene}{$indel_type}{$chr}{$indel_start}{$indel_end}{$bases}}, $sample_name;		
					}
															     
				}
				
				
			}
		} else {
			modules::Exception->warning("No production runs for $sample_name");
		}
	}
	
	
	
}
open(GENE,">gene.tsv") || modules::Exception->throw("Can't open file to write gene.tsv\n");

if (keys %gene_output) {
	print join("\t",
					'VAR_TYPE','COUNT', 'GENE', 'CHR', 'COORD', 'BASES', 'SAMPLES'
				) . "\n\n";	
	print GENE join("\t",
					'VAR_TYPE','COUNT', 'GENE', 'CHR', 'COORD', 'BASES', 'SAMPLES'
				) . "\n\n";					
}

for my $gene (sort keys %gene_output) {
	
	if (exists $gene_output{$gene}{snvs}) {
		for my $snv_chr (sort keys %{$gene_output{$gene}{snvs}}) {
			for my $coord (sort {$a<=>$b} keys %{$gene_output{$gene}{snvs}{$snv_chr}}) {
				for my $change (keys %{$gene_output{$gene}{snvs}{$snv_chr}{$coord}}) {
					my $var_count = @{$gene_output{$gene}{snvs}{$snv_chr}{$coord}{$change}};
					my $samples = join(",",@{$gene_output{$gene}{snvs}{$snv_chr}{$coord}{$change}});
					print join("\t",
								'SNV',$var_count, $gene, $snv_chr, $coord, $change, $samples
								) . "\n";	
					print GENE join("\t",
								'SNV',$var_count, $gene, $snv_chr, $coord, $change, $samples
								) . "\n";	
				}
			}
		}
	}
	
	if (exists $gene_output{$gene}{DEL}) {
		for my $del_chr (sort keys %{$gene_output{$gene}{DEL}}) {
			for my $start_coord (sort {$a<=>$b} keys %{$gene_output{$gene}{DEL}{$del_chr}}) {
				for my $end_coord (sort {$a<=>$b} keys %{$gene_output{$gene}{DEL}{$del_chr}{$start_coord}}) {
					for my $bases (keys %{$gene_output{$gene}{DEL}{$del_chr}{$start_coord}{$end_coord}}) {
						my $var_count = @{$gene_output{$gene}{DEL}{$del_chr}{$start_coord}{$end_coord}{$bases}};
						my $samples = join(",",@{$gene_output{$gene}{DEL}{$del_chr}{$start_coord}{$end_coord}{$bases}});
						print join("\t",
									'DEL',$var_count, $gene, $del_chr, $start_coord, $bases, $samples
									) . "\n";	
						print GENE join("\t",
									'DEL',$var_count, $gene, $del_chr, $start_coord, $bases, $samples
									) . "\n";
					}
				}
			}
		}

	}
	
	if (exists $gene_output{$gene}{INS}) {
		for my $ins_chr (sort keys %{$gene_output{$gene}{INS}}) {
			for my $start_coord (sort {$a<=>$b} keys %{$gene_output{$gene}{DEL}{$ins_chr}}) {
				for my $end_coord (sort {$a<=>$b} keys %{$gene_output{$gene}{INS}{$ins_chr}{$start_coord}}) {
					for my $bases (keys %{$gene_output{$gene}{INS}{$ins_chr}{$start_coord}{$end_coord}}) {
						my $var_count = @{$gene_output{$gene}{INS}{$ins_chr}{$start_coord}{$end_coord}{$bases}};
						my $samples = join(",",@{$gene_output{$gene}{INS}{$ins_chr}{$start_coord}{$end_coord}{$bases}});
						print join("\t",
									'INS',$var_count, $gene, $ins_chr, $start_coord, $bases, $samples
									) . "\n";	
						print GENE join("\t",
									'INS',$var_count, $gene, $ins_chr, $start_coord, $bases, $samples
									) . "\n";	 			
						
					}
				}
			}
		}	
	}
	
	
}



#print Dumper \%gene_output;




#Get the latest directory based on formats DDMMYY
sub GetLatest {
	my ($name) = @_;
	opendir(DIR,"$conf_dir/$name/") || modules::Exception->throw("ERROR: Cannot open directory $conf_dir/$name/");
	my @files = grep {/^\d/} readdir DIR;
	closedir DIR;
	my ($dir_tmp) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @files);
	return "$conf_dir/$name/$dir_tmp";
}
