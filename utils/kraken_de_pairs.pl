#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);
use modules::Exception;
use Algorithm::Combinatorics qw(combinations);
use Bio::SeqIO;


GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
		"biom=s",
		"r_dir=s",
		"qiime_map=s",
		"outfile=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{biom});



=pod

=head1 SYNOPSIS

pairs.pl -biom biom_file -qiime_map qiime_map(default='./qiime_map.txt') r_dir R_directory(default='./R/) -outfile outfile(default='./Pairwise_summary.tsv')[options]

Required flags: -biom

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

blast_template.pl -> Parse blast output and append hit info to original fasta header

=head1 DESCRIPTION

Feb 5th, 2015

Blast parsing template

=head1 AUTHOR

./andreas/sample_set/pairs.pl
Matthew Field

=head1 EXAMPLE

>blast_template.pl -fasta_in file.fa -blast_in blast.bls -fasta_out fasta_new.fa

=cut


#Need to set this each time to address specific pairs
my %map = (
			#Group => ['LungNaiveCtrl','LungNaiveT2D','LungInfCtrl','LungInfT2D'],
			Group => ['Control_v2','Control_v5','Infected_v2','Infected_v5']
			);


my $r_dir = defined $OPT{r_dir}?$OPT{r_dir}:'./R';

my $qiime_map = defined $OPT{qiime_map}?$OPT{qiime_map}:'./qiime_map.txt';

if ( !-e $qiime_map ) {
	modules::Exception->throw("Qiime map $qiime_map doesn't exist");
}

my $biom = defined $OPT{biom}?$OPT{biom}:'./kraken.biom';

if ( !-e $biom ) {
	modules::Exception->throw("Biom file $biom doesn't exist");
}

if ( !-e $r_dir ) {
	mkdir($r_dir);
}

my $outfile = defined $OPT{outfile}?$OPT{outfile}:'./Pairwise_summary.tsv';

my $feature_file = $r_dir.'/otu_map.tsv';
my $r_script = $r_dir.'/Pair_commands.R';

open(R_INPUT,">$r_script") || modules::Exception->throw("Can't open file to write\n");
print R_INPUT 'library("phyloseq")' ."\n";
print R_INPUT 'library("DESeq2")' ."\n";
print R_INPUT 'biom <- import_biom("'.$biom.'")' ."\n";
print R_INPUT 'tax_col <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")' ."\n";
print R_INPUT 'map <- import_qiime_sample_data("'.$qiime_map.'")' ."\n";
print R_INPUT 'phylo <- merge_phyloseq(biom,map)' ."\n";
print R_INPUT 'colnames(tax_table(phylo)) <- tax_col' ."\n";



my %files = ();

for my $field ( keys %map ) {
	my $iter = combinations(\@{$map{$field}},2);
	while (my $c = $iter->next) {
		$field =~ s/\d$//;
		my $file = "./R/".$c->[0].'_vs_'.$c->[1].'.tsv';
		$files{$file}++;
		print R_INPUT $field.'_phy_deseq <- phyloseq_to_deseq2(phylo,~'.$field.')' ."\n";
		print R_INPUT $field.'_phy_deseq <- DESeq('.$field.'_phy_deseq, test="Wald", fitType="parametric",sfType="poscounts")' . "\n";
		print R_INPUT 'res <- results('.$field.'_phy_deseq, cooksCutoff = FALSE, contrast = c("'.$field.'","'.$c->[0].'","'.$c->[1].'"))' . "\n";
		print R_INPUT 'full_map <- cbind(as(res, "data.frame"), as(tax_table(phylo)[rownames(res), ], "matrix"))'. "\n";
		print R_INPUT 'write.table(full_map,sep="\t",file="'.$feature_file.'")' ."\n";
		print R_INPUT 'res <- res[which(res$padj<0.1),]' ."\n";
		print R_INPUT 'res <- res[order(res$padj,decreasing = FALSE),]' ."\n";
		print R_INPUT 'write.table(res,file="'.$file.'",sep="\t")'. "\n";
	}
}			
		

system("Rscript $r_script");

open(GG,$feature_file) || modules::Exception->throw("Can't open file $feature_file\n");
my %gg_map = ();

while (<GG>) {
	chomp;
	my @fields = split("\t");
	$fields[0] =~ s/"//g;
	$gg_map{$fields[0]}{tax} = join(";",
								$fields[7],
								$fields[8],
								$fields[9],
								$fields[10],
								$fields[11],
								$fields[12],
								$fields[13]
								);
}


#my $fasta_obj_in = Bio::SeqIO->new(-file => "/drive3/work/microbiome/AIP/final/unmapped_contigs.fa", 
#								-format=>'fasta');
#
#my %fasta_lookup =  ();
#
##Read the fasta entries
#while( my $fasta_seq = $fasta_obj_in->next_seq ) {
#	$fasta_lookup{$fasta_seq->id} = $fasta_seq->seq;
#}

open(SUMMARY,">$outfile") || modules::Exception->throw("Can't open file to write\n");

print SUMMARY
			join("\t",
			"PAIR_COMPARE",
			"SPECIES",
			"OTU_ID",
			"BASEMEAN",
			"FOLD_CHANGE",
			"lfcSE",
			"Stat",
			"Pvalue",
			"PADJ",
				) ."\n";



for my $file ( keys %files ) {
    open(FILE,"$file") || modules::Exception->throw("Can't open file $file\n");
    
    
    my ($pair) = basename($file);
    $pair =~ s/.tsv//;
    
    print SUMMARY "\n$pair\n\n";
    
    while (<FILE>) {
    	next if /basemean/i;
    	chomp;
    	my @fields = split;
    	$fields[0] =~ s/\"//g;
    	my $tax = $gg_map{$fields[0]}{tax};
    	
    	
    	print SUMMARY join("\t",
    						$pair,
    						$tax,
    						$fields[0],
    						$fields[1],
    						$fields[2],
    						$fields[3],
    						$fields[4],
    						$fields[5],
    						$fields[6],
    						) ."\n";
    }
    
}
			
