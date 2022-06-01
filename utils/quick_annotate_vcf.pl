#! /usr/bin/perl -w

use strict;
use modules::Exception;
use modules::Vcf;
use modules::SystemCall;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Data::Dumper;
use Cwd 'abs_path';
use vars qw(%OPT);

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vcf_in=s",
	   "outdir=s",
	   "ref=s",
	   "no_zyg",
	   "gene_anno_file=s",
	   "gene_coord_file=s",
	   "gnomad_version=s",
	   "no_run",
	   "skip_vep",
	   "sample_file=s",
	   "control_file=s",
	   "somatic",
	   "outfile=s",
	   "max_nocall_count=i",
	   "min_sample_count=i",
       "chr=s",
       "vartrix_summary=s",
       "sc_min_var=i",
       "sc_min_total=i",
       "sc_min_portion=s",
       "missionbio"
   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf_in});

=pod

=head1 SYNOPSIS

quick_annotate_vcf.pl -vcf_in output_dir -outdir outdir(default=cwd) -min_sample_count minimum_number_samples_with_variant -ref ref_genome(default=GRCh38) -no_zyg no_zyg_info -no_run list_commands_and_quit -gene_anno_file gene_annotation_file -gene_coord_file gene_coordinate_file -skip_vep vep_already_run -control_file for_mixed_datasets_where_specific_controls_are_needed -sample_file only_include_vars_in_these_samples -somatic with_sample_file_vars_are_somatic_and_exclusive_to_sample_list  -outfile output_summary_name(in_outdir) -gnomad_version version(default=2.0.1) -max_nocall_count don't_include_variants_with_this_many_nocalls -chr run_for_chr -vatrix_summary vartrix_file(from_vartrix_count.pl) -sc_min_var sc_min_variant_cells -sc_min_total sc_min_cells_with_data -sc_min_portion portion_sc_with_data_that_are_variant

Required flags: -vcf_in 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

quick_annotate_vcf.pl -> annotate a vcf from multiple samples

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $rare_cutoff = 0.02;

if (!-d $svndir) {
	modules::Exception->throw("ERROR: $svndir doesn't exist\n");
}

my $ref = defined $OPT{ref}?$OPT{ref}:"GRCh38";


my $gnomad_version = defined $OPT{gnomad_version}?$OPT{gnomad_version}:"2.0.1";

my $gnomad_base = $svndir.'/conf/human/'.$ref.'/gnomad/'.$gnomad_version.'/'.$ref.'.gnomAD.overlap';

my $vcf = $OPT{vcf_in};
$vcf = abs_path($vcf);

if ( !-e $vcf ) {
	modules::Exception->throw("File $vcf doesn't exist");	
}

#Vartrix flag
my $vartrix_summary = defined $OPT{vartrix_summary}?$OPT{vartrix_summary}:0;

#MissionBio flag
my $mb = defined $OPT{missionbio}?1:0;

my ($vcf_short) = basename($vcf);
(my $vcf_out = $vcf_short) =~ s/.vcf/.txt/;

my $chr_filter = defined $OPT{chr}?$OPT{chr}:"all";

my $outdir = defined $OPT{outdir}?$OPT{outdir}:`pwd`;
chomp $outdir;
$outdir = abs_path($outdir);

if (!-d $outdir) {
        modules::Exception->throw("Dir $outdir doesn't exist");
}

#Either pass sample file with -somatic flag (everything else treated as control) or pass in control AND sample file for more complex subsetting

if ($OPT{somatic} && $OPT{control_file}) {
	modules::Exception->throw("Only use control_file or somatic flag. Somatic flag treats everything not in sample as control");
}

if ($OPT{somatic} && !$OPT{sample_file}) {
	modules::Exception->throw("ERROR: Can't use somatic option without sample_file\n");
}

if ($OPT{control_file} && !$OPT{sample_file}) {
	modules::Exception->throw("ERROR: Can't use control option without sample_file\n");
}


my $parse_vcf = "$svndir/utils/parse_vcf.pl";
my $overlap_bin = "$svndir/utils/overlap_files.pl";
my $vep_wrapper = "$svndir/utils/vep_wrapper.pl";
my $conf_dir = "$svndir/conf/human/".$ref;


my $pipe_config = modules::Pipeline::get_pipe_conf();
my @chrs = split(" ",$pipe_config->read('human_single_gatk','annotation_version','chr'));

my $gene_dir = &GetLatest('gene');

my @var_types = qw(snv indel);

my $gene_coord_file = my $gene_anno_file;

my $zyg = defined $OPT{no_zyg}?0:1;

#First parse the annotations file -> this is joined on ENSEMBL gene name

if ($OPT{gene_anno_file}) {
	$gene_anno_file = $OPT{gene_anno_file};
} else {
	my $file = $ref . '.gene.all';
	$gene_anno_file = $gene_dir . '/' . $file;
}

if ( !-e $gene_anno_file ) {
	modules::Exception->throw("File $gene_anno_file doesn't exist");	
}

my $subset = 0;


#Look at specific samples
my %samples = ();
if ($OPT{sample_file}) {
	open(SAMPLE,"$OPT{sample_file}") || modules::Exception->throw("Can't open file $OPT{sample_file}\n");
	while (<SAMPLE>) {
		chomp;
		my ($sample) = $_ =~ /(\S+)/;
		$samples{$sample} = 1;
	}
}

my $somatic = defined $OPT{somatic}?1:0;

#Look at specific samples
my %controls = ();
if ($OPT{control_file}) {
	open(CONTROL,"$OPT{control_file}") || modules::Exception->throw("Can't open file $OPT{control_file}\n");
	while (<CONTROL>) {
		chomp;
		my ($sample) = $_ =~ /(\S+)/;
		$controls{$sample} = 1;
	}
}


my @anno_headers = ();
my $first_line = 1;
my $anno_count;
my %gene_anno = ();

open(ANNO,"$gene_anno_file") || modules::Exception->throw("Can't open file $gene_anno_file\n");

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
	    $gene_anno{$ens_gene} = \@fields;
    }
}

my @no_anno_line = ();
push @no_anno_line, 'NO_ANNO' for (1..$anno_count); 

#Gene coord file is used for overlap variants -> vep isn't enough as not run on indels

if ($OPT{gene_coord_file}) {
	$gene_coord_file = $OPT{gene_coord_file};
} else {
	my $file = $ref . '.gene.overlap.all';
	$gene_coord_file = $gene_dir . '/'.$file;
}
	
if ( !-e $gene_coord_file ) {
	modules::Exception->throw("File $gene_coord_file doesn't exist");	
}

#Build up command list
my @commands = ();
my $zyg_str = $zyg?' -keep_zyg ':'';

my $vep_in_command ="cat $outdir/$vcf_out.snv". ' | sed -e "s/:/ /g" -e "s/;/ /g" -e "s/->/ /" | awk \'{print $1,$2,$3,$7,$8,"+"}'."' > $outdir/vep.in"; 

push @commands, "$parse_vcf -vcf $vcf $zyg_str -out $outdir/$vcf_out";
push @commands, "grep SNV $outdir/$vcf_out > $outdir/$vcf_out.snv";
push @commands, "grep -v SNV $outdir/$vcf_out > $outdir/$vcf_out.indel";
push @commands, $vep_in_command;
push @commands, "$vep_wrapper -vep_bin $svndir/ext/bin/vep -vep_in $outdir/vep.in > $outdir/$vcf_out.vep.exon" unless $OPT{skip_vep};
push @commands, "$vep_wrapper -vep_bin $svndir/ext/bin/vep -vep_in $outdir/vep.in -all > $outdir/$vcf_out.vep.all" unless $OPT{skip_vep};
push @commands, "$overlap_bin -ref $outdir/$vcf_out -coord $gene_coord_file -just_overlap -all > $outdir/$vcf_out.gene_coord";



#Gnomad needs to be split by chr
for my $chr (@chrs) {
	next if $chr eq 'Y';
	
	if ($chr_filter =~ /[0-9X]/) {
    	next unless $chr_filter eq $chr;
  	}
	
	push @commands, "grep -w ^$chr $outdir/$vcf_out.snv > $outdir/$vcf_out.snv.gnomad.$chr";
	push @commands, "grep -w ^$chr $outdir/$vcf_out.indel > $outdir/$vcf_out.indel.gnomad.$chr";
	push @commands, "$overlap_bin -ref $outdir/$vcf_out.snv.gnomad.$chr -coord $gnomad_base.snv.$chr -just_overlap -all > $outdir/gnomad.snv.$chr";
	push @commands, "$overlap_bin -ref $outdir/$vcf_out.indel.gnomad.$chr -coord $gnomad_base.indel.$chr -just_overlap -all > $outdir/gnomad.indel.$chr";
}

push @commands, "cat $outdir/gnomad.snv.* >> $outdir/$vcf_out.gnomad";
push @commands, "cat $outdir/gnomad.indel.* >> $outdir/$vcf_out.gnomad";

my $sys_call = modules::SystemCall->new();

#Run the commands
for my $command (@commands) {
	print "$command\n";
	`$command` unless $OPT{no_run};
}


#Now parse the files for the final report
my @samples = ();

open(VCF,"$vcf") || modules::Exception->throw("Can't open file\n");
while (<VCF>) {
	chomp;
	next unless /^#CHROM/;
	my @fields = split("\t",$_);
	for my $field (@fields) { 
		next if $field eq '#CHROM';
		next if $field eq 'POS';
		next if $field eq 'ID';
		next if $field eq 'REF';
		next if $field eq 'ALT';
		next if $field eq 'QUAL';
		next if $field eq 'FILTER';
		next if $field eq 'INFO';
		next if $field eq 'FORMAT';
		push @samples, $field;
	}
  last;
}

my $sample_count = @samples;
my $max_nocall_count = defined $OPT{max_nocall_count}?$OPT{max_nocall_count}:$sample_count;
my $min_sample_count = defined $OPT{min_sample_count}?$OPT{min_sample_count}:1;


my %data = ();

open(PARSED,"$outdir/$vcf_out") || modules::Exception->throw("Can't open file\n");

my $line_count = 0;
#vartrix doesn't report variant base so need to record it
my %vartrix_lookup = ();

#multiple allele handling for counting
my %mult_allele = ();
my %total_alleles = (); #Use for generating average score per variant cell

while (<PARSED>) {
    $_ =~ s/^chr//;
	chomp;
    my @fields = split ("\t");
	my ($chr,$start,$end,$data,@genotypes) = split("\t");

  	if ($chr_filter =~ /[0-9X]/) {
    	next unless $chr_filter eq $chr;
  	}

	my ($var_type,$var_base_str,$qual,$var_count) = $data =~ /([A-Z]+);.*:(\S+);Q=(\S+);AC=(\d+)/;
	my $var_base = my $ref_base;
	if ($var_type eq 'SNV') {
		($ref_base,$var_base) = split('->',$var_base_str); 
	} else {
		$var_base = $var_base_str;
		$ref_base = 'N/A';
	} 
	my $key = "$chr:$start:$end:$var_base";
	
	my $allele_count = 1;
	
	#Search for next allele (i.e. 0/2 not 0/1)
	if (exists $mult_allele{"$chr:$start:$end"}) {
		$mult_allele{"$chr:$start:$end"}++;
		$allele_count = $mult_allele{"$chr:$start:$end"};
	} else {
		$mult_allele{"$chr:$start:$end"} = 1;
	}
	
	if ($vartrix_summary) {
		$vartrix_lookup{"$chr:$start"} = $key;
	}
	$data{$key}{var_type} = $var_type;
	$data{$key}{ref} = $ref_base;
	$data{$key}{var} = $var_base; 
	$data{$key}{qual} = $qual;
	$data{$key}{var_count} = $var_count;

	my $zyg;
	my $sample;


	#only if zygosity is included
	for (my $count = 0; $count < @genotypes; $count++) {
		my @geno_fields = split(':',$genotypes[$count]);
		$sample = $samples[$count];
		my ($allele1,$allele2);
		if ($geno_fields[0] =~ /\//) {
			($allele1,$allele2) = split('/',$geno_fields[0]);
		} elsif ($geno_fields[0] =~ /\|/) {
			($allele1,$allele2) = split('\|',$geno_fields[0]);
		} else {
			modules::Exception->throw("ERROR: Can't handle genotype $geno_fields[0]\n");
		}
		if ($allele1 eq '0' && $allele2 eq '0') {
			$zyg = 'ref';
			$data{$key}{ref_count}++;
		} elsif ($geno_fields[0] eq './.' || $geno_fields[0] eq '.|.') {
			$zyg = 'no_call';
			$data{$key}{no_data_count}++;
		} elsif ($allele1 == $allele2) {
			if ($allele_count == $allele1) {
				$zyg = 'hom';
				$data{$key}{hom_count}++;
				push @{$data{$key}{var_samples}},$sample;
			}
		} elsif ($allele1 != $allele2) {
			if ($allele1 == $allele_count || $allele2 == $allele_count) {
				$zyg = 'het';
				$data{$key}{het_count}++;
				push @{$data{$key}{var_samples}},$sample;
			} 
		}  else {
			modules::Exception->throw("ERROR with $genotypes[$count]\n");
		}
		
		$data{$key}{zyg}{$sample} = $zyg;
		
	}
	
	$total_alleles{"$chr:$start:$end"} += $var_count;
	#if (exists $data{$key}{var_count}) {
	#	$total_alleles{"$chr:$start:$end:$var_type"} += $data{$key}{var_count};
	#} else {
	#	$total_alleles{"$chr:$start:$end:$var_type"} = $data{$key}{var_count}
	#}
	$line_count++;

  if ($line_count % 100000 == 0) {
    print "Parsing vcf $chr $start\n";
  }
}

print "Parsed vcf...\n";

open(VEPALL,"$outdir/$vcf_out.vep.all") || modules::Exception->throw("Can't open file $outdir/$vcf_out.vep.all\n");

while (<VEPALL>) {
    $_ =~ s/^chr//;
    next unless /^[0-9XY]+\s/;
    
    chomp;
    my @fields = split("\t");
    my $key = $fields[0].':'.$fields[1].':'.$fields[1] .':'.$fields[4];
    
    
    if ($chr_filter =~ /[0-9X]/) {
        next unless $fields[0] eq $chr_filter;
    }
    if (!exists $data{$key}) {
    	next;
    	#modules::Exception->throw("ERROR: Key $key doesn't exist\n");
    }
    $data{$key}{rs} = $fields[5];
    $data{$key}{gmaf} = $fields[6];
    $data{$key}{domain} = $fields[7];
    $data{$key}{pubmed} = $fields[8];
    $data{$key}{clin} = $fields[9];
    $data{$key}{exon_str} = $fields[10]; 
    $data{$key}{ens_gene} = $fields[11];
    $data{$key}{ens_trans} = $fields[12];
}

#print Dumper \%data;

open(VEPEXON,"$outdir/$vcf_out.vep.exon") || modules::Exception->throw("Can't open file $outdir/$vcf_out.vep.exon\n");

while (<VEPEXON>) {
    $_ =~ s/^chr//;
    next unless /^[0-9XY]+\s/;
    chomp;
    my @fields = split("\t");
    my $key = $fields[0].':'.$fields[1].':'.$fields[1] .':'.$fields[4];
    if ($chr_filter =~ /[0-9X]/) {
      	next unless $chr_filter eq $fields[0];
   	}
	if (!exists $data{$key}) {
    	next;
    	#modules::Exception->throw("ERROR: Key $key doesn't exist\n");
    }
    my ($poly_score) = $fields[9] =~ /([0-9\.]+)/;
    my ($sift_score) = $fields[11] =~ /([0-9\.]+)/;
    
    $data{$key}{aa_change} = $fields[5];
    $data{$key}{ens_gene} = $fields[6];
    $data{$key}{ens_trans} = $fields[7];
    $data{$key}{poly_cat} = $fields[8];
    $data{$key}{poly_score} = $poly_score;
	$data{$key}{sift_cat} = $fields[10];
    $data{$key}{sift_score} = $sift_score;
    $data{$key}{cadd_phred} = $fields[12];
}

print "Parsed VEP...\n";

open(GENE,"$outdir/$vcf_out.gene_coord") || modules::Exception->throw("Can't open file\n");

while (<GENE>) {
    chomp;
    $_ =~ s/^chr//;
    my @fields = split;
    if ($chr_filter =~ /[0-9X]/) {
    	next unless $chr_filter eq $fields[0];
    }

    my @annos = split(';',$fields[3]);
    my $var_base;
    my $ens_gene;
    
    if ($fields[-1] =~ /(ENSG\d+)/) {
    	$ens_gene = $1;
    } elsif ($fields[-1] == 1) {
    	$ens_gene = 'NO_GENE';
    } else {
    	modules::Exception->throw("ERROR: Can't have no gene entry $_\n");
    }
    
    if ($annos[0] eq 'SNV') {
    	($var_base) = $annos[1] =~ /->([ACTG])/;
    } elsif ($annos[0] eq 'DEL') {
    	($var_base) = $annos[1] =~ /(\-[ATGC]+)/;
    } elsif ($annos[0] eq 'INS') {
    	($var_base) = $annos[1] =~ /(\+[ATGC]+)/;
    }
    my $key = $fields[0].':'.$fields[1].':'.$fields[2].':'.$var_base;
    
	if (!exists $data{$key}) {
    	modules::Exception->throw("ERROR: Key $key doesn't exist\n");
    }
    if (!exists $data{$key}{ens_gene}) {
	    $data{$key}{ens_gene} = $ens_gene;
    }
    
}

open(GNOMAD,"$outdir/$vcf_out.gnomad") || modules::Exception->throw("Can't open file\n");

while (<GNOMAD>) {
    chomp;
 	$_ =~ s/^chr//;
    my @fields = split("\t");
    if ($chr_filter =~ /[0-9X]/) {
      next unless $chr_filter eq $fields[0];
    }

    my @annos = split(';',$fields[3]);
    my $var_base;
    if ($annos[0] eq 'SNV') {
    	($var_base) = $annos[1] =~ /->([ACTG])/;
    } elsif ($annos[0] eq 'DEL') {
    	($var_base) = $annos[1] =~ /(\-[ATGC]+)/;
    } elsif ($annos[0] eq 'INS') {
    	($var_base) = $annos[1] =~ /(\+[ATGC]+)/;
    }
    my $key = $fields[0].':'.$fields[1].':'.$fields[2].':'.$var_base;
    if (!exists $data{$key}) {
    	modules::Exception->throw("ERROR: Key $key doesn't exist\n");
    }
    my $match = $fields[-1] eq '1'?'NO_GNOMAD':$fields[-1];
    $data{$key}{gnomad} = $match;
    $data{$key}{gnomad} =~ s/^\d+://;
}
close GNOMAD;

print "Parsed GNOMAD...\n";

if ($vartrix_summary) {
	open(VARTRIX,"$vartrix_summary") || modules::Exception->throw("Can't open file $vartrix_summary\n");
	my $sc_min_var = defined $OPT{sc_min_var}?$OPT{sc_min_var}:0;
 	my $sc_min_total = defined $OPT{sc_min_total}?$OPT{sc_min_total}:0;      
   	my $sc_min_portion = defined $OPT{sc_min_portion}?$OPT{sc_min_portion}:0;

	while (<VARTRIX>) {
		chomp;
		next unless $_ =~ /:/;
		$_ =~ s/"//g;
		$_ =~ s/^chr//;
		my @fields = split("\t");
		my ($chr,$coord) = split(':',$fields[0]);
		if ($chr_filter =~ /[0-9X]/) {
	      next unless $chr_filter eq $chr;
	    }
	    
	    
	    if (exists $vartrix_lookup{"$chr:$coord"}) {
		    my $key = $vartrix_lookup{"$chr:$coord"};
		   
		    $data{$key}{sc_nd} = $fields[1];
		    $data{$key}{sc_ref} = $fields[2];
		    $data{$key}{sc_variants} = $fields[3] + $fields[4] . ' ( '. $fields[5] .'% )';
	    }
	}
	print "Parsed VARTRIX...\n";
}




my $out = defined $OPT{outfile}?"$outdir/".$OPT{outfile}:"$outdir/${vcf_out}.annotated.tsv";

if ($out !~ /tsv$/) {
	$out .= '.tsv';
}

(my $out_short = $out) =~ s/.tsv//;
my $out_priority = $out_short . '_rare_missense_nonsense.tsv';

open(OUT,">$out") || modules::Exception->throw("Can't open file to write\n");
open(PRIORITY,">$out_priority") || modules::Exception->throw("Can't open file to write $out_priority\n");

my @fhs = (*OUT,*PRIORITY);


for my $fh ( @fhs ) {
	print $fh join("\t",
					'chr',
					'start',
					'end',
					'total_sample_quals',
					'average_qual_per_sample',
					'variant_count (het/hom)',
					'ref_count',
					'no_data_count') ."\t";
    
}




if ($vartrix_summary) {
	#Extra colums reported
	for my $fh ( @fhs ) {
		print $fh join("\t",
						'singlecell_nodata',
						'singlecell_ref',
						'singlecell_variants (singlecell_variant_%)'
						) ."\t";
	}
} 

for my $fh ( @fhs ) {	
	print $fh join("\t",			
					'variant_samples',
					'var_type',
					'ref_base',
					'var_base',
					'ens_gene',
					'ens_trans',
					'dbsnp',
					'gmaf',
					'gnomad',
					'aa_change',
					'poly_cat',
					'poly_score',
					'sift_cat',
					'sift_score',
					'cadd_phred',
					'domain',
					'pubmed',
					'clinical',
					'ensembl_link',
					'ensembl_canonical_transcript',
					'gene_name',
					'cosmic',
					'vogelstein_gene',
					'uniprot',
					'ccds',
					'refseq',
					'gene_desc',
					'omim',
					'go_term',
					'mammalian_phenotype',
					@anno_headers
					);
}



if ($zyg) {
	for my $fh ( @fhs ) {
		print $fh "\t";
		print $fh join("\t",@samples);
	}
} 
for my $fh ( @fhs ) {
	print $fh "\n\n";
}

my @keys = sort { my ($a_chr,$a_coord) = $a =~ /([0-9X]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9X]+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord } keys(%data);

for my $key (@keys) {
	my ($chr,$start,$end,$var_base) = split(":",$key);
	if ($chr_filter =~ /[0-9X]/) {
    		next unless $chr_filter eq $chr;
  	}

  	next unless $chr =~ /^[0-9X]/;
	my $aa_change = exists $data{$key}{aa_change}?$data{$key}{aa_change}:'NO_AA_CHANGE';
	my $ens_trans = exists $data{$key}{ens_trans}?$data{$key}{ens_trans}:'NO_ENS_TRANS';
	my $ens_gene = exists $data{$key}{ens_gene}?$data{$key}{ens_gene}:'NO_ENS_GENE';
	my $poly_cat = exists $data{$key}{poly_cat}?$data{$key}{poly_cat}:'NO_POLY_CAT';
	my $poly_score = exists $data{$key}{poly_score}?$data{$key}{poly_score}:'NO_POLY_SCORE';
	my $sift_cat = exists $data{$key}{sift_cat}?$data{$key}{sift_cat}:'NO_SIFT_CAT';
	my $sift_score = exists $data{$key}{sift_score}?$data{$key}{sift_score}:'NO_SIFT_SCORE';
	my $cadd_phred = exists $data{$key}{cadd_phred}?$data{$key}{cadd_phred}:'NO_CADD_SCORE';
	my $gnomad = exists $data{$key}{gnomad}?$data{$key}{gnomad}:'NO_GNOMAD';
	my $domain = !exists $data{$key}{domain} || $data{$key}{domain} eq 'N/A'?'NO_DOMAIN':$data{$key}{domain};
	my $pubmed = !exists $data{$key}{pubmed} || $data{$key}{pubmed} eq 'N/A'?'NO_PUBMED':$data{$key}{pubmed};
	my $clin = !exists $data{$key}{clin} || $data{$key}{clin} eq 'N/A'?'NO_CLIN':$data{$key}{clin};
	my $rs = !exists $data{$key}{rs} || $data{$key}{rs} eq 'N/A'?'NO_DBSNP':$data{$key}{rs};
	my $gmaf = !exists $data{$key}{gmaf} || $data{$key}{gmaf} eq 'N/A'?'NO_GMAF':$data{$key}{gmaf};	
	my $var_samples;
	if (!exists $data{$key}{var_samples}) {
		$var_samples = "Complex overlapping event";
	} else {
		$var_samples = join(",",@{$data{$key}{var_samples}});
	}
	
	my $het_count = exists $data{$key}{het_count}?$data{$key}{het_count}:0;
	my $hom_count = exists $data{$key}{hom_count}?$data{$key}{hom_count}:0;
	my $ref_count = exists $data{$key}{ref_count}?$data{$key}{ref_count}:0;
	my $rogue_count = exists $data{$key}{rogue_count}?$data{$key}{rogue_count}:0;
	my $control_count = exists $data{$key}{control_count}?$data{$key}{control_count}:0;
	my $nd_count = exists $data{$key}{no_data_count}?$data{$key}{no_data_count}:0;
	my $var_count = 0;
	if (exists $data{$key}{var_count}) {
		$var_count =  $data{$key}{var_count};
	}
	my $var_str = $var_count . ' ('.$het_count . '/'. $hom_count .')';
	my $average_score = 'COMPLEX EVENT';

	my $alleles_key = "$chr:$start:$end";

	if (exists  $total_alleles{$alleles_key}) {
		if ($total_alleles{$alleles_key} > 0) {
			$average_score = sprintf("%.2f",$data{$key}{qual} / $total_alleles{$alleles_key});
		} else {
			print "No >0 allele key $key $alleles_key $total_alleles{$alleles_key}\n";
		}	
	} else {
		print "Doesn't exist allele key $key $alleles_key $total_alleles{$alleles_key}\n";
	}
	
	
	if ($nd_count > $max_nocall_count) {
		next;
	}
	
	if (exists $data{$key}{var_count} && $min_sample_count > $data{$key}{var_count}) {
		next;
	}
	
	#Only keep rare nonsense/missense
	my $priority_flag = 1;
	
	if ($aa_change eq 'NO_AA_CHANGE') {
		$priority_flag = 0;
	}
	
	if ($gmaf =~ /\d/) {
		if ($gmaf > $rare_cutoff) {
			$priority_flag = 0;
		}
	}
	
	if ($gnomad =~ /\d/) {
		my @fields = split(':',$gnomad);
		if ($fields[1] > $rare_cutoff) {
			$priority_flag = 0;
		}
	}
	
	
	
	
	my @anno = exists $gene_anno{$ens_gene}?@{$gene_anno{$ens_gene}}:@no_anno_line;
	
	#Here we filter for specific samples
	if (keys %samples) {
		my $report = 0;
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant if at least once required sample
			if (exists $samples{$var_sample}) {
				$report = 1;
			}
		}
		#Means we didn't find any sample containing the variant
		next if $report == 0;
	}
	
	#Now we check if it's somatic
	if ($somatic) {
		my $somatic = 1;
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant in at least once control sample so don't keep report it
			if (!exists $samples{$var_sample}) {
				$somatic = 0;
			}
		}
		#Means a sample not in the list had the variant so it's not somatic
		next unless $somatic == 1;
	} elsif (%controls) {
		my $control = 0;
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant in at least once control sample so don't keep report it
			if (exists $controls{$var_sample}) {
				$control = 1;
			}
		}
		#Skip as found in control
		next unless $control == 0;
		
		#Special handling for sample/control and min_sample_count
		my $sample_include_count = 0;
		
		for my $var_sample ( @{$data{$key}{var_samples}} ) {
			#Here we found the variant if at least once required sample
			if (exists $samples{$var_sample}) {
				$sample_include_count++;
			}
		}
		if ($min_sample_count > $sample_include_count) {
			next;
		}
	}
	
	if ($aa_change =~ /Stop/) {
		$poly_score = 'N/A';
	}
	
	if ($sift_cat eq 'N/A') {
		$sift_score = 'N/A';
	}
	
	if ($vartrix_summary) {
		my $sc_nd = defined $data{$key}{sc_nd}?$data{$key}{sc_nd}:"N/A";
		my $sc_ref = defined $data{$key}{sc_ref}?$data{$key}{sc_ref}:"N/A";
		my $sc_var = defined $data{$key}{sc_variants}?$data{$key}{sc_variants}:"N/A";
		
		#Vartrix specific filters to check; set to -1 to account for 0 cases
		my $sc_min_var = defined $OPT{sc_min_var}?$OPT{sc_min_var}:-1;
 		my $sc_min_total = defined $OPT{sc_min_total}?$OPT{sc_min_total}:-1;      
   		my $sc_min_portion = defined $OPT{sc_min_portion}?$OPT{sc_min_portion}:-1;
		
		
		#Skip filters for unreported SC variants 
		if ($sc_nd ne 'N/A') {
			
			my ($sc_var_num,$sc_var_percent) = $sc_var =~ /(\d+).*\(\s(\d+)/; 
			my $sc_var_portion =  sprintf("%.2f",$sc_var_percent / 100);
			#print "SC num $sc_var_num / SC percent $sc_var_percent / SC portion $sc_var_portion\n";
			my $sc_total = $sc_var_num + $sc_ref;

			if ($sc_total <= $sc_min_total) {
				next;
			}
			
			if ($sc_var_num <= $sc_min_var) {
				next;
			}
			
			if ($sc_var_portion <= $sc_min_portion) {
				next;
			}
			
		}
		 
		 
		 print OUT join("\t",
						$chr,
						$start,
						$end,
						$data{$key}{qual},
						$average_score,
						$var_str,
						$ref_count,
						$nd_count,
						$sc_nd,
						$sc_ref,
						$sc_var,
						$var_samples,
						$data{$key}{var_type},
						$data{$key}{ref},
						$var_base,
						$ens_gene,
						$ens_trans,
						$rs,
						$gmaf,
						$gnomad,
						$aa_change,
						$poly_cat,
						$poly_score,
						$sift_cat,
						$sift_score,
						$cadd_phred,
						$domain,
						$pubmed,
						$clin,
						@anno
				);
				
			if ($priority_flag) {
				print PRIORITY join("\t",
									$chr,
									$start,
									$end,
									$data{$key}{qual},
									$average_score,
									$var_str,
									$ref_count,
									$nd_count,
									$sc_nd,
									$sc_ref,
									$sc_var,
									$var_samples,
									$data{$key}{var_type},
									$data{$key}{ref},
									$var_base,
									$ens_gene,
									$ens_trans,
									$rs,
									$gmaf,
									$gnomad,
									$aa_change,
									$poly_cat,
									$poly_score,
									$sift_cat,
									$sift_score,
									$cadd_phred,
									$domain,
									$pubmed,
									$clin,
									@anno
									);
			}
	} else {
		print OUT join("\t",
						$chr,
						$start,
						$end,
						$data{$key}{qual},
						$average_score,
						$var_str,
						$ref_count,
						$nd_count,
						$var_samples,
						$data{$key}{var_type},
						$data{$key}{ref},
						$var_base,
						$ens_gene,
						$ens_trans,
						$rs,
						$gmaf,
						$gnomad,
						$aa_change,
						$poly_cat,
						$poly_score,
						$sift_cat,
						$sift_score,
						$cadd_phred,
						$domain,
						$pubmed,
						$clin,
						@anno
						);
		if ($priority_flag) {
				print PRIORITY join("\t",
									$chr,
									$start,
									$end,
									$data{$key}{qual},
									$average_score,
									$var_str,
									$ref_count,
									$nd_count,
									$var_samples,
									$data{$key}{var_type},
									$data{$key}{ref},
									$var_base,
									$ens_gene,
									$ens_trans,
									$rs,
									$gmaf,
									$gnomad,
									$aa_change,
									$poly_cat,
									$poly_score,
									$sift_cat,
									$sift_score,
									$cadd_phred,
									$domain,
									$pubmed,
									$clin,
									@anno
									);
		}
	}
	
				
	if ($zyg) {
		print OUT "\t";
		print PRIORITY "\t" if $priority_flag;
		my @sample_zyg;
		for my $sample (@samples) {
			push @sample_zyg, $data{$key}{zyg}{$sample};
		}
		print OUT join("\t",@sample_zyg);
		print PRIORITY join("\t",@sample_zyg) if $priority_flag;
	}
	
	print OUT "\n";			
	print PRIORITY "\n" if $priority_flag;			
}


#Get the latest directory based on formats DDMMYY
sub GetLatest {
	my ($name) = @_;
	opendir(DIR,"$conf_dir/$name/") || modules::Exception->throw("ERROR: Cannot open directory $conf_dir/$name/");
	my @files = grep {/^\d/} readdir DIR;
	closedir DIR;
	my ($dir_tmp) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @files);
	return "$conf_dir/$name/$dir_tmp";
}



