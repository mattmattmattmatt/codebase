#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "outdir=s",
	   "gene_file=s",
	   "outfile_base=s",
	   "joint"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

script_name [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

script_name.pl -> One line description

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut


my %data = ();

my %gene_map = ();

my $outdir = defined $OPT{outdir}?$OPT{outdir}:`pwd`;
chomp $outdir;

my $outbase = defined $OPT{outfile_base}?$OPT{outfile_base}:'manta';

my $somatic = defined $OPT{joint}?0:1; #default is somatic


my $gene_file = defined $OPT{gene_file}?$OPT{gene_file}:"$outdir/ensembl_all";

open(ALL,$gene_file) || modules::Exception->throw("Can't open file $gene_file\n");


while (<ALL>) {
	chomp;
    my @fields = split("\t");
    my ($gene) = $fields[0] =~ /g=(\S+)$/;
    shift @fields;
	$gene_map{$gene} = \@fields;    
}



open(SV,"$outdir/${outbase}.all") || modules::Exception->throw("Can't open file $outdir/${outbase}.all");


while (<SV>) {
	chomp;
	my @fields = split("\t");
  	next if /hs37/;
	my @info = ();
	my $id;

	if ($fields[3] =~ /^\d+$/) {
		($id) = $fields[4] =~ /ID=([a-zA-Z]+[0-9:]+)/;
		$data{$id}{coord1} = "$fields[0]:$fields[1]";
		$data{$id}{coord2} = "$fields[2]:$fields[3]";
		@info = split('\^\^\^',$fields[4]);
	} else {
		($id) = $fields[3] =~ /ID=([a-zA-Z:]+[0-9:]+)/;
		$data{$id}{coord1} = "$fields[0]:$fields[1]";
		$data{$id}{coord2} = "$fields[0]:$fields[2]";
		@info = split('\^\^\^',$fields[3]);
	}
	
	$info[1] =~ s/SVTYPE=//;
	$info[3] =~ s/PE=//;
	$info[4] =~ s/SR=//;
	$info[5] =~ s/LEN=//;
	$info[6] =~ s/QUAL=//;
	
	$data{$id}{svtype} = $info[1];
	$data{$id}{pe} = $info[3];
	$data{$id}{sr} = $info[4];
	$data{$id}{length} = $info[5];
	$data{$id}{qual} = $info[6];
	
	if (!$somatic) {
		$info[8] =~ s/SAMPLE_COUNT=//;
		$data{$id}{sample_count} = $info[8];
	}
	
	#shift @info;
	#shift @info;
	#shift @info;
	#shift @info;
	#$data{$id}{rest} = \@info;
}

open(GENE,"$outdir/${outbase}.all.gene") || modules::Exception->throw("Can't open file $outdir/${outdir}.all.gene");

while (<GENE>) {
	chomp;
  	next if /hs37/;

	my @fields = split;
	my @sv_fields = split('\^\^\^',$fields[3]);
	my @genes = split('\^\^\^',$fields[4]);
	my $id;
	if ($fields[3] =~ /^\d+$/) {
		($id) = $fields[4] =~ /ID=([a-zA-Z]+[0-9:]+)/;
	} else {
		($id) = $fields[3] =~ /ID=([a-zA-Z:]+[0-9:]+)/;
	}
  
  	#print Dumper \@genes;	
	my $gene_str;
	if ($genes[0] eq 'NO_GENE') {
			$data{$id}{genes}{"NO_GENE"} = 1;
	} elsif (@genes > 5) {
			$data{$id}{genes}{">5 genes"} = 1;	
	} else {
		for my $gene (@genes) {
			$gene =~ s/-DUP//;
			$data{$id}{genes}{$gene} = \@{$gene_map{$gene}};
		}
	}
	
	
	
}
#print Dumper \%data;

if ($somatic) {
	print join("\t",
				'chr1',
				'coord1',
				'chr2',
				'coord1',
				'ID',
				'somatic_score',
				'svtype',
				'supporting_paired_end_reads',
				'supporting_split_reads',
				'length',
				'gene_info'
				) ."\n\n";
} else {
	print join("\t",
				'chr1',
				'coord1',
				'chr2',
				'coord1',
				'ID',
				'quality_score',
				'sample_zyg',
				'svtype',
				'supporting_paired_end_reads',
				'supporting_split_reads',
				'length',
				'gene_info'
				) ."\n\n";
}

for my $id ( sort keys %data ) {
	my ($chr1,$coord1) = split(':',$data{$id}{coord1});
	my ($chr2,$coord2) = split(':',$data{$id}{coord2});
	
	my $length = 0;
	if ($data{$id}{length} =~ /\d/) {
		$length = $data{$id}{length};
	}
	
	if ($somatic) {
    	print join ("\t",
    			$chr1,
    			$coord1,
    			$chr2,
    			$coord2,
    			$id,
    			$data{$id}{qual},
    			$data{$id}{svtype},
    			$data{$id}{pe},
    			$data{$id}{sr},
    			$length
    			#@{$data{$id}{rest}}
    			);
	} else {
		print join ("\t",
    			$chr1,
    			$coord1,
    			$chr2,
    			$coord2,
    			$id,
    			$data{$id}{qual},
    			$data{$id}{sample_count},
    			$data{$id}{svtype},
    			$data{$id}{pe},
    			$data{$id}{sr},
    			$length
    			#@{$data{$id}{rest}}
    			);
	}
    
    
   my $useful_gene = 0;
    my $five_genes = 0;	
    #print Dumper 	$data{$id}{genes};
    my $gene_count = 0;	
   	for my $gene (keys %{$data{$id}{genes}}) {
   		if ($gene eq 'NO_GENE') {
   			next;
   		} elsif ($gene =~ /genes/) {
        	$five_genes = 1;
        	next;
      	} else {
   			$useful_gene = 1;
   			if ($somatic) {
	   			if ($gene_count > 0) {
	   				print "\n\t\t\t\t\t\t\t\t\t";
	   			}
   			} else {
   				if ($gene_count > 0) {
   					print "\n\t\t\t\t\t\t\t\t\t\t";
   				}
   			}
   			print "\t". join("\t",
   				@{$data{$id}{genes}{$gene}}
   				);
   		}
   		$gene_count++;
    }
    
    if ($five_genes) {
      print "\t>5 Genes in Region\n";
      next;
    }

    if (!$useful_gene) {
    	print "\tNO_GENE";
    }
    
    
   			
   	print "\n";
    
}

#print Dumper \%data;








