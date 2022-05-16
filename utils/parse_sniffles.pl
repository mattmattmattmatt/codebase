#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);
use Bio::SearchIO;
use modules::Exception;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "sniffles_vcf=s",
	   "output=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sniffles_vcf});



=pod

=head1 SYNOPSIS

parse_sniffles.pl -sniffles_vcf vcf_in -output output_file_name(default=sniffles.txt) [options]

Required flags: -sniffles_vcf

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_sniffles.pl -> Parse sniffles vcf and generate text coord files

=head1 DESCRIPTION

Feb 5th, 2015

sniffles parser

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

>blast_template.pl -fasta_in file.fa -blast_in blast.bls -fasta_out fasta_new.fa

=cut


my $sniffles_vcf = $OPT{sniffles_vcf};

my $output = defined $OPT{output}?$OPT{output}:'sniffles';

if ( !-e $sniffles_vcf ) {
	modules::Exception->throw("File $sniffles_vcf doesn't exist");	
}

my $sniffles_sv = &parse_sniffles_result($sniffles_vcf);
&print_files($sniffles_sv,'./'.$output,'sniffles');


sub print_files {
	my ($results,$file,$sv_caller) = @_;
	open(FILE,">$file.txt") || modules::Exception->throw("Can't open file $file\n");
	open(TRA,">$file.tra") || modules::Exception->throw("Can't open file $file\n");
	for my $sv_type (qw(del dup inv ins bnd)) {
		for my $coord_str1 (sort {my ($a_chr,$a_coord) = $a =~ /(\S+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /(\S+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord }  keys %{$results->{$sv_type}}) {
				for my $coord_str2 ( keys %{$results->{$sv_type}{$coord_str1}}) {
					my ($chr1,$coord1) = split(':',$coord_str1);
					my ($chr2,$coord2) = split(':',$coord_str2);
					my $anno_str = join("^^^",'SVTYPE='.$sv_type,
												'ID='.$results->{$sv_type}{$coord_str1}{$coord_str2}{id},
												'VAR_READS='.$results->{$sv_type}{$coord_str1}{$coord_str2}{var_reads},
											  'REF_READS='.$results->{$sv_type}{$coord_str1}{$coord_str2}{ref_reads},
											  'LEN='.$results->{$sv_type}{$coord_str1}{$coord_str2}{length},
											   'ZYG='.$results->{$sv_type}{$coord_str1}{$coord_str2}{zyg},
											  'COORD='.$coord_str1.'+'.$coord_str2,
											  'sniffles_GFF='.$results->{$sv_type}{$coord_str1}{$coord_str2}{gff}
											 );
					$anno_str .= '^^^ins='.$results->{$sv_type}{$coord_str1}{$coord_str2}{seq} if $sv_type eq 'ins';						 
					
					
					if ($sv_type eq 'bnd') {
						print TRA join ("\t",
									$chr1,
									$coord1,
									$chr2,
									$coord2,
									$anno_str
									) ."\n";
					}  else {
							print FILE join("\t",
									$chr1,
									$coord1,
									$coord2,
									$anno_str
									) ."\n";
					}
					
							
				}
		}
	}
}

sub parse_sniffles_result {
    my ($output_vcf) = @_;
    
    
    open(my $OUTPUT, $output_vcf) || modules::Exception->throw("ERROR: Can't open output file $output_vcf");
	my %sv_data = ();
	
	while (<$OUTPUT>) {
		chomp;
		next if /^#/;
		my ($chr1,$start1,$id,$ref_base,$var_base,undef,undef,$gff_str,undef,$evidence_str) = split("\t");
		my $chr2 = $chr1;
		my ($start2) = $gff_str =~ /END=(\d+);/; 
		my ($length) = $gff_str =~ /SVLEN=\-?(\d+)/;
		
		my ($zyg,$ref_reads,$var_reads) = split(':',$evidence_str);
		
		my ($event_type_uc) = $gff_str =~ /SVTYPE=([A-Z]+)/;
		my $event_type = lc $event_type_uc;
		
		
		if ($event_type eq 'bnd') {
			($chr2,$start2) = $var_base =~ /[\[\]](.+):(\d+)[\[\]]/;
			#Don't record SVs twice
			
			if ($chr2 !~ /\S/) {
				modules::Exception->throw("ERROR with $var_base\n");
			}
			
			if (exists $sv_data{$event_type}{"$chr2:$start2"}) {
				next;
			}
		}
		
		
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{ref_reads} = $ref_reads;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{var_reads} = $var_reads;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{length} = $length;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{gff} = $gff_str;	
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{id} = $id;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{gff} = $gff_str;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{zyg} = $zyg;
				
		
		if ($event_type eq 'ins') {
			if (length($var_base) > 1000) {
				$var_base = '>1000bp';
			} 	
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{seq} = $var_base;
			
		}
	}

	
    return \%sv_data;
}


