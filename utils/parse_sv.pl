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
	   "sv_vcf=s",
	   "output=s",
	   "sv_caller=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{sv_vcf});



=pod

=head1 SYNOPSIS

parse_sniffles.pl -sniffles_vcf vcf_in -output output_file_name(default=sv_calls) -sv_caller sv_tool(sniffles(default) or pbsv)[options]

Required flags: -sv_vcf

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_sv.pl -> Parse sniffles vcf and generate text coord files

=head1 DESCRIPTION

Feb 5th, 2015

sniffles parser

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

>blast_template.pl -fasta_in file.fa -blast_in blast.bls -fasta_out fasta_new.fa

=cut


my $sv_vcf = $OPT{sv_vcf};

my $output = defined $OPT{output}?$OPT{output}:'sv_calls';

my $sv_caller = defined $OPT{sv_caller}?$OPT{sv_caller}:'sniffles';

if ( !-e $sv_vcf ) {
	modules::Exception->throw("File $sv_vcf doesn't exist");	
}

my $sniffles_sv = &parse_sv_result($sv_vcf);
&print_files($sniffles_sv,'./'.$output,$sv_caller);


sub print_files {
	my ($results,$file,$sv_caller) = @_;
	
	for my $sv_type (qw(del dup inv ins bnd)) {
		
		my $header;
	
		if ($sv_type eq 'bnd') {
			$header = join("\t",
						"CHR1",
						"COORD1",
						"CHR2",
						"COORD2",
						"SVTYPE",
						"VCF_ID",
						"PASS",
						"VAR_READS",
						"REF_READS",
						"VAR_FREQ",
						"LENGTH",
						"ZYG_CALL",
						"COORD_RANGE",
						"GFF"
						) . "\n\n";
		
		} elsif ($sv_type eq 'ins') {
			$header = join("\t",
						"CHR1",
						"COORD",
						"INSERT_BASES",
						"SVTYPE",
						"VCF_ID",
						"PASS",
						"VAR_READS",
						"REF_READS",
						"VAR_FREQ",
						"LENGTH",
						"ZYG_CALL",
						"COORD_RANGE",
						"GFF"
						) . "\n\n";
		} else {
			$header = join("\t",
						"CHR",
						"COORD1",
						"COORD2",
						"SVTYPE",
						"VCF_ID",
						"PASS",
						"VAR_READS",
						"REF_READS",
						"VAR_FREQ",
						"LENGTH",
						"ZYG_CALL",
						"COORD_RANGE",
						"GFF"
						) . "\n\n";
		
		}
		
		open(FILE,">${file}_${sv_type}.tsv") || modules::Exception->throw("Can't open file $file\n");
		print FILE $header;
		for my $coord_str1 (sort {my ($a_chr,$a_coord) = $a =~ /(\S+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /(\S+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord }  keys %{$results->{$sv_type}}) {
				for my $coord_str2 ( keys %{$results->{$sv_type}{$coord_str1}}) {
					my ($chr1,$coord1) = split(':',$coord_str1);
					my ($chr2,$coord2) = split(':',$coord_str2);
					
					my $ins_bases = $sv_type eq 'ins'?$results->{$sv_type}{$coord_str1}{$coord_str2}{seq}:"NA";
					my $var_freq = sprintf("%.2f",$results->{$sv_type}{$coord_str1}{$coord_str2}{var_reads}/($results->{$sv_type}{$coord_str1}{$coord_str2}{var_reads}+$results->{$sv_type}{$coord_str1}{$coord_str2}{ref_reads}));
					
					my $anno_str = join("\t",$sv_type,
											$results->{$sv_type}{$coord_str1}{$coord_str2}{id},
											$results->{$sv_type}{$coord_str1}{$coord_str2}{status},
											$results->{$sv_type}{$coord_str1}{$coord_str2}{var_reads},
											$results->{$sv_type}{$coord_str1}{$coord_str2}{ref_reads},
											$var_freq,
											$results->{$sv_type}{$coord_str1}{$coord_str2}{length},
											$results->{$sv_type}{$coord_str1}{$coord_str2}{zyg},
											$coord_str1.'-'.$coord_str2,
											$results->{$sv_type}{$coord_str1}{$coord_str2}{gff}
											 );
					
					
					if ($sv_type eq 'bnd') {
						print FILE join ("\t",
									$chr1,
									$coord1,
									$chr2,
									$coord2,
									$anno_str
									) ."\n";
					}  elsif ($sv_type eq 'ins') {
						print FILE join ("\t",
									$chr1,
									$coord1,
									$ins_bases,
									$anno_str
									) ."\n";
					} else {
							print FILE join("\t",
									$chr1,
									$coord1,
									$coord2,
									$anno_str
									) ."\n";
					}
					
							
				}
		}
		close FILE;
	}
}

sub parse_sv_result {
    my ($output_vcf) = @_;
    
    
    open(my $OUTPUT, $output_vcf) || modules::Exception->throw("ERROR: Can't open output file $output_vcf");
	my %sv_data = ();
	
	while (<$OUTPUT>) {
		chomp;
		next if /^#/;
		my ($chr1,$start1,$id,$ref_base,$var_base,undef,$pass,$gff_str,undef,$evidence_str) = split("\t");
		my $chr2 = $chr1;
		my ($start2) = $gff_str =~ /END=(\d+);/; 
		my ($length) = $gff_str =~ /SVLEN=\-?(\d+)/;
		
		my ($zyg,$ref_reads,$var_reads);
		if ($sv_caller eq 'sniffles') {
			($zyg,$ref_reads,$var_reads) = split(':',$evidence_str);
		} elsif ($sv_caller eq 'pbsv') {
			my $all_reads;
			($zyg,$all_reads) = split(':',$evidence_str);
			($ref_reads,$var_reads) = split(',',$all_reads);
		} else {
			modules::Exception->throw("ERROR with sv_caller $sv_caller\n");
		}
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
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{status} = $pass;				
		
		if ($event_type eq 'ins') {
			if (length($var_base) > 1000) {
				$var_base = '>1000bp';
			} 	
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{seq} = $var_base;
			
		}
	}

	
    return \%sv_data;
}


