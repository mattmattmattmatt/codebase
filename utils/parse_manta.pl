#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);
use Bio::SearchIO;
use modules::Exception;
use Cwd 'abs_path';

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "manta_vcf=s",
	   "no_filter",
	   "outfile_base=s",
	   "outdir=s",
	   "joint",
	   "no_human"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{manta_vcf});



=pod

=head1 SYNOPSIS

parse_manta.pl -manta_vcf vcf_in -no_filter do_not_filter_failed_results -outfile_base basename_for_outputs -joint joint_vcf(not_somatic) [options]

Required flags: -manta_vcf

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_manta.pl -> Parse manta vcf and generate text coord files

=head1 DESCRIPTION

Feb 5th, 2015

Manta parser

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

>blast_template.pl -fasta_in file.fa -blast_in blast.bls -fasta_out fasta_new.fa

=cut


my $manta_vcf = $OPT{manta_vcf};
my $somatic = defined $OPT{joint}?0:1; #default is somatic
my $outbase = defined $OPT{outfile_base}?$OPT{outfile_base}:'manta';

my $outdir = defined $OPT{outdir}?$OPT{outdir}:`pwd`;
chomp $outdir;
$outdir = abs_path($outdir);

if ( !-e $manta_vcf ) {
	modules::Exception->throw("File $manta_vcf doesn't exist");	
}

my $manta_sv = &parse_manta_result($manta_vcf);

&print_files($manta_sv,$outbase,'manta');



sub print_files {
	my ($results,$file,$sv_caller) = @_;
	open(FILE,">$outdir/${file}.txt") || modules::Exception->throw("Can't open file $file\n");
	open(TRA,">$outdir/${file}.tra") || modules::Exception->throw("Can't open file $file\n");
	for my $sv_type (qw(del dup inv ins bnd)) {
		for my $coord_str1 (sort {my ($a_chr,$a_coord) = $a =~ /([0-9XYM]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9XYM]+):(\d+)/; $a_chr <=> $b_chr || $a_coord <=> $b_coord }  keys %{$results->{$sv_type}}) {
				for my $coord_str2 ( keys %{$results->{$sv_type}{$coord_str1}}) {
					my ($chr1,$coord1) = split(':',$coord_str1);
					my ($chr2,$coord2) = split(':',$coord_str2);
					
					my $anno_str;
					
					
					if ($somatic) {
						
							$anno_str = join("^^^",'SVCALLER='.$sv_caller,
												'SVTYPE='.$sv_type,
												'ID='.$results->{$sv_type}{$coord_str1}{$coord_str2}{id},
												'PE='.$results->{$sv_type}{$coord_str1}{$coord_str2}{pe},
											  'SR='.$results->{$sv_type}{$coord_str1}{$coord_str2}{sr},
											  'LEN='.$results->{$sv_type}{$coord_str1}{$coord_str2}{length},
											  'QUAL='.$results->{$sv_type}{$coord_str1}{$coord_str2}{qual},
											  'MANTA_GFF='.$results->{$sv_type}{$coord_str1}{$coord_str2}{gff},
											  'REF_EVIDENCE='.$results->{$sv_type}{$coord_str1}{$coord_str2}{norm},
											  'TUMOUR_EVIDENCE='.$results->{$sv_type}{$coord_str1}{$coord_str2}{tumour},
											  'COORD='.$coord_str1.'+'.$coord_str2
											 );
					} else {
						$anno_str = join("^^^",'SVCALLER='.$sv_caller,
												'SVTYPE='.$sv_type,
												'ID='.$results->{$sv_type}{$coord_str1}{$coord_str2}{id},
												'PE='.$results->{$sv_type}{$coord_str1}{$coord_str2}{pe},
											  'SR='.$results->{$sv_type}{$coord_str1}{$coord_str2}{sr},
											  'LEN='.$results->{$sv_type}{$coord_str1}{$coord_str2}{length},
											  'QUAL='.$results->{$sv_type}{$coord_str1}{$coord_str2}{qual},
											  'MANTA_GFF='.$results->{$sv_type}{$coord_str1}{$coord_str2}{gff},
											  'SAMPLE_COUNT='.$results->{$sv_type}{$coord_str1}{$coord_str2}{sample_count},
											  'SAMPLE_STR='.$results->{$sv_type}{$coord_str1}{$coord_str2}{sample_str},
											  'COORD='.$coord_str1.'+'.$coord_str2
											  );
					}
					
					
					$anno_str .= '^^^ins='.$results->{$sv_type}{$coord_str1}{$coord_str2}{seq} if $sv_type eq 'ins';						 
					
					my $second_coord = $sv_type eq 'tra'?$coord1:$coord2;
					if (!$OPT{no_human}) {
						$chr1 =~ s/chr//;
						$chr2 =~ s/chr//;
					}
					
					
					
					if ($sv_type eq 'bnd') {
						print TRA join ("\t",
									$chr1,
									$coord1,
									$chr2,
									$coord2,
									$anno_str,
									
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

sub parse_manta_result {
    my ($output_vcf) = @_;
    
    
    open(my $OUTPUT, $output_vcf) || modules::Exception->throw("ERROR: Can't open output file $output_vcf");
	my %sv_data = ();
	
	while (<$OUTPUT>) {
		chomp;
		next if /^#/;
		
		my ($chr1,$start1,$id,$ref_base,$var_base,$score,$pass_status,$gff_str,$evidence_fields,$evidence_str_ref,$evidence_str_tumour);
		my @sample_gts = ();
		my $hom_count = my $ref_count = my $het_count = 0;
		if ($somatic) {
			($chr1,$start1,$id,$ref_base,$var_base,undef,$pass_status,$gff_str,$evidence_fields,$evidence_str_ref,$evidence_str_tumour) = split("\t");
		} else {
			($chr1,$start1,$id,$ref_base,$var_base,$score,$pass_status,$gff_str,$evidence_fields,@sample_gts) = split("\t");
		}
		
		if ($pass_status ne 'PASS') {
			next unless $OPT{no_filter};
		}
		
		
		$chr1 =~ s/chr// unless $OPT{no_human};
		my $chr2 = $chr1;
		
		if (!$OPT{no_human}) {
			next unless $chr1 =~ /^([0-9XYM]+)/;
			next unless $chr2 =~ /^([0-9XYM]+)/;
			next if $chr1 =~ /random/;
			next if $chr2 =~ /random/;
		}
		
		my ($start2) = $gff_str =~ /END=(\d+);/; 
		my ($length) = $gff_str =~ /SVLEN=\-?(\d+)/;
		
		my $sr =  my $pr = my $sr_count = my $pe_count = 0;
		my $qual;
		
		if ($somatic) {
			if ($evidence_fields eq 'PR:SR') {
				($pr,$sr) = split(':',$evidence_str_tumour);
				(undef, $sr_count) = split(",",$sr);
				(undef, $pe_count) = split(",",$pr);				
			} else {
				(undef, $pe_count) = split(",",$evidence_str_tumour);
			}
			($qual) = $gff_str =~ /SOMATICSCORE=(\d+)/;
			
		} else {
			
			
			
			my @pr_counts;
			my @sr_counts;
			for my $sample_str (@sample_gts) {
				my @fields = split(":",$sample_str);
				if ($fields[0] eq '0/0') {
					$ref_count++;
				} elsif ($fields[0] eq '1/1') {
					$hom_count++;
				} else {
					$het_count++;
				}
					
				if ($evidence_fields =~ /PR:SR/) {
					my (undef, $sr_tmp) = split(",",$fields[-1]);
					my (undef, $pe_tmp) = split(",",$fields[-2]);
					push @sr_counts, $sr_tmp;
					push @pr_counts, $pe_tmp;					
				} else {
					my (undef, $pe_tmp) = split(",",$fields[-1]);
					push @pr_counts, $pe_tmp;
				}
			}
			
			$sr_count = join(",",@sr_counts);
			$pe_count = join(",",@pr_counts);
			
			$qual = $score;
			
		}
		
		
		my ($event_type_uc) = $gff_str =~ /SVTYPE=([A-Z]+)/;
		my $event_type = lc $event_type_uc;
		
		
		if ($event_type eq 'ins') {
			$start2 = $start1;
		}
		
		if ($event_type eq 'bnd') {
			($chr2,$start2) = $var_base =~ /([0-9a-zA-Z\._]+):(\d+)/;
			
			#print "VB $var_base C $chr2 S $start2\n";
			
			if ($chr1 eq $chr2) {
				$event_type = "inv";
				$length = $start2-$start1;
			}
			
			if ($chr2 !~ /^([0-9XYM]+)/) {
				next unless $OPT{no_human};
			}
			next if $chr2 =~ /random/;
			
			#Don't record SVs twice
			if (exists $sv_data{$event_type}{"$chr2:$start2"}) {
				next;
			}
		}
		
		
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{pe} = $pe_count;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{sr} = $sr_count;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{length} = $length;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{gff} = $gff_str;	
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{id} = $id;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{qual} = $qual;
		$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{gff} = $gff_str;
		
		if ($somatic) {
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{norm} = $evidence_str_ref;
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{tumour} = $evidence_str_tumour;
		} else {
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{sample_count} = "Ref=".$ref_count.":Hom=".$hom_count.":Het=".$het_count;
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{sample_str} = join('^^^',@sample_gts);
		}
		
		
		
		if ($event_type eq 'ins') {
			my ($seq) = $gff_str =~ /SEQ=([A-Z]+)/;
			$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{length} = length($seq) if $sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{length} !~ /\d/;
			$seq = 'TOO_LONG' if length($seq) > 1000;
			if ($seq =~ /[AGCT]/) {
				$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{seq} = '+'.$seq;
			} else {
				$sv_data{$event_type}{"$chr1:$start1"}{"$chr2:$start2"}{seq} = '+'.$var_base;
			}
		}
	}


	
    return \%sv_data;
}


