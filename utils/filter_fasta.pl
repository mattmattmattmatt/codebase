#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use File::Basename;
use modules::SystemCall;
use Cwd 'abs_path';
use File::Basename;
use Bio::SeqIO;
use Bio::SeqUtils;
use modules::Exception;

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
		"fasta=s",
		"header=s",
		"keep",
		"range=s",
		"out_file=s",
		"out_remove=s",
		"replaceNs"
	   );
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{fasta});

	   
=pod

=head1 SYNOPSIS

./filter_fasta.pl -fasta fasta_file -range contig_coords(default=whole_contig) -replaceNs replace_seq_with_Ns -out_remove file_name_for_removed -keep keep_sequence(default=remove_sequence) -header only_process_this_contig(default=all) -out_file new_fasta(default=current.new.fasta) [options]

Required flags: -fasta -range

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

bam_to_fastq.pl -> Convert BAM file to fastq to allow consistent analysisSychronise v2.1 and v2.2

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

#Remove chr1
/drive2/variantdb/trunk/utils/filter_fasta.pl -fasta test.fa -header chr1 

#Same but save chr1 in file
/drive2/variantdb/trunk/utils/filter_fasta.pl -fasta test.fa -header chr1 -keep

#Remove chr1:2-5 only
/drive2/variantdb/trunk/utils/filter_fasta.pl -fasta test.fa -header chr1 -range 2-5

#Same but save chr1:2-5 in file
/drive2/variantdb/trunk/utils/filter_fasta.pl -fasta test.fa -header chr1 -range 2-5 -keep

#Same but save chr1:2-5 in file and replace with Ns in original file
/drive2/variantdb/trunk/utils/filter_fasta.pl -fasta test.fa -header chr1 -range 2-5 -keep -replaceNs

=cut


#Read in input files
my $fasta_in = $OPT{fasta};

if (!-e $fasta_in) {
	die("Can't find file $fasta_in");
}

my $header = defined $OPT{header}?$OPT{header}:0;

my $range = defined $OPT{range}?1:0;
my ($start,$end,$length);

if ($range) {
	($start,$end) = $OPT{range} =~ /^(\d+)[,-](\d+)$/;
	$length = $end-$start+1;
} 

#Object for parsing fastas
my $sequences = Bio::SeqIO->new(-file => $fasta_in, 
								-format=>'fasta');

my $replaceNs = defined $OPT{replaceNs}?1:0;

my $keep = defined $OPT{keep}?1:0;

my $outfasta = defined $OPT{out_file}?$OPT{out_file}:$fasta_in.'.filter.fasta';

my $removedfasta = defined $OPT{out_remove}?$OPT{out_remove}:$fasta_in.'.removed';

open(OUTFASTA,">$outfasta") || modules::Exception->throw("Can't open file to write $outfasta\n");

if ($keep) {
	open(REMOVED,">$removedfasta") || modules::Exception->throw("Can't open file to write $removedfasta\n");
}


my $header_found = 0;


while ( my $dna = $sequences->next_seq ){
  	my $desc = ' ';
   	if ($dna->desc) {
   		$desc = $dna->desc;
   	} 
	
   	if ($header) {
   		if ($header eq $dna->id) {
   			my $removed_seq = my $new_seq;
   			$header_found = 1;
   			
   			
   			if ($range) {
	   			($new_seq,$removed_seq) = &process_seq($dna->seq,$start,$end,$keep);
	   			print OUTFASTA join ("\n",
									'>'.$dna->id. ' '.$desc,
									$new_seq
									) ."\n";
				if ($keep) {
					print REMOVED join ("\n",
										'>'.$dna->id. ' '.$desc ." position $start-$end was removed",
										$removed_seq
										) ."\n";
				}
   			} else {
   				$removed_seq = $dna->seq;
   				if ($keep) {
					print REMOVED join ("\n",
										'>'.$dna->id. ' '.$desc ." whole contig was removed",
										$removed_seq
										) ."\n";
				}
   			}
   		} else {
   			#Don't touch these ones
   			print OUTFASTA join ("\n",
								'>'. $dna->id. ' '.$desc,
								$dna->seq
								) ."\n"
   		}
   		
   	} else {
   		#Process every sequence here
   		my $removed_seq = my $new_seq;
   		if ($range) {
	   		($new_seq,$removed_seq) = &process_seq($dna->seq,$start,$end,$keep);
	   		print OUTFASTA join ("\n",
								'>'. $dna->id . ' '. $desc,
								$new_seq
								) ."\n";
			if ($keep) {
				print REMOVED join ("\n",
									'>'.$dna->id. ' '.$desc ." position $start-$end was removed",
									$removed_seq
									) ."\n";
			}
   		} else {
   			$removed_seq = $dna->seq;
   			if ($keep) {
				print REMOVED join ("\n",
									'>'.$dna->id. ' '.$desc ." whole contig was removed",
									$removed_seq
									) ."\n";
			}
   		}
   	}
  	
}

if ($header && !$header_found) {
	modules::Exception->throw("ERROR: Didn't find header $header in file so nothing change\n");
}

sub process_seq {
	my ($seq,$start_coord,$end_coord,$keep_flag) = @_;
	my $new_seq;
	my $removed_seq;
	if (!$keep_flag) {
		$new_seq = substr($seq,$start_coord-1,$length)
	} else {
		my $seq1 = substr($seq,0,$start_coord-1);
		my $seq2 = substr($seq,$end_coord);
		
		if ($replaceNs) {
			$new_seq = $seq1.'N'x ($length).$seq2;			
		} else {
			$new_seq = $seq1.$seq2;
		}
		if ($keep) {
			$removed_seq = substr($seq,$start-1,$length);
		}
	}
	return ($new_seq,$removed_seq);

}









