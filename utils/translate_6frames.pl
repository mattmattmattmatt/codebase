#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use modules::Exception;
use vars qw(%OPT);
use Bio::SearchIO;
use Bio::SeqIO;
use List::Util;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "nt_in=s"
	   
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{nt_in});



=pod

=head1 SYNOPSIS

blast_template.pl -blast_in input_blast_in_bsl -fasta_in fasta_file -fasta_out new_fasta_name(deault=fasta_in.new.fa) [options]

Required flags: -blast_in -fasta_in

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

blast_template.pl -> Parse blast output and append hit info to original fasta header

=head1 DESCRIPTION

Feb 5th, 2015

Blast parsing template

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

>blast_template.pl -fasta_in file.fa -blast_in blast.bls -fasta_out fasta_new.fa

=cut

#Read in input files
my $nt_in = $OPT{nt_in};

#Check the files exist
if (!-e $nt_in) {
	die("Can't find file $nt_in");
}

#Object for parsing fastas
my $nt_obj = Bio::SeqIO->new(-file => $nt_in, 
								-format=>'fasta');

(my $nt_base = $nt_in) =~ s/\.fa*//;

open(FILE, ">${nt_base}_translate.fasta") || modules::Exception->throw("Can't open file");
open(FILE2,">${nt_base}_translate_best_Mstart.fasta") || modules::Exception->throw("Can't open file $\n");

my %seq = ();

#Read the protein entries
while( my $fasta_seq = $nt_obj->next_seq ) {
	my $header = $fasta_seq->id;
	my $seq = $fasta_seq->seq;
	(my $seq2 = $seq) =~ s/^[A-Z]//;
	(my $seq3 = $seq2) =~ s/^[A-Z]//;
	my $revseq = revcom($seq);
	(my $revseq2 = $revseq) =~ s/^[A-Z]//;
	(my $revseq3 = $revseq2) =~ s/^[A-Z]//;
	
	my $seq_f1 = Bio::Seq->new(-seq => $seq, -alphabet => 'dna');
	my $prot_f1 = $seq_f1->translate(); 
		
	my $seq_f2 = Bio::Seq->new(-seq => $seq2, -alphabet => 'dna');
	my $prot_f2 = $seq_f2->translate();
	
	my $seq_f3 = Bio::Seq->new(-seq => $seq3, -alphabet => 'dna');
	my $prot_f3 = $seq_f3->translate();
	
	my $revseq_f1 = Bio::Seq->new(-seq => $revseq, -alphabet => 'dna');
	my $revprot_f1 = $revseq_f1->translate(); 
		
	my $revseq_f2 = Bio::Seq->new(-seq => $revseq2, -alphabet => 'dna');
	my $revprot_f2 = $revseq_f2->translate();
	
	my $revseq_f3 = Bio::Seq->new(-seq => $revseq3, -alphabet => 'dna');
	my $revprot_f3 = $revseq_f3->translate();
	
	print FILE join("\n",
					">".$header.'_f1',
					$prot_f1->seq,
					">".$header.'_f2',
					$prot_f2->seq,
					">".$header.'_f3',
					$prot_f3->seq,
					">".$header.'_r1',
					$revprot_f1->seq,
					">".$header.'_r2',
					$revprot_f2->seq,
					">".$header.'_r3',
					$revprot_f3->seq,
					) ."\n";	
	
	my @f1_count = $prot_f1->seq =~ /\*/g;
	my @f2_count = $prot_f2->seq =~ /\*/g;
	my @f3_count = $prot_f3->seq =~ /\*/g;
	my @r1_count = $revprot_f1->seq =~ /\*/g;
	my @r2_count = $revprot_f2->seq =~ /\*/g;
	my @r3_count = $revprot_f3->seq =~ /\*/g;
	
	my $best_seq = $prot_f1->seq;
	my $min = scalar @f1_count;
	my $header_tag = '_f1';
	
	if (scalar @f2_count < $min) {
		$best_seq = $prot_f2->seq;
		$min = scalar @f2_count;
		$header_tag = '_f2';
	}	
	
	if (scalar @f3_count < $min) {
		$best_seq = $prot_f3->seq;
		$min = scalar @f3_count;
		$header_tag = '_f3';
	}	
	
	if (scalar @r1_count < $min) {
		$best_seq = $revprot_f1->seq;
		$min = scalar @r1_count;
		$header_tag = '_r1';
	}	
		
	if (scalar @r2_count < $min) {
		$best_seq = $revprot_f2->seq;
		$min = scalar @r2_count;
		$header_tag = '_r2';
	}	
	
	if (scalar @r3_count < $min) {
		$best_seq = $revprot_f3->seq;
		$min = scalar @r3_count;
		$header_tag = '_r3';
	}	
	
	my ($m_seq) = $best_seq =~ /(M.*)/;
	
	my $mlength = length($m_seq);
	my $full_length = length($best_seq);
	
	if (length($m_seq) > 1) {
		print FILE2 join ("\n",
						">".$header.$header_tag,
						$m_seq
						) . "\n";
	}
	
	
	
}


sub revcom {
	my ( $seq ) = @_;
    my $revcomp = reverse($seq);
  	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
  	return $revcomp;
}
