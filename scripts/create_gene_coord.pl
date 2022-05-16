#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use vars qw(%OPT);
use Bio::SearchIO;
use List::MoreUtils qw(uniq);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "exon_file=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || $OPT{man} || !$OPT{exon_file});



=pod

=head1 SYNOPSIS

create_gene_coord.pl -exon_file exon_file [options]

Required flags: -exon_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

create_gene_coord.pl -> Create gene coord file for overlapping SVs

=head1 DESCRIPTION

Feb 5th, 2015


=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

>create_gene_coord.pl -exon_file  /drive2/variantdb/trunk/conf/human/GRCh37/exon/040614/GRCh37.exon.overlap.1 > /drive2/variantdb/trunk/conf/human/GRCh37/gene/040614/GRCh37.gene.overlap.1

=cut

my $gene;
my $start;
my $end;
open(FILE,$OPT{exon_file}) || modules::Exception->throw("Can't open directory $OPT{exon_file}");


while (<FILE>) {
  chomp;
  my @fields = split();
  my ($new_gene) = $fields[3] =~ /(ENSG\d+)/;
  if ($new_gene ne $gene) {
    print join(" ",
                $fields[0],
                $start,
                $end,
                $gene) ."\n" if $gene =~ /ENS/;
    $start = $fields[1];
  }
  $end = $fields[2];
  $gene = $new_gene;
}

