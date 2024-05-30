#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use modules::Utils;
use File::Basename;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	"bam=s",
	   	"cellmap=s",
	   	"celltype=s",
	   	"samtools=s",
	   	"threads=i",
	   	"no_run"
	   	);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{bam} || !$OPT{cellmap});



=pod

=head1 SYNOPSIS

split_bam_by_celltype.pl -bam bam_with_barcode_readgroups -cellmap tsv_with_barcode_cell_type -celltype run_only_on_this_celltype(default=all) -samtools samtools_path(default=samtools) -no_run print_commands_to_run -threads threads_to_use(default=16)-help -man [options]

Required flags: 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

split_bam_by_celltype.pl 

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./split_bam_by_celltype.pl -bam /drive3/work/Tapestri/1912/1912.tube1.cells.bam -cellmap /drive3/work/Tapestri/1912/group_file_cat1.tsv -celltype Rogue

=cut

my $bam = $OPT{bam};
my $cellmap = $OPT{cellmap};
my $basename = dirname($cellmap);
my $threads = defined $OPT{threads}?$OPT{threads}:16;

if ( !-e $bam ) {
	modules::Exception->throw("File $bam doesn't exist");
}

if ( !-e $cellmap ) {
	modules::Exception->throw("File $cellmap doesn't exist");
}

my $samtools = defined $OPT{samtools}?$OPT{samtools}:'samtools';

my $celltype = defined $OPT{celltype}?$OPT{celltype}:'all';

open(MAP,"$cellmap") || modules::Exception->throw("Can't open file $cellmap\n");

my %cellmap = ();

while (<MAP>) {
	chomp;
	my ($barcode,$type) = split("\t");
	$type =~ s/ //g;
	$barcode =~ s/ //g;
	push @{$cellmap{$type}}, $barcode;
}


my @commands = ();


for my $type ( keys %cellmap ) {
    if ($celltype eq 'all') {
    	my $file = $basename.'/'.$type."_rg.txt";
    	open(FILE,">$file") || modules::Exception->throw("Can't open file to write $file\n");
    	print FILE join("\n",
    					@{$cellmap{$type}}
    					) . "\n";
    	(my $out = $bam) =~ s/.bam/_${type}.bam/; 
    	push @commands, "$samtools view $bam -@ $threads -R $file -b -o $out";
    	push @commands, "$samtools index -@ $threads $out";
    } elsif ($celltype eq $type) {
    	my $file = $basename.'/'.$type."_rg.txt";
    	open(FILE,">$file") || modules::Exception->throw("Can't open file to write $file\n");
    	print FILE join("\n",
    					@{$cellmap{$type}}
    					) . "\n";
    					
    	(my $out = $bam) =~ s/.bam/_${type}.bam/; 
    	push @commands, "$samtools view $bam -@ $threads -R $file -b -o $out";
    	push @commands, "$samtools index -@ $threads $out";
    }
}

for my $command (@commands) {
	print $command ."\n";
	`$command` unless $OPT{no_run};
}

#print Dumper \@commands;
