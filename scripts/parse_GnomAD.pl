#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "input=s",
	   "debug=s",
	   "outdir=s",
	   "chr=s",
	   "ref=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{input} || !$OPT{chr});



=pod

=head1 SYNOPSIS

parse_ExAC.pl -input gnomAD_file -chr chr -outdir output_directory(default=pwd) -ref reference_genome(default=hs37d5) [options]

Required flags: -input -chr

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

parse_gnomAD.pl -> Script to parse ExAC data

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./parse_gnomAD.pl -input gnomad.genomes.r2.0.2.sites.chr1.vcf -chr 1 -outdir /home/matt/work/pipeline/conf/human/GRCh37/gnomAD/2.0.2


=cut

my $gnomad = $OPT{input};
my $outdir = defined $OPT{outdir}?`readlink -f $OPT{outdir}`:`pwd`;
chomp $outdir;
my $count = 0;
my %final_snps = ();
my $ref = defined $OPT{ref}?$OPT{ref}:'hs37d5';
my $chr  = $OPT{chr};

open(SNV,">$outdir/$ref.gnomAD.overlap.snv.$chr") || die "Can't open file to snv write $outdir\n";
open(INDEL,">$outdir/$ref.gnomAD.overlap.indel.$chr") || die "Can't open indel file to write $outdir\n";


my @fh = (\*SNV,\*INDEL);


#now parse the entries
open(GNOMAD,"$gnomad") || die "Can't open file $gnomad\n";

while (<GNOMAD>) {
	next if /^#/;
	next unless /\w/;
	chomp;
	$count++;
    if ($OPT{debug}) {
		next unless $_ =~ /$OPT{debug}/;
		print "$_\n";
	}
    
    	
    	
    my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest) = split;
    
	if ($count %1000000 == 0) {
		print "$count $chr $first_coord\n";
	}
	
    if ($chr eq 'MT') {
		$chr = 'M';
	}
    	

    if ($qual !~ /\d/ && $qual ne '.') {
    	modules::Exception->throw("ERROR: Error with qual $qual format at line $_");
    }

	if ($ref =~ /N/ || $var_str =~ /N/) {
		next;
	}

	if ($var_str eq '.') {
		next;
	}

	my @vars = split(",",$var_str);
	my ($var_count_str) = $rest =~ /AC=([0-9,]+);/;
	my ($allele_count) = $rest =~ /AN=([0-9]+);/;
	
	#not sure why these cases are present
	if ($allele_count == 0) {
		next;
	}
	my @var_counts = split(",",$var_count_str);
	
	if (@var_counts != @vars) {
		modules::Exception->throw("ERROR: Problem with counts for line $_\n");
	}		
	
	my %coord_vars = ();
		
	for ( my $var_count = 0 ; $var_count < @vars ; $var_count++ ) {
	    
		my ($var_key,$var_type) = _get_variant_key(-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$vars[$var_count]);

		if ($var_key == -1) {
			print "Skip variant $chr $first_coord $vars[$var_count]\n";
			next;
		}
		my ($chr,$start,$end) = $var_key =~ /(.+):(\d+)-(\d+)/;

		my $var_count = $var_counts[$var_count];

		my $coord_key = "$chr:$start-$end";
		
		push @{$coord_vars{$var_type}{$coord_key}}, join('~',$var_key,$var_count,$allele_count,$qual);
	}
		
	for my $local_var_type (keys %coord_vars) {
		for my $coord_str (sort {my ($a_end) = $a =~ /(\d+)$/; my ($b_end) = $b =~ /(\d+)$/; $a_end <=> $b_end } keys %{$coord_vars{$local_var_type}}) {
			
			
			my @exac_strs;
			my ($chr,$start,$end) = $coord_str =~ /(.+):(\d+)-(\d+)/;
			my $bases;
			for my $exac_entry (@{$coord_vars{$local_var_type}{$coord_str}}) {
				my ($var_key,$var_count,$allele_count,$qual) = split ('~',$exac_entry);
				(undef,undef,undef,$bases) = $var_key =~ /([0-9XYM]+):(\d+)-(\d+):(.+)$/;
				my $freq = sprintf("%.6f",$var_count/$allele_count);
				push @exac_strs, $bases.':'.$freq.':'.$qual;
			}
			#need to put start here to make them unique; otherwise the overlap takes too long -> we strip it out later in processing
			my $exac_str = $start.':'.join(';',@exac_strs);
			my $fh;
			if ($local_var_type eq 'SNV') {
				$fh = $fh[0];  
			} elsif ($local_var_type eq 'DEL') {
				$fh = $fh[1];
			} elsif ($local_var_type eq 'INS') {
				$fh = $fh[1];
			} else {
				modules::Exception->throw("ERROR: Var type $local_var_type doesn't match");
			}
			#print join ("\t", $chr,$start,$end,$exac_str) ."\n"; 
			print $fh join ("\t", $chr,$start,$end,$exac_str) ."\n"; 
		}
	}
}
for my $fh (@fh) {
	close $fh;
}

#Gets variant key from either vcf or vep; standardises naming for loading into data structure
sub _get_variant_key {
	 my @args = @_;
	
	 my %args = @args;


    my @required_args = (
    					-chrom,
    					-first,
    					-ref_seq,
    					-var_seq,
    					);

    foreach my $required_arg (@required_args){

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set $args{$required_arg}");
		} 
    }
    
    my $ref = $args{-ref_seq};
    my $var = $args{-var_seq};
    my $first_coord = $args{-first};
    my $chr = $args{-chrom};
    
    my $start_coord = my $end_coord = my $bases;
    my $length_ref = length($ref);
    my $length_var = length($var);
    my $var_type;
    
	if ($length_ref > $length_var) {
		$var_type = 'DEL';
		$start_coord = $first_coord + 1;
		$end_coord = $start_coord + $length_ref - $length_var - 1;				
		my $del_length = $length_ref - $length_var;
		
		$bases = '-'. substr($ref,1,$del_length);
		#print "VCF R $ref L $del_length B $bases\n";
	} elsif ($length_ref < $length_var) {
		#Add the ref length and var length difference to the coord 
		#$start_coord = $end_coord = $first_coord + 1;
		$var_type = 'INS';
		$start_coord = $end_coord = $first_coord;
		my $ins_length = $length_var - $length_ref;
		$bases = '+'.substr($var,1,$ins_length);
	} elsif ($length_ref == 1) {
		$var_type = 'SNV';
		$start_coord = $end_coord = $first_coord;
		$bases = $ref . '->' .$var;
	} else {
		#Here we have special case where snvs are >1 length as indels and snvs are combined into single line
		my @var_bases = split('',$var);
		my @ref_bases = split('',$ref);
		my $same = 1;
		for ( my $count = 1 ; $count < @var_bases ; $count++ ) {
		    $same = 0 if $ref_bases[$count] != $var_bases[$count];
		}
		if ($same) {
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			$bases = $ref_bases[0] . '->' .$var_bases[0];
		} else {
			return (-1,-1); #don't handle this case (e.g. AT->GC)
		}
    	
    } 
    
	my $var_key = $chr . ':'.$start_coord .'-'.$end_coord .':' .$bases;
	return($var_key,$var_type);
}




