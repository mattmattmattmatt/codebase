#! /usr/bin/perl -w

use strict;
use modules::Adaptors::SNV;
use modules::Adaptors::Filter;
use modules::Adaptors::SNV_Filter;
use modules::Adaptors::Variant;
use modules::Adaptors::Variant_Filter;
use modules::Adaptors::BulkInsert;
use modules::VariantXML;
use modules::Exception;
use modules::Pipeline;
use modules::Utils;
use modules::QualityEncoding;
use modules::SystemCall;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(
			\%OPT,                   
			"help|h",
			"man|m",                 
			"runid=i",
			"snv_pileup=s",          
			"indel_pileup=s",
			"overlap_outfile_snv=s", 
			"overlap_outfile_indel=s",
			"infile_snv=s",          
			"infile_indel=s",
			"chr=s",
			"writeDB=i",
			"software=s",             
			"coord=i",
			"debug",
);
pod2usage( -verbose => 2 ) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{snv_pileup} || !$OPT{indel_pileup} || !$OPT{overlap_outfile_snv} || !$OPT{overlap_outfile_indel} || !$OPT{infile_snv} || !$OPT{infile_indel} );


=pod

=head1 SYNOPSIS

db_variants.pl -runid runid -snv_pileup snv_pileup_file  -indel_pileup indel_pileup_file  -software software_used(default=samtools) -writeDB(default=1) -overlap_outfile_snv output_overlap_snv_file -infile_snv input_overlap_file -debug output_commands [options]

Required flags: -runid -snv_pileup -indel_pileup -overlap_outfile_snv -overlap_outfile_indel -infile_snv -infile_indel 

=head1 OPTIONS

	-help  brief help message
	
	-man   full documentation

=head1 NAME

db_variants.pl -> Add all variants to the db

=head1 DESCRIPTION

Feb 15, 2011

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./db_variants.pl

=cut

my $delim             = '^^^';
my $runid             = $OPT{runid};
my $write_to_db       = defined $OPT{writeDB} ? $OPT{writeDB} : 1;
my $pileup_file_snv   = $OPT{snv_pileup};
my $pileup_file_indel = $OPT{indel_pileup};
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $software = defined $OPT{software}?$OPT{software}:'samtools';
my $tumour = 0;

if (modules::Pipeline::get_tumour_flag(-sample_name=>$sample_name)) {
	$tumour = 1;
} 


my $pileup_normal_file_snv;
my $pileup_normal_file_indel;

if ($tumour) {
	($pileup_normal_file_snv = $pileup_file_snv) =~ s/pileup/normal.pileup/;
	($pileup_normal_file_indel = $pileup_file_indel) =~ s/pileup/normal.pileup/;
	if ( !-e $pileup_normal_file_snv ) {
		modules::Exception->throw("File $pileup_normal_file_snv doesn't exist");	
	}
	if ( !-e $pileup_normal_file_indel ) {
		modules::Exception->throw("File $pileup_normal_file_indel doesn't exist");	
	}
} 

my $output_snv = $OPT{overlap_outfile_snv};
my $input_snv  = $OPT{infile_snv};
if ( !-e $input_snv ) {
	modules::Exception->throw("File $input_snv doesn't exist");
}

my $output_indel = $OPT{overlap_outfile_indel};
my $input_indel  = $OPT{infile_indel};
if ( !-e $input_indel ) {
	modules::Exception->throw("File $input_indel doesn't exist");
}

#Get the quality stores from a previous step
my %snv_info = ();
my $run;
my $filter_snv;
my $filter_deletion;
my $filter_insertion;
my $filter_mutation;

if ($write_to_db) {
	($run) = modules::Adaptors::Run->search( 'id' => $runid );
	unless (    defined $run    )
		{
			modules::Exception->throw(
						 "Unable to retrieve run and filter objects from database");
		}
}

my $chr_flag = defined $OPT{chr}?$OPT{chr}:'all';

#Open filter_vcf snv file to get snv details like CLR score (eg C013_sg1_tumour1_78.filter_vcf.match.snv)
open( SNV, "$input_snv" ) || modules::Exception->throw("Can't open file $input_snv\n");

while (<SNV>) {
	#1	77636	77636	SNV^^^C->T^^^Q214^^^REF_NORM207^^^SOMATIC^^^D36	
	chomp;

  	# Then, strip out chr from chromsome column at source, if it exists
    s/^chr//i;

	next unless /^[0-9XYM]/;
	my @fields = split("\t");

	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $fields[0];
	}

	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $fields[1];
	}
	
	my ( $type, $base_change, $quality, $clr, $snv_class, $depth ) = split( '\^\^\^', $fields[3] );

	if ($base_change =~ /N/i) {
		next;
	}

    # tidy the base change, as GATK output sometimes records this like AGT->CGT.  Take the first base of each substring - this seems to always be the variant one.
    $base_change =~ s/^(.{1})[^\-]+/$1/;
    $base_change =~ s/^([^\-]+)\-\>(.{1}).*/$1\-\>$2/;

	$quality =~ s/Q//;
	$depth   =~ s/D//;

	if ($tumour) {	
		$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{clr} = $clr;
		$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{class} = $snv_class;
	}
	$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{quality} = $quality;
	$snv_info{ $fields[0] }{ $fields[1] }{$base_change}{quality_reads} = $depth;
	
}

my $sys = modules::SystemCall->new();

#If it's the first chromosome remove the file as we're appending; should happen with db cleanup but if run fails in this step it wouldn't work
if ($chr_flag eq '1' || $chr_flag eq 'all') {
	$sys->run("rm $output_snv") if (-e $output_snv);
	$sys->run("rm $output_indel") if (-e $output_indel);
}

#Append since we're running on per chr basis
open( SNVOUT, ">>$output_snv" ) || modules::Exception->throw("Can't open file $output_snv\n");
open( INDELOUT, ">>$output_indel" ) || modules::Exception->throw("Can't open file $output_indel\n");

#Sometimes we may not have snvs and indels for a particular chr
if ( !-e $pileup_file_snv && !-e $pileup_file_indel ) {
	modules::Exception->warning("SKIP STEP db_variants: input pileup file $pileup_file_snv and $pileup_file_indel doesn't exist\n");
	exit;
}


my $count         = 0;
my $print_counter = 0;

open( my $PF_SNV, $pileup_file_snv ) || modules::Exception->throw("Unable to open file [$pileup_file_snv]");

#my @snv_filter_inserts = ();

my @snv_rows;
my @snvfilter_rows;
my @mutant_snv_rows;
my %mutant_snv_data;
my %snv_dbid;
my %seen_chr;

my $homozygous_min_snv_frequency   = 0.8;
my $heterozygous_min_snv_frequency = 0.3;

#Get the remaining snv infor from the pileup file generated in merge_vcf (eg C013_sg1_tumour1_78.merge_vcf.out.snv.pileup)
while (<$PF_SNV>) {
 
	#10	61419	G	11	,,.aa.a,,A,	CFCFFDFFFCJ
	print STDERR "Processing base " . $print_counter . "\n" if ( $print_counter % 10000000 == 0 && $print_counter != 0 ); 
	$print_counter++;

  	chomp;


    # Then, strip out chr from chromosome column at source, if it exists
  	s/^chr//i;

	my @cols = split /\t/;
	
	if ( $OPT{chr} ) {
		next unless $OPT{chr} =~ /$cols[0]/;
	}
	
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $cols[1];
	}
	
	#Skip if it failed the filter_vcf step
	next unless exists $snv_info{ $cols[0] }{ $cols[1] } ;	

	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $cols[1];
	}
	my $ref_base = uc( $cols[2] );

	if ( $ref_base eq 'n' || $ref_base eq 'N' || $ref_base eq '*' || $cols[3] =~ /N/i) {    # skip this base if it is an N or an indel
		    #	print STDERR "-- skip N base\n";
		next;
	}

	my $orig_snv_string = $cols[4];

	if (!defined $orig_snv_string) { #Avoid indel realign cases where there is no data
		next;
	}

	#Convert ref base symbols to the reference base for pileup object
	$cols[4] =~ s/[\.\,]/$ref_base/g;
	$cols[4] =~ s/[\~\$]//g;
	$cols[4] =~ s/\^.//g;

	#Get rid of indel calls 
	while ( $cols[4] =~ /([\+\-])(\d+)/ )
	{
		my $sign   = "\\" . $1;
		my $length = $2;

		$cols[4] =~ s/$sign\d+.{$length}//;
	}

	my $read_depth = $cols[3];

	$cols[4] =~ s/(.)/$1\$delim/g;
	my @bases = split /\$delim/, $cols[4];
	my %bases;

	# call a consensus base for this snv
	foreach my $base (@bases)
	{
		next if $base eq '*';    #Skip deletions annotated this way
		$bases{ uc($base) }++;
	}
	my @sorted_bases = sort { $bases{$b} <=> $bases{$a} } keys %bases;

	my $ref_base_freq =
	  defined $bases{$ref_base} ? sprintf( "%.3f",$bases{ uc($ref_base) } / $read_depth) : '0.0';

	my @variant_bases;
	foreach my $sorted_base (@sorted_bases) {
		if ( defined $ref_base && ( $sorted_base eq '.' || uc($sorted_base) eq $ref_base ) ) {
			next;
		} else {
			push @variant_bases, uc($sorted_base);
		}
	}

	if ( scalar @variant_bases == 0 ) {    # Not a variant
		next;
	}

	my $commonest_variant_base      = $variant_bases[0];
	my $next_commonest_variant_base = $variant_bases[1] if scalar @variant_bases > 1;

	my $commonest_variant_base_freq = sprintf( "%.3f", $bases{$commonest_variant_base} / $read_depth );
	my $next_commonest_variant_base_freq = sprintf( "%.3f", $bases{$next_commonest_variant_base} / $read_depth ) if defined $next_commonest_variant_base;

	my $commonest_base = $sorted_bases[0];
	my $var_base       = uc($commonest_variant_base);

	$cols[4] = $orig_snv_string;

	$count++;

	$cols[5] =~ s/(.)/$1\t/g;
	my @chars = split /\t/, $cols[5];

	#    pop @chars;

	my @qualities;
	for ( my $i = 0 ; $i < scalar @chars ; $i++ ) {
		$qualities[$i] = ord( $chars[$i] ) - 33;
	}

	$cols[5] = join( "\:", @qualities );

	#     Homozygous snvs

	my $zygosity;

	if ( defined $commonest_variant_base_freq && $commonest_variant_base_freq >= $homozygous_min_snv_frequency ) {
		$zygosity = 'hom';
	} elsif ((defined $ref_base_freq
			&& defined $commonest_variant_base_freq
			&& ( $ref_base_freq + $commonest_variant_base_freq ) >= $homozygous_min_snv_frequency
			&& $ref_base_freq >= $heterozygous_min_snv_frequency
			&& $commonest_variant_base_freq >= $heterozygous_min_snv_frequency )
		    || (defined $next_commonest_variant_base_freq && defined $commonest_variant_base_freq
			&& ( $commonest_variant_base_freq + $next_commonest_variant_base_freq ) >= $homozygous_min_snv_frequency
			&& $next_commonest_variant_base_freq >= $heterozygous_min_snv_frequency 
			&& $commonest_variant_base_freq >= $heterozygous_min_snv_frequency)) {
		#Complex 'het' definition
		$zygosity = 'het';
	} else {
		$zygosity = 'other';
	}

	#Get the quality score
	my $base_change_key = $cols[2] . '->' . $var_base;
	if ( !exists $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key} ) {
		#Here mpileup calls a different varbase then our pileup module so go with the mpileup change
		for my $mpileup_base_change (keys %{ $snv_info{ $cols[0] }{ $cols[1] } } ) {
			($var_base) = $mpileup_base_change =~ /(\S)$/;
		}
		$base_change_key = $cols[2] . '->' . $var_base;
	}

	if (! exists $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{quality}) {
#	    modules::Exception->warning("ERROR: $cols[0] $cols[1] $base_change_key ----- THERE SHOULD NOT BE _ANY_ OF THESE, but warning for now for debug\n");
		modules::Exception->throw("ERROR: $cols[0] $cols[1] $base_change_key\n");
	}

	my $quality_score = $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{quality};
	my $clr_str = defined $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{clr}?$snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{clr}:'N/A';
	my $snv_class = defined $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{class}? $snv_info{ $cols[0] }{ $cols[1] }{$base_change_key}{class}:'N/A';

	### Write to database
	my $median = modules::Utils->median( \@qualities );

	#Different depending on whether tumour or normal
	my $normal_base_str = "N/A"; #set to N/A by default in case of NO_CLR 
	my $tumour_base_str;
	
	#Slightly different for normal and tumour samples; if tumour have to get the normal pileup string 
	if ($tumour) {
		#Get the normal base str
		$tumour_base_str = $cols[4];
		my @lines = split("\n",`grep -w $cols[1] $pileup_normal_file_snv`);
		
		if (@lines) {
			my @fields = split("\t",$lines[0]);
			$normal_base_str = $fields[4];
		} 
		
		
		
	} else {		
		$tumour_base_str = 'N/A';
		$normal_base_str = $cols[4];
	
	}
	
	my  $snv_str = join( "\t",
						$cols[0],
						$cols[1],
						$cols[1],
						'SNV'. $delim . $cols[2] . '->'. $var_base
						  . $delim
						  . $snv_class
						  . $delim . 'Q'
						  . $quality_score 
						  . $delim 
						  . $clr_str 
						  . $delim . 'D' 
						  . $read_depth 
						  . $delim 
						  . $zygosity
						  . $delim
						  . $bases{$commonest_variant_base} . '/'
						  . $read_depth
						  . $delim . 'MQ' 
						  . $median 
						  . $delim . 'TBS:'
						  . $tumour_base_str
						  . $delim . 'NBS:'
						  . $normal_base_str
						  . $delim
						  . $cols[5]
						  );

	print SNVOUT $snv_str . "\n";

	#default
	my $clr_score = 0;
			
	#if it's a tumour sample
	if ($clr_str =~ /(\d+)/) {
		$clr_score = $1;
	}
			
	#Make sure there's no single quotes as it will break the mysql inserts
	$normal_base_str =~ s/'//g;
	$normal_base_str =~ s/"//g;
	$tumour_base_str =~ s/'//g;
	$tumour_base_str =~ s/"//g;

	#Now create the data struct for the xml
	my $mutant_snv_key = $cols[0] . ':' . $cols[1];
	my $end = $cols[1];
	$normal_base_str =~ s/&/%/g; #Doesn't work in xml so hack as we don't use the ASCII quality
	$normal_base_str =~ s/</=/g;
	$normal_base_str =~ s/>/=/g;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{ref_bases} = $ref_base;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{var_bases} = $var_base;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{ref_base_freq} = $ref_base_freq;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{var_base_freq} = $commonest_variant_base_freq;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{median_quality_score} = $median;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{snv_score} = $quality_score;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{normal_base_string} = $normal_base_str;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{tumour_base_string} = $tumour_base_str;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{clr_score} = $clr_score;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{snv_class} = $snv_class;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{read_depth} = $read_depth;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{var_depth} = $bases{$commonest_variant_base};
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{zygosity} = $zygosity;
	$mutant_snv_data{$mutant_snv_key}{$end}{$var_base}{pass} = 1;
		
}


#Whole genome objects
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');


my $file_name_snv = $sample_name . '_' . $runid . '.db_variants.snv.xml';
$var_xml->create_var_xml(-file_name=>$file_name_snv,-data=>\%mutant_snv_data,-chr=>$chr_flag);
if ($chr_flag eq 'Y' || $chr_flag eq 'all') {
	my $full_xml_file = $run_dir . '/conf/'.$file_name_snv;
	$var_xml->split_xml(-file_name=>$full_xml_file);
}


close($PF_SNV);
close(SNVOUT);

#Now indels...
	
#First parse the pileup to get the base info; do this first as we don't always see support from mpileup output
my %indel_info = ();
open( my $PF_INDEL, $pileup_file_indel ) || modules::Exception->throw("Unable to open file [$pileup_file_indel]");
while(<$PF_INDEL>) {

	#10	82307	C	2	.-1G,-1g	JB
	print STDERR "Processing base " . $print_counter . "\n" if ( $print_counter % 10000000 == 0 && $print_counter != 0 ); 
	$print_counter++;

       chomp;


       # Then, strip out chr from chromosome column at source, if it exists
       s/^chr//i;

	my ($chr,$start,$ref,$total_depth,$base_str,$qual_str) = split /\t/;

	next if $total_depth == 0;
	
	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $chr;
	}
	
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $start;
	}
	
	#Account for start of indel being one more than vcf coordinate
	my $indel_start = $start + 1;
	$indel_info{$chr}{$indel_start}{base_str} = $base_str;
	$indel_info{$chr}{$indel_start}{qual_str} = $qual_str;
	$indel_info{$chr}{$indel_start}{total_reads} = $total_depth;
}




my @insertion_rows        = ();
my @deletion_rows = ();
my @mutant_indel_rows = ();
my %mutant_indel_data = ();
my @indel_filter_rows = ();
my %seen_chr_indel          = ();
my %indel_ids         = ();


open( INDEL, "$input_indel" ) || modules::Exception->throw("Can't open file $input_indel\n");

#now open the filter vcf indel file 
while (<INDEL>) {
	#1	1301711	1301714	DEL^^^TGGG^^^Q113^^^REF_NORM123^^^SOMATIC^^^D56

       chomp;

       # Then, strip out chr from chromosome column at source, if it exists
       s/^chr//i;

	next unless /^[0-9XYM]/;
	my ($chr,$start,$end,$rest) = split("\t");

	if ( $OPT{chr} ) {
		next unless $OPT{chr} eq $chr;
	}
	if ( $OPT{coord} ) 	{
		next unless $OPT{coord} == $start;
	}
	
	$seen_chr_indel{$chr}++;
	
	my ( $type, $bases, $quality, $clr_str, $indel_class, $depth ) = split( '\^\^\^', $rest );
	if ($bases =~ /N/i) {
		next;
	}
	$quality =~ s/Q//;
	$depth =~ s/D//;
	
	my $zygosity = 'N/A'; 
	my $mpileup_str = 'N/A';
	my $quality_str = 'N/A';
	my $ref_freq = 0;
	my $var_freq = 0;
	my $var_count = 0;
	my $median = 0;
	my $total_reads = 0;

	#adjust coordinate for insertions
	my $new_start;
	if ($type eq 'DEL') {
		$new_start = $start;
	} else {
		$new_start = $start + 1;
	}

	if (exists $indel_info{$chr}{$new_start}) {
		$mpileup_str = $indel_info{$chr}{$new_start}{base_str};
		$quality_str = $indel_info{$chr}{$new_start}{qual_str};
		
		my $mpileup_pattern;
		my $length = length($bases);
		if ($type eq 'DEL') {
			$mpileup_pattern = '-'.$length.$bases;
		} else {
			$mpileup_pattern = '+'.$length.$bases;
		}
		
		my $pattern = quotemeta($mpileup_pattern);
		
		#Count up the indel occurrences
		while ($mpileup_str =~ /$pattern/gi) { 
			$var_count++
		}

		my $ref_count = 0;
		for my $char (split("",$mpileup_str)) {
			$ref_count++ if $char eq ',' || $char eq '.';
		}
		
		#refs are listed as .-2TT or ,-2TT so subtract these
		$ref_count = $ref_count - $var_count;
		
		$total_reads = $indel_info{$chr}{$new_start}{total_reads};
		$ref_freq = $ref_count == 0?0.0:sprintf( "%.3f",$ref_count/$total_reads);
		$var_freq = $var_count == 0?0.0:sprintf( "%.3f",$var_count/$total_reads);

		if ($var_freq > $homozygous_min_snv_frequency) {
			$zygosity = 'hom';
		} elsif ($var_freq > $heterozygous_min_snv_frequency) {
			$zygosity = 'het';
		}  else {
			$zygosity = 'other';
		}

		my @chars = split("",$indel_info{$chr}{$new_start}{qual_str});

	#    pop @chars;

		my @qualities;
		for ( my $i = 0 ; $i < scalar @chars ; $i++ ) {
			$qualities[$i] = ord( $chars[$i] ) - 33;
		}
		$median = modules::Utils->median( \@qualities );
		$quality_str = join( "\:", @qualities );
		
		
	} else {
		modules::Exception->warning("No pileup info for $chr $start $rest\n");
	}

	#Different depending on whether tumour or normal
	my $normal_base_str = "N/A";
	my $tumour_base_str;
	
	#Slightly different for normal and tumour samples
	if ($tumour) {
		#Get the normal base str
		$tumour_base_str = $mpileup_str;
		my $pileup_start = $start -1 ; 
		my @lines = split("\n",`grep -w $pileup_start $pileup_normal_file_indel`);
		
		if (@lines) {
			my @fields = split("\t",$lines[0]);
			$normal_base_str = $fields[4];
		} 
		
	} else {		
		$tumour_base_str = 'N/A';
		$normal_base_str = $mpileup_str;
	
	}
	
	
	my $indel_str = join( "\t",
								$chr,
								$start,
								$end,
								$rest . 
								$delim . $zygosity . 
								$delim . $var_count .'/' . $total_reads .
								$delim . 'MQ' . $median . 
								$delim . 'TBS:' . $tumour_base_str .
								$delim . 'NBS' . $normal_base_str .
								$delim . $quality_str
								) . "\n";
	
	print INDELOUT "$indel_str";
	
	my $ref_base = my $var_base;
	if ($type eq 'DEL') {
		$ref_base = $bases;
		$var_base = "-".$bases;
	} else {
		$ref_base = "N/A";
		$var_base = "+".$bases;
	}
	
	#default
	my $clr_score = 0;
	
	if ($clr_str =~ /(\d+)/) {
		$clr_score = $1;
	}	
			
	my $mutant_indel_key = $chr . ':' . $start;
	
	$normal_base_str =~ s/&/%/g; #Doesn't work in xml so hack as we don't use the ASCII quality
	$normal_base_str =~ s/</=/g;
	$normal_base_str =~ s/>/=/g;
	
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{ref_bases} = $ref_base;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{var_bases} = $var_base;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{ref_base_freq} = $ref_freq;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{var_base_freq} = $var_freq;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{median_quality_score} = $median;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_score} = $quality;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{tumour_base_string} = $tumour_base_str;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{clr_score} = $clr_score;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{variant_class} = $indel_class;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{normal_base_string} = $normal_base_str;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{read_depth} = $depth;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{var_depth} = $var_count;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{zygosity} = $zygosity;
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_type} = $type; 
	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{pass} = 1;
}


my $file_name_indel = $sample_name . '_' . $runid . '.db_variants.indel.xml';
$var_xml->create_var_xml(-file_name=>$file_name_indel,-data=>\%mutant_indel_data,-chr=>$chr_flag);								   
if ($chr_flag eq 'Y' || $chr_flag eq 'all') {
	my $full_xml_file = $run_dir . '/conf/'.$file_name_indel;
	$var_xml->split_xml(-file_name=>$full_xml_file);
}


close INDELOUT;
close $PF_INDEL;
