package modules::Utils;

use strict;
use modules::Exception;
use Data::Dumper;

sub median() {
    my ($self, $ra_values) = @_;

    my @sorted = sort {$a <=> $b} @$ra_values;

    my $length = scalar @sorted;

    if ($length%2 == 0 && $length >2){
    	return ($sorted[($length/2 - 1)] + $sorted[($length/2)])/2;
    } elsif ($length == 2) {
    	return ($sorted[0] + $sorted[1])/2;
    } else {
    	return $sorted[$length/2];
    }
}

sub quartiles() {
    my ($self, $ra_values) = @_;

    my @sorted = sort {$a <=> $b} @$ra_values;

    my $length = scalar @sorted;

    my $second_quartile = $self->median($ra_values);

    my $first_quartile;
    my $third_quartile;

    # Calculating quartiles excluding the value(s) used to calculate the overall median
    if ($length%2 == 0 && $length > 4){
		my @subset = @sorted[0..(($length/2)-2)];
		$first_quartile = $self->median(\@subset);
		@subset = @sorted[($length/2)..$length-1];
		$third_quartile = $self->median(\@subset);
    } elsif ($length <= 4) {
	  	modules::Exception->warning("Choosing not to calculate quartiles on four or fewer values");
	  	return [undef, undef, undef, undef, undef];
    } else {
		my @subset = @sorted[0..int($length/2)-1];
		$first_quartile = $self->median(\@subset);
		@subset = @sorted[(int($length/2)+1)..$length-1];
		$third_quartile = $self->median(\@subset);
    }

    my $min_value = $sorted[0];
    my $max_value = $sorted[-1];

    return [$min_value, 
	    $first_quartile, 
	    $second_quartile, 
	    $third_quartile, 
	    $max_value];
}


sub min() {
	my ($self, $ra_values) = @_;
	my @sorted = sort {$a <=> $b} @$ra_values;
	return $sorted[0];
}

sub max() {
	my ($self, $ra_values) = @_;
	my @sorted = sort {$a <=> $b} @$ra_values;
	return $sorted[-1];
}



sub GetTime {
	my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) = localtime time;
	my $timestr = sprintf( "%02d-%02d-%d_%02d:%02d:%02d",
						   $mday, 
						   $mon + 1, 
						   $year + 1900,
						   $hour, 
						   $min,
						   $sec );
	return $timestr;
	
}

sub translate {
	my ($self, $rna) = @_;
	my %genetic_code = (
			'TCA' => 'S', # Serine
			'TCC' => 'S', # Serine
			'TCG' => 'S', # Serine
			'TCT' => 'S', # Serine
			'TTC' => 'F', # Phenylalanine
			'TTT' => 'F', # Phenylalanine
			'TTA' => 'L', # Leucine
			'TTG' => 'L', # Leucine
			'TAC' => 'Y', # Tyrosine
			'TAT' => 'Y', # Tyrosine
			'TAA' => '_', # Stop
			'TAG' => '_', # Stop
			'TGC' => 'C', # Cysteine
			'TGT' => 'C', # Cysteine
			'TGA' => '_', # Stop
			'TGG' => 'W', # Tryptophan
			'CTA' => 'L', # Leucine
			'CTC' => 'L', # Leucine
			'CTG' => 'L', # Leucine
			'CTT' => 'L', # Leucine
			'CCA' => 'P', # Proline
			'CAT' => 'H', # Histidine
			'CAA' => 'Q', # Glutamine
			'CAG' => 'Q', # Glutamine
			'CGA' => 'R', # Arginine
			'CGC' => 'R', # Arginine
			'CGG' => 'R', # Arginine
			'CGT' => 'R', # Arginine
			'ATA' => 'I', # Isoleucine
			'ATC' => 'I', # Isoleucine
			'ATT' => 'I', # Isoleucine
			'ATG' => 'M', # Methionine
			'ACA' => 'T', # Threonine
			'ACC' => 'T', # Threonine
			'ACG' => 'T', # Threonine
			'ACT' => 'T', # Threonine
			'AAC' => 'N', # Asparagine
			'AAT' => 'N', # Asparagine
			'AAA' => 'K', # Lysine
			'AAG' => 'K', # Lysine
			'AGC' => 'S', # Serine
			'AGT' => 'S', # Serine
			'AGA' => 'R', # Arginine
			'AGG' => 'R', # Arginine
			'CCC' => 'P', # Proline
			'CCG' => 'P', # Proline
			'CCT' => 'P', # Proline
			'CAC' => 'H', # Histidine
			'GTA' => 'V', # Valine
			'GTC' => 'V', # Valine
			'GTG' => 'V', # Valine
			'GTT' => 'V', # Valine
			'GCA' => 'A', # Alanine
			'GCC' => 'A', # Alanine
			'GCG' => 'A', # Alanine
			'GCT' => 'A', # Alanine
			'GAC' => 'D', # Aspartic Acid
			'GAT' => 'D', # Aspartic Acid
			'GAA' => 'E', # Glutamic Acid
			'GAG' => 'E', # Glutamic Acid
			'GGA' => 'G', # Glycine
			'GGC' => 'G', # Glycine
			'GGG' => 'G', # Glycine
			'GGT' => 'G'  # Glycine
			);
	my ($protein) = '';
	for(my $i=0;$i<length($rna)-2;$i+=3)
	{
		my $codon = substr($rna,$i,3);
		$protein .= $genetic_code{$codon};
	}
	return $protein;
}

sub revcom {
	my ( $self, $seq ) = @_;
    my $revcomp = reverse($seq);
  	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
  	return $revcomp;
}

sub pileup_string {
	my ( $self, $pileup_str ) = @_;
	my @bases = split("",uc($pileup_str));
    my $count = 0;
    my $allele_count = 0;
    my $base_count = @bases;
    my %final_bases = ();
    
    while ($count < $base_count) {
     	if (defined $bases[$count+2] && ($bases[$count+1] eq '+' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
    		my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		my $indel_end = $length_indel+1+length($length_indel)+$count;
    		my $indel_bases = join("",@bases[$count+1..$indel_end]);
    		$count +=  $length_indel+2+length($length_indel);
    		$final_bases{$indel_bases}++;
    		$allele_count++;
    		next;
    	} elsif (defined $bases[$count+2] && ($bases[$count+1] eq '-' && $count+2 != $base_count && $bases[$count+2] =~ /(\d)/)) {
    		my $rest_indel = join("",@bases[$count+2..$#bases]);
			my ($length_indel) = $rest_indel =~ /(\d+)/;
    		my $indel_end = $length_indel+1+length($length_indel)+$count;
	    	my $indel_bases = join("",@bases[$count+1..$indel_end]);
    		$count +=  $length_indel+2+length($length_indel);
    		$final_bases{$indel_bases}++;
    		$allele_count++;
    		next;
    		
   		} elsif ($bases[$count] =~ /[ATCGN]/) {
	    	$final_bases{$bases[$count]}++;
	    	$allele_count++;
	    } elsif ($bases[$count] eq '.' || $bases[$count] eq ',') {
	    	$final_bases{'ref'}++;
	    	$allele_count++;
	    } elsif ($bases[$count] eq '*') {
	    	$final_bases{'other_del'}++;
	    	$allele_count++;
	    }
	    $count++;
	}
    
    my $zyg = 'het'; #default is het
    my @pileup_strings;
    for my $allele (sort {$final_bases{$b}<=>$final_bases{$a}} keys %final_bases) {
    	#Simply hom test (90% non reference allele) and REF test
    	if ($final_bases{$allele}/$allele_count >= 0.9 && $allele ne 'ref') {
    		$zyg = 'hom';
    	} elsif ($final_bases{$allele}/$allele_count >= 0.9 && $allele eq 'ref') {
    		$zyg = 'ref';
    	}
    	
    	push @pileup_strings, $allele.':'.$final_bases{$allele};
    }
    my $pileup_string = join('_',@pileup_strings);
    
    return ($pileup_string,$zyg);    
	
	
	
}

return 1;
