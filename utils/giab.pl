#! /usr/bin/perl -w

use strict;
use modules::Exception;
use modules::Vcf;
use modules::SystemCall;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vcf1=s",
	   "vcf2=s",
	   "vcf3=s",
	   "vcf_cutoff=i",
	   "variant_type=s",
	   "out=s",
	   "venn",
	   "tmpdir=s",
	   "sample_name1=s",
	   "sample_name2=s",
	   "sample_name3=s",
	   "genotype_qual=i"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf1} || !$OPT{vcf2});

=pod

=head1 SYNOPSIS

giab.pl -vcf1 vcf_file1 -vcf2 vcf_file2 -vcf3 vcf_file3  -vcf_cutoff vcf_cutoff_score(default=40) -genotype_qual GQ_Cutoff(ignores for hom calls) -sample_name1 sample_name_for_group_called_vcf1(same_for_2_and_3) -venn generate_venns -tmpdir tmp_dir_for_files(default=cwd/tmp) -variant_type report_only_this_variant_type(del, ins, or snv) -out outfile(default=summary.out)

Required flags: -vcf1 -vcf2

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

giab.pl -> compare giab vcfs

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

=cut

my $vcf1 = $OPT{vcf1};
my $vcf2 = $OPT{vcf2};
#Get all the giab info
my $svndir = $ENV{'SVNDIR'};

my $giab_base = "$svndir/conf/human/hs37d5/giab/2.19";

if (!-e $giab_base) {
	modules::Exception->throw("ERROR: dir $giab_base doesn't exist");
}
my $giab_region_file = $giab_base.'/giab.highqual.bed';
#Get the total number for FP rate
#my $num_snv = 2790044; #From wc $giab_base/*snv*
#my $num_del = 195463;
#my $num_ins = 177556;
my $num_snv = 26520; #From wc $giab_base/*snv*
my $num_del = 2016;
my $num_ins = 1972;
my $wd = `pwd`;
chomp $wd;
my $tmpdir = defined $OPT{tmpdir}?$OPT{tmpdir}:"$wd/tmp";
mkdir($tmpdir) if !-d $tmpdir; 

my $compare_number = defined $OPT{vcf3}?3:2;


if ( !-e $vcf1 || !-e $vcf2) {
	modules::Exception->throw("File $vcf1 or $vcf2 doesn't exist");	
}

my $geno_qual = defined $OPT{genotype_qual}?$OPT{genotype_qual}:0;
my $out = defined $OPT{out}?$OPT{out}:'summary.txt';
my $vcf_cutoff = defined $OPT{vcf_cutoff}?$OPT{vcf_cutoff}:40;
my $sample_name1 = defined $OPT{sample_name1}?$OPT{sample_name1}:0;
my $sample_name2 = defined $OPT{sample_name2}?$OPT{sample_name2}:0;
my $sample_name3 = defined $OPT{sample_name3}?$OPT{sample_name3}:0;


my $variant_type = 'all';
if (defined $OPT{variant_type}) {
	if ($OPT{variant_type} ne 'SNV' && $OPT{variant_type} ne 'DEL' && $OPT{variant_type} ne 'INS') {
		modules::Exception->throw("ERROR: variant_type must be INS, DEL, or SNV");
	}
	$variant_type = $OPT{variant_type};
} 



if (!-e "$tmpdir/$vcf1.snv") {
	my $vcf_data1;
	if ($sample_name1) {
		$vcf_data1 = parse_vcf(-vcf_file=>$vcf1,-sample_name=>$sample_name1);				
	} else {
		$vcf_data1 = parse_vcf(-vcf_file=>$vcf1);
	}
	open(SNV1,">$tmpdir/$vcf1.snv") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf1.snv\n");
	open(DEL1,">$tmpdir/$vcf1.del") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf1.del\n");
	open(INS1,">$tmpdir/$vcf1.ins") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf1.ins\n");

	#Write out text files for each variant type
	for my $var_key (sort {my ($a_chr,$a_start) = $a =~ /([0-9XYMT]+):(\d+)/; my ($b_chr,$b_start) = $b =~ /([0-9XYMT]+):(\d+)/; $a_chr cmp $b_chr || $a_start <=> $b_start} keys %{$vcf_data1}) {
	
		if ($OPT{variant_type}) {
			next unless $OPT{variant_type} eq $vcf_data1->{$var_key}{type};
		}
		my ($chr,$coord,$event) = split(':',$var_key);
		my ($start,$end) = split('-',$coord);
		my $line = join("\t", 
						$chr,
						$start,
						$end,
						$vcf_data1->{$var_key}{type}.';'.$var_key.';Q='.$vcf_data1->{$var_key}{qual}
						)."\n";
		if ($vcf_data1->{$var_key}{type} eq 'SNV') {
			print SNV1 $line;
		} elsif ($vcf_data1->{$var_key}{type} eq 'DEL') {
			print DEL1 $line;	
		} elsif ($vcf_data1->{$var_key}{type} eq 'INS') {
			print INS1 $line;
		}
	}
}


if (!-e "$tmpdir/$vcf2.snv") {
	my $vcf_data2;
	if ($sample_name2) {
		$vcf_data2 = parse_vcf(-vcf_file=>$vcf2,-sample_name=>$sample_name2);				
	} else {
		$vcf_data2 = parse_vcf(-vcf_file=>$vcf2);
	}
	open(SNV2,">$tmpdir/$vcf2.snv") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf2.snv\n");
	open(DEL2,">$tmpdir/$vcf2.del") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf2.del\n");
	open(INS2,">$tmpdir/$vcf2.ins") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf2.ins\n");



	for my $var_key (sort {my ($a_chr,$a_start) = $a =~ /([0-9XYMT]+):(\d+)/; my ($b_chr,$b_start) = $b =~ /([0-9XYMT]+):(\d+)/; $a_chr cmp $b_chr || $a_start <=> $b_start} keys %{$vcf_data2}) {
	
		if ($OPT{variant_type}) {
			next unless $OPT{variant_type} eq $vcf_data2->{$var_key}{type};
		}
		my ($chr,$coord,$event) = split(':',$var_key);
		my ($start,$end) = split('-',$coord);
		my $line = join("\t", 
						$chr,
						$start,
						$end,
						$vcf_data2->{$var_key}{type}.';'.$var_key.';Q='.$vcf_data2->{$var_key}{qual}
						)."\n";
		if ($vcf_data2->{$var_key}{type} eq 'SNV') {
			print SNV2 $line;
		} elsif ($vcf_data2->{$var_key}{type} eq 'DEL') {
			print DEL2 $line;	
		} elsif ($vcf_data2->{$var_key}{type} eq 'INS') {
			print INS2 $line;
		}
	}
}
my $vcf3 = $OPT{vcf3};

if ($compare_number == 3) {
	my $vcf_data3;
	if ($sample_name3) {
		$vcf_data3 = parse_vcf(-vcf_file=>$vcf3,-sample_name=>$sample_name3);				
	} else {
		$vcf_data3 = parse_vcf(-vcf_file=>$vcf3);
	}
	open(SNV3,">$tmpdir/$vcf3.snv") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf3.snv\n");
	open(DEL3,">$tmpdir/$vcf3.del") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf3.del\n");
	open(INS3,">$tmpdir/$vcf3.ins") || modules::Exception->throw("Can't open file to write $tmpdir/$vcf3.ins\n");

	for my $var_key (sort {my ($a_chr,$a_start) = $a =~ /([0-9XYMT]+):(\d+)/; my ($b_chr,$b_start) = $b =~ /([0-9XYMT]+):(\d+)/; $a_chr cmp $b_chr || $a_start <=> $b_start} keys %{$vcf_data3}) {
	
		if ($OPT{variant_type}) {
			next unless $OPT{variant_type} eq $vcf_data3->{$var_key}{type};
		}
		my ($chr,$coord,$event) = split(':',$var_key);
		my ($start,$end) = split('-',$coord);
		my $line = join("\t", 
						$chr,
						$start,
						$end,
						$vcf_data3->{$var_key}{type}.';'.$var_key.';Q='.$vcf_data3->{$var_key}{qual}
						)."\n";
		if ($vcf_data3->{$var_key}{type} eq 'SNV') {
			print SNV3 $line;
		} elsif ($vcf_data3->{$var_key}{type} eq 'DEL') {
			print DEL3 $line;	
		} elsif ($vcf_data3->{$var_key}{type} eq 'INS') {
			print INS3 $line;
		}
	}
}

my @var_types = qw(snv del ins);

#now build up the overlap commnads
my $overlap_bin = "$svndir/utils/overlap_files.pl";
my @commands = ();
my @first_commands = ();
for my $var_type ( @var_types ) {
	push @first_commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf1.$var_type -coord $giab_region_file > $tmpdir/$vcf1.$var_type.regions";
	push @first_commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf2.$var_type -coord $giab_region_file > $tmpdir/$vcf2.$var_type.regions";
}

if ($compare_number == 3) {
	for my $var_type ( @var_types ) {
		push @first_commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf3.$var_type -coord $giab_region_file > $tmpdir/$vcf3.$var_type.regions";
	}
}

my $sys_call = modules::SystemCall->new();
my %counts = ();

for my $command (@first_commands) {
    my (undef,$output) = split("> ",$command);
    $sys_call->run($command) unless -e $output;
    my ($vcf,$var_type) = $output =~ /\/(.+)\.(\w{3})\.regions/;
    if ($output !~ /[0-9X]$/) {
    	$counts{$var_type}{$vcf}{'total'}{'regions'} = count_lines($output);
    }
}

my @chrs = (1..22);
push @chrs, 'X';

for my $chr (@chrs) {
	for my $var_type (@var_types) {
		my $out1 = `grep -w ^$chr $tmpdir/$vcf1.$var_type.regions 2>/dev/null | head`;
	    my $out2 = `grep -w ^$chr $tmpdir/$vcf2.$var_type.regions 2>/dev/null | head`;
		if ($out1 ne "") {
			push @commands, "grep -w ^$chr $tmpdir/$vcf1.$var_type.regions > $tmpdir/$vcf1.$var_type.regions.$chr";
			push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf1.$var_type.regions.$chr -coord $giab_base/giab_${var_type}.$chr > $tmpdir/$vcf1.$var_type.regions.overlap.$chr";
		}
	
		if ($out2 ne "") {
	    	push @commands, "grep -w ^$chr $tmpdir/$vcf2.$var_type.regions > $tmpdir/$vcf2.$var_type.regions.$chr";
	        push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf2.$var_type.regions.$chr -coord $giab_base/giab_${var_type}.$chr > $tmpdir/$vcf2.$var_type.regions.overlap.$chr";
	 	}
	  	
	}
}

for my $var_type (@var_types) {
	push @commands, "cat $tmpdir/$vcf1*$var_type*overlap* >> $tmpdir/$vcf1.$var_type.regions.overlap";
	push @commands, "cat $tmpdir/$vcf2*$var_type*overlap* >> $tmpdir/$vcf2.$var_type.regions.overlap";
}



if ($compare_number == 3) {
	for my $chr (@chrs) {
		for my $var_type (@var_types) {
			my $out3 = `grep -w ^$chr $tmpdir/$vcf3.$var_type.regions 2>/dev/null | head`;
			if ($out3 ne "") {
		    	push @commands, "grep -w ^$chr $tmpdir/$vcf3.$var_type.regions > $tmpdir/$vcf3.$var_type.regions.$chr";
		        push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf3.$var_type.regions.$chr -coord $giab_base/giab_${var_type}.$chr > $tmpdir/$vcf3.$var_type.regions.overlap.$chr";
		 	}
		}
	}
	for my $var_type (@var_types) {
		push @commands, "cat $tmpdir/$vcf3*$var_type*overlap* >> $tmpdir/$vcf3.$var_type.regions.overlap";
	}
}

#Get tool-specific numbers
if ($compare_number == 2) {
	#Easy case here as we just use the two files to get the non-matches
	for my $var_type (@var_types) {
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf1.$var_type.regions -coord $tmpdir/$vcf2.$var_type.regions -fail > $tmpdir/$vcf1.$var_type.regions.tool";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf2.$var_type.regions -coord $tmpdir/$vcf1.$var_type.regions -fail > $tmpdir/$vcf2.$var_type.regions.tool";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf1.$var_type.regions.tool -coord $giab_base/giab_${var_type}.all > $tmpdir/$vcf1.$var_type.regions.tool.overlap";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf2.$var_type.regions.tool -coord $giab_base/giab_${var_type}.all > $tmpdir/$vcf2.$var_type.regions.tool.overlap";
	}
} else {
	#Here we need to concat the other files to get the unique variants
	for my $var_type (@var_types) {
		push @commands, "cat $tmpdir/$vcf1.$var_type.regions $tmpdir/$vcf2.$var_type.regions >> $tmpdir/$vcf1.$vcf2.$var_type.tmp";
		push @commands, "cat $tmpdir/$vcf1.$var_type.regions $tmpdir/$vcf3.$var_type.regions >> $tmpdir/$vcf1.$vcf3.$var_type.tmp";
		push @commands, "cat $tmpdir/$vcf2.$var_type.regions $tmpdir/$vcf3.$var_type.regions >> $tmpdir/$vcf2.$vcf3.$var_type.tmp";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf1.$var_type.regions -coord $tmpdir/$vcf2.$vcf3.$var_type.tmp -fail > $tmpdir/$vcf1.$var_type.regions.tool";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf2.$var_type.regions -coord $tmpdir/$vcf1.$vcf3.$var_type.tmp -fail > $tmpdir/$vcf2.$var_type.regions.tool";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf3.$var_type.regions -coord $tmpdir/$vcf1.$vcf2.$var_type.tmp -fail > $tmpdir/$vcf3.$var_type.regions.tool";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf1.$var_type.regions.tool -coord $giab_base/giab_${var_type}.all > $tmpdir/$vcf1.$var_type.regions.tool.overlap";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf2.$var_type.regions.tool -coord $giab_base/giab_${var_type}.all > $tmpdir/$vcf2.$var_type.regions.tool.overlap";
		push @commands, "$overlap_bin -just_overlap -ref $tmpdir/$vcf3.$var_type.regions.tool -coord $giab_base/giab_${var_type}.all > $tmpdir/$vcf3.$var_type.regions.tool.overlap";
	}
}



for my $command (@commands) {
  	#Now get the counts
  	my (undef,$output) = split("> ",$command);
	$sys_call->run($command) unless -e $output;
	my $total = $output =~ /\.tool/?'tool':'total';
	my $overlap = $output =~ /overlap$/?'overlap':'regions';
	my ($vcf,$var_type) = $output =~ /\/(.+)\.(\w{3})\.regions/;
	if ($output !~ /[0-9X]$/ && $output !~ /tmp$/) {
		$counts{$var_type}{$vcf}{$total}{$overlap} = count_lines($output);
	}
}

print Dumper \%counts;
open(OUT,">$out") || modules::Exception->throw("Can't open output file $out\n");
print join("\t",'VAR_TYPE','FILE','TOTAL_VAR_COUNT','TOTAL_FP','TOTAL_FN','TOOL_SPECIFIC_CALLS','TOOL_SPECIFIC_FP')."\n\n";
print OUT join("\t",'VAR_TYPE','FILE','TOTAL_VAR_COUNT','TOTAL_FP','TOTAL_FN','TOOL_SPECIFIC_CALLS','TOOL_SPECIFIC_FP')."\n\n";

for my $var_type (@var_types) {
	for my $vcf (sort keys %{$counts{$var_type}}) {
		my $region = $counts{$var_type}{$vcf}{'total'}{'regions'};
		my $overlap =  $counts{$var_type}{$vcf}{'total'}{'overlap'};
		my $tp = $overlap/$region*100;
		my $fp = sprintf("%.2f",100-$tp);  
		my $fn;
		
		if ($var_type eq 'snv') {
			my $tn = $overlap/$num_snv*100;
			$fn = sprintf("%.2f",100-$tn); 
		} elsif ($var_type eq 'del') {
			my $tn = $overlap/$num_del*100;
			$fn = sprintf("%.2f",100-$tn); 
		} elsif ($var_type eq 'ins') {
			my $tn = $overlap/$num_ins*100;
			$fn = sprintf("%.2f",100-$tn); 
		}
		
		my $tool_region = $counts{$var_type}{$vcf}{'tool'}{'regions'};
		my $tool_overlap =  $counts{$var_type}{$vcf}{'tool'}{'overlap'};
		my $tool_tp = $tool_overlap/$tool_region*100;
		my $tool_fp = sprintf("%.2f",100-$tool_tp);  
		
		print join("\t",$var_type,$vcf,$region,$fp,$fn,$tool_region,$tool_fp)."\n";
		print OUT join("\t",$var_type,$vcf,$region,$fp,$fn,$tool_region,$tool_fp)."\n";
	}
	print OUT "\n";
	print "\n";
}


if ($OPT{venn}) {
	my @venn_commands = ();
	my $venn_bin = "$svndir/utils/venn.pl";
	for my $var_type ( @var_types ) {
	    (my $s1 = $vcf1) =~ s/.vcf//;
	    (my $s2 = $vcf2) =~ s/.vcf//;
	    
		if ($compare_number == 2) {
			my $files = join(",",$tmpdir.'/'.$vcf1.'.'.$var_type,$tmpdir.'/'.$vcf2.'.'.$var_type);
			my $title = $s1 . '_' . $s2 . '_' . $var_type;
			push @venn_commands, join(" ",$venn_bin,'-files',$files,'-s1',$s1,'-s2',$s2,'-title',$title);
		} else {
			(my $s3 = $vcf3) =~ s/.vcf//;
			my $files = join(",",$tmpdir.'/'.$vcf1.'.'.$var_type,$tmpdir.'/'.$vcf2.'.'.$var_type,$tmpdir.'/'.$vcf3.'.'.$var_type);
			my $title = $s1 . '_' . $s2 . '_' . $s3 . '_'. $var_type;
			push @venn_commands, join(" ",$venn_bin,'-files',$files,'-s1',$s1,'-s2',$s2,'-s3',$s3,'-title',$title);
		}
	}
	for my $command (@venn_commands) {
		$sys_call->run($command);
	}
	
	#./venn.pl -files all_chr22//anu_samtools_chr22.vcf.snv,all_chr22//anu_gatk_chr22.vcf.snv,all_chr22//garvan_gatk_chr22.vcf.snv -s1 anu_samtools_chr22 -s2 anu_gatk_chr22 -s3 garvan_gatk_chr22 -title anu_samtools_chr22_anu_gatk_chr22_garvan_gatk_chr22_venn_SNV
	
}


sub count_lines {
	my ($file) = @_;
	if (!-e $file) {
		modules::Exception->throw("Can't open file file\n");
	}
	my ($out) = split(' ',`wc -l $file`);
	chomp $out;
	return $out;
}



sub parse_vcf {
    
   	my (@args) = @_;

    my @required_args = (
			             -vcf_file
						 );

	my %args = @args;

    foreach my $required_arg (@required_args) {

		if (! defined $args{$required_arg}){
		    modules::Exception->throw("Required argument [$required_arg] not set");
		}
    }
    
    if ( !-e $args{-vcf_file} ) {
    	modules::Exception->throw("File $args{-vcf_file} doesn't exist");	
    }
    
    
    my $vcf = $args{-vcf_file};
    my $sample_name = defined $args{-sample_name}?$args{-sample_name}:0;   
    my %vcf_data = ();
    
    open(VCF,$vcf) || modules::Exception->throw("Can't open file $vcf\n");
    
   	#Genotype quality index
   	my $gq_index = 0;
   	my $sample_index = 0;

    while (<VCF>) {
    	if ($sample_name) {
			unless (/#CHROM/) {
	    		next if /^#/;				
			}    		
    	} else {
	    	next if /^#/;
    	}
    	
    	next unless /\w/;
    	chomp;
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,@alleles) = split;
    	my $line = $_;
   		my @fields = split;
    	
    	if (/CHROM/ && $sample_name) {
    		#Get the sample_index
    		for ( my $count = 0 ; $count < @fields ; $count++ ) {
    		    if ($fields[$count] eq $sample_name) {
    		    	$sample_index = $count;
    		    }
    		}
    		
    		if (!$sample_index) {
    			modules::Exception->throw("ERROR: Can't find sample $sample_name in line $_");
    		}
    		next;
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



		my @vars;
		if ($sample_name) {
			#Get only the relevant genotype
	    	my $var_index = 0;
			my ($zyg_str) = split(':',$fields[$sample_index]);
    		next if $zyg_str eq '0/0';
    		next if $zyg_str eq './.';
    		my @nums = split('/',$zyg_str);
    		#Get the index from the gt string (eg 0/1 or 0/2 etc)
    		for my $num (@nums) {
    			if ($num != 0) {
    				$var_index = $num;
    			}
    		}
    		#print "sample index $sample_index $zyg_str VI $var_index\n";
    		
    		if (!$var_index) {
    			modules::Exception->throw("Can't get var_index for $_");
    		}
    		
    		#Only add the relevant genotype
    		my @tmp_vars = split(",",$var_str);
    		
    		push @vars,$tmp_vars[$var_index-1];
    		
    		if ($geno_qual && !$gq_index) {
    			#only do this once with the first non-header line; get the index of the GQ field
    			my @gt_fields = split(':',$gt_fields);
    			for ( my $count = 0 ; $count < @gt_fields ; $count++ ) {
    			    if ($gt_fields[$count] eq 'GQ') {
    			    	$gq_index = $count;
    			    }
    			}
    			if (!$gq_index) {
    				modules::Exception->throw("ERROR: Can't get GQ index for $gt_fields");
    			}
    		}
		} else {
			@vars = split(",",$var_str);
		}

		for my $var ( @vars ) {
			next if $var eq '*'; #Due to upstream deletion
			my ($var_key,$var_type) = _get_variant_key(-type=>'vcf',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var);


			my ($start,$end) = $var_key =~ /(\d+)-(\d+)/;

			

			#If fails quality test; only apply if quality score available
			if ($qual ne '.' && $qual <= $vcf_cutoff) {
				next;
			}

			if ($sample_name && $geno_qual) {
				#Here we additionally screen on GQ
				my @gt_fields = split(':',$fields[$sample_index]);
				my $gq = $gt_fields[$gq_index];
				
				if ($gq < $geno_qual) {
					#Make sure it's not a het or hom decision -> this is still a variant
					my @nums = split('/',$gt_fields[0]);
					if ($nums[0] eq $nums[1]) {
						next;
					}
				}
				
			}

			
			if ($qual eq '.') {
				$vcf_data{$var_key}{qual} = "N/A";				
			} else {
				$vcf_data{$var_key}{qual} = $qual;
			}
			$vcf_data{$var_key}{type} = $var_type;
			
			if ($OPT{keep_zyg}) {
				if ($sample_name) {
					$vcf_data{$var_key}{zyg} = \split(':',$fields[$sample_index]);
				} else {
					$vcf_data{$var_key}{zyg} = \@alleles;
				}
			}
		}
    }
    return \%vcf_data;
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
    					-type
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
    my $type = $args{-type};
    
    my $start_coord = my $end_coord = my $bases;
    my $length_ref = length($ref);
    my $length_var = length($var);
    my $var_type;
    
    if ($type eq 'vcf') {
		if ($length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord + 1;
			$end_coord = $start_coord + $length_ref - $length_var - 1;				
			my $del_length = $length_ref - $length_var;
			#print "VCF R $ref L $del_length\n";
			
			$bases = '-'. substr($ref,1,$del_length);
		} elsif ($length_ref < $length_var) {
			#Add the ref length and var length difference to the coord 
			#$start_coord = $end_coord = $first_coord + 1;
			$var_type = 'INS';
			$start_coord = $end_coord = $first_coord;
			my $ins_length = $length_var - $length_ref;
			$bases = '+'.substr($var,1,$ins_length);
		} elsif ($length_ref == $length_var && $length_ref != 1) {
			#Handling for cases like AT->GC ot ATATA->CTATA
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			my ($first_ref) = split("",$ref);
			my ($first_var) = split("",$var);
			if ($first_var eq $first_ref) {
				modules::Exception->throw("ERROR: No handling for event $chr . ':'.$start_coord .'-'.$end_coord .':' .$ref -> $var ");
			}
			$bases = $first_ref . '->' .$first_var;
		} else {
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			$bases = $ref . '->' .$var;
		}
    	
    } else {
    	if ($ref eq '-' || $length_ref < $length_var) {
	    	$start_coord = $end_coord = $first_coord - 1;
			$var_type = 'INS';	
			my $ins_length = $length_var - $length_ref;
			$ins_length++ if $ref eq '-';
			$bases = '+'.substr($var,0,$ins_length);
		}  elsif ($var eq '-' || $length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord;
			my $del_length = $length_ref - $length_var;
			$end_coord = $start_coord + $length_ref - $length_var - 1;
			$del_length++ if $var eq '-';
			$end_coord++ if $var eq '-';
			$bases = '-'. substr($ref,0,$del_length);	
		} elsif ($length_ref == $length_var && $length_ref == 1) {
			#single snvs
			$var_type = 'SNV'; 
			$bases = $ref .'->'.$var;
			$start_coord = $end_coord = $first_coord;
		} else {
			modules::Exception->warning("ERROR: Can't identify var_type doesn't match any var type\n");
			#next;
		}
    }
    
	my $var_key = $chr . ':'.$start_coord .'-'.$end_coord .':' .$bases;
	return($var_key,$var_type);
	
}



=comment

sub parse_vcf {
    my ($file) = @_;
    my %vcf_data = ();
    
    open(VCF,$file) || modules::Exception->throw("Can't open file file\n");
    
    while (<VCF>) {
    	next if /^#/;
    	next unless /\w/;
    	chomp;
    	
    	
    	my ($chr,$first_coord,undef,$ref,$var_str,$qual,undef,$rest,$gt_fields,@ped_alleles) = split;
    	$chr =~ s/chr//;	
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
		
		for my $var ( @vars ) {
			my ($var_key,$var_type) = _get_variant_key(-type=>'vcf',-chrom=>$chr,-first=>$first_coord,-ref_seq=>$ref,-var_seq=>$var);


			my ($start,$end) = $var_key =~ /(\d+)-(\d+)/;


			#If fails quality test; only apply if quality score available
			if ($qual ne '.' && $qual <= $vcf_cutoff) {
				next;
			}


			my $length_ref = length($ref);
			my $length_var = length($var);
			
			if ($length_ref == $length_var && $length_ref != 1) {
				#modules::Exception->warning("Skip vcf entry $ref -> $var");
			}
			
			if ($qual eq '.') {
				$vcf_data{$var_key}{qual} = "N/A";				
			} else {
				$vcf_data{$var_key}{qual} = $qual;
			}
			$vcf_data{$var_key}{type} = $var_type;
		}
    }
    return \%vcf_data;
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
    					-type
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
    my $type = $args{-type};
    
    my $start_coord = my $end_coord = my $bases;
    my $length_ref = length($ref);
    my $length_var = length($var);
    my $var_type;
    
    if ($type eq 'vcf') {
		if ($length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord + 1;
			$end_coord = $start_coord + $length_ref - $length_var - 1;				
			my $del_length = $length_ref - $length_var;
			#print "VCF R $ref L $del_length\n";
			
			$bases = '-'. substr($ref,1,$del_length);
		} elsif ($length_ref < $length_var) {
			#Add the ref length and var length difference to the coord 
			#$start_coord = $end_coord = $first_coord + 1;
			$var_type = 'INS';
			$start_coord = $end_coord = $first_coord;
			my $ins_length = $length_var - $length_ref;
			$bases = '+'.substr($var,1,$ins_length);
		} else {
			$var_type = 'SNV';
			$start_coord = $end_coord = $first_coord;
			$bases = $ref . '->' .$var;
		}
    	
    } else {
    	if ($ref eq '-' || $length_ref < $length_var) {
	    	$start_coord = $end_coord = $first_coord - 1;
			$var_type = 'INS';	
			my $ins_length = $length_var - $length_ref;
			$ins_length++ if $ref eq '-';
			$bases = '+'.substr($var,0,$ins_length);
		}  elsif ($var eq '-' || $length_ref > $length_var) {
			$var_type = 'DEL';
			$start_coord = $first_coord;
			my $del_length = $length_ref - $length_var;
			$end_coord = $start_coord + $length_ref - $length_var - 1;
			$del_length++ if $var eq '-';
			$end_coord++ if $var eq '-';
			$bases = '-'. substr($ref,0,$del_length);	
		} elsif ($length_ref == $length_var && $length_ref == 1) {
			#single snvs
			$var_type = 'SNV'; 
			$bases = $ref .'->'.$var;
			$start_coord = $end_coord = $first_coord;
		} else {
			modules::Exception->warning("ERROR: Can't identify var_type doesn't match any var type\n");
			#next;
		}
    }
    
	my $var_key = $chr . ':'.$start_coord .'-'.$end_coord .':' .$bases;
	return($var_key,$var_type);
	
}
=cut