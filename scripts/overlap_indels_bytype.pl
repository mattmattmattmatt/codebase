#! /usr/bin/perl -w

use strict;
use modules::Adaptors::Variant;
use modules::Adaptors::Filter;
use modules::Adaptors::Variant_Filter;
use modules::Overlap;
use modules::Exception;
use modules::Pipeline;
use modules::Adaptors::Source_Group;
use modules::Adaptors::Sample;
use modules::Adaptors::BulkInsert;
use modules::VariantXML;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use vars qw(%OPT);

GetOptions(\%OPT, 
		   "help|h",
		   "man|m",
		   "writeDB=i",
		   "coord_file=s",
		   "coord_filtername=s",
		   "overlap_outfile=s",
		   "ref_file=s",
		   "indel_type=s",
		   "runid=i",
		   "chr=s",
		   "match=s",
		   "ref=s"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{coord_file} || !$OPT{overlap_outfile} || !$OPT{ref_file});

	   
=pod

=head1 SYNOPSIS

overlap_indels_bytype.pl -match determine_if_match_means_pass_or_fail(options=PASS/FAIL;default=FAIL) -indel_type comma_delim_list_of_Var_types(default=INS,DEL)-runid runid -ref_filtername db_filtername -overlap_outfile output_overlap_file -ref_file input_overlap_file -coord_file formatted_snp_file -chr chrname -writeDB write_results_to_DB(default=0) [options]

Required flags: -coord_filtername -coord_file -runid -overlap_outfile

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

overlap_indels_bytype.pl -> Script to run a overlaps for indels for dbsnp

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE



=cut

my $coord_filter = defined $OPT{coord_filtername}?$OPT{coord_filtername}:'filter_dbsnp_indel';

#Dbsnp snv flag requires us to insert the allele freq into the database if the info is there
my $dbsnp_indel_flag = 0;
if ($coord_filter =~ /filter_dbsnp_indel/) {
	$dbsnp_indel_flag = 1;
}

my $filter_common = 0;
my $runid = $OPT{runid};
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');


my $ignore_string;
if ($coord_filter =~ /filter_common/) {
	#here we get the bioreps from the database to pass as -ignore arguments to overlap
	if (!defined $OPT{sample_name}) {
		modules::Exception->throw("ERROR: Require args -sample_name with filter_common");
	}
	my $source_group_name = modules::Pipeline::get_source_group_name($sample_name);
	my ($source_group_obj) = modules::Adaptors::Source_Group->search(source_group_name=>$source_group_name);
	
	if (!defined $source_group_obj) {
		modules::Exception->throw("ERROR: Can't get source group object from name $source_group_name");
	}
	my @sample_objs = modules::Adaptors::Sample->search(source_group_id=>$source_group_obj->id);
	if (!@sample_objs) {
		modules::Exception->throw("ERROR: Can't get sample objects from $source_group_name");
	}
	
	my @sample_names = ();
	for my $sample_obj ( @sample_objs ) {
	    push @sample_names, $sample_obj->sample_name;
	}
	$ignore_string = join(",",@sample_names);
	$filter_common = 1;
}
my $filter_exac = 0;
if ($coord_filter =~ /filter_exac/) {
	$filter_exac = 1;
}

my $filter_gnomad = 0;
if ($coord_filter =~ /filter_gnomad/) {
	$filter_gnomad = 1;
}


my $filter_clinvar = 0;
if ($coord_filter =~ /filter_clinvar/) {
        $filter_clinvar = 1;
}

my $coord_file = $OPT{coord_file};
my $chromosome = defined $OPT{chr}?$OPT{chr}:'';

my $last_chr = 0;
if (defined $OPT{chr} && $OPT{chr} eq 'Y') {
	$last_chr = 1;
}

my $match = defined $OPT{match}?$OPT{match}:'FAIL';
my $filterpass;
if ($match ne 'PASS' && $match ne 'FAIL') {
	modules::Exception->throw("Match parameter must be set to PASS or FAIL (default=FAIL)");
} else {
	$filterpass = $match eq 'FAIL'?0:1;
}

my $write_to_db = defined $OPT{writeDB}?$OPT{writeDB}:0;

my $indel_type = defined $OPT{indel_type}?$OPT{indel_type}:'INS,DEL';

my $overlap_args;

my $overlap = modules::Overlap->new();

#Need to create tmp files for overlap mo
my $outfile_match = $OPT{overlap_outfile};
#Now call for nonmatching; need to change some of the arguments
(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;

my $ref_file =  $OPT{ref_file};

#print "REF $ref_file OUT $outfile_match\n";
my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
#Build up the arguments for the overlap call
my %args = (
			-ref=>$ref_file,
			-coord=>$coord_file,
			-silent=>1, #Don't need the file
			-exact=>1 #Only allow exact coordinate matches
			);

#Here we have to get bioreps from the database and pass them to overlap object
if ($filter_common && length($ignore_string)) {
	$args{-ignore} = $ignore_string;
}

if ($OPT{chr}) {
	$args{-chr} = $chromosome;
	#If step died in the middle remove the appended file
	if ($OPT{chr} eq 1) {
		
		system("rm $outfile_match") if -e $outfile_match;
		system("rm $outfile_nomatch") if -e $outfile_nomatch;
	}
} 



#Call for matching
($overlap_match,$match_count) = $overlap->overlap(%args);

$args{-fail} = 1;
($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);



# Fetch INDELS from database

my %search_params = ('run_id' => $runid);

if ($OPT{chr}){
	$search_params{'chr'} = $chromosome;	
}

my @db_indels;
my @indel_types = split(",",$indel_type);
my @sorted_indels;

my $indels;
if ($OPT{chr}) {
	my $file_name = $sample_name . '_' . $runid . '.db_variants.indel.'.$OPT{chr}.'.xml';
	my $full_file = $var_xml->outdir().'/'.$file_name;
	if (-e $full_file) {
		$var_xml->load_xml_keys(-file_name=>$file_name);
		$indels = $var_xml->search_xml_keys(-chr=>$OPT{chr},-var_type=>'indel');
	}
} else {
	my $file_name = $sample_name . '_' . $runid . '.db_variants.indel.xml';
	$var_xml->load_xml_keys(-file_name=>$file_name);
	$indels = $var_xml->search_xml_keys(-var_type=>'indel');
}
if ($indels && @{$indels}) {
	@sorted_indels = sort {$a->{start_coord} <=> $b->{start_coord}} @{$indels};
}


# Sort INDELs

my @indel_filter_inserts;

#Here we overwrite the default match/nomatch files as coordinate matching isn't enough with indel
open(MATCH,">>$outfile_match") || modules::Exception->throw("Can't open file to write $outfile_match\n");
open(NOMATCH,">>$outfile_nomatch") || modules::Exception->throw("Can't open file to write $outfile_nomatch\n");

my %mutant_indel_data;

foreach my $newindel  (@sorted_indels) {
	my %indel_filter;
	
		
	my $newindel_chr = $newindel->{chr};
	my $newindel_start = $newindel->{start_coord};
	my $newindel_end = $newindel->{end_coord};
	my $newindel_varbases = $newindel->{var_bases};
	$newindel_varbases =~ s/[+-]//;
		
	my $matched = 0;	
		
	if (exists $overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}) {
		
		#If the indel is potentially matched set as non matched by default and change once we find the exact match
		$indel_filter{'filtermatch'} = 0;
		$indel_filter{'filterpass'} = $filterpass == 0?1:0;
		
		if ($dbsnp_indel_flag) {
			
			#Check if allele freq is available for snvs
		    for my $indel_entry ( sort keys %{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}} ) {
	    		my ($indel_type_file,$indel_bases,$qual) = split('\^\^\^',$indel_entry);
	    		next unless $indel_bases eq $newindel_varbases; #only look at relevant match (i.e. don't match +A and +AA insertions at same base)
		    	 
	    		#If there are multiple entries		
				my @matches = ();
				my %matches = ();
				my %rs = ();
				my $rs;
	
		    	for my $match_entry (sort @{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}{$indel_entry}}) {
		    		#Need to confirm that indel matches wrt all indel type and bases
		    		my ($match_type,$local_rs,$match_bases) = split('\^\^\^',$match_entry);
		    		next unless $match_type eq $indel_type_file; #Don't 'match' insertions and deletions
	    			#For dels we check the bases match
	    			my ($just_match_bases) = $match_bases =~ /^([ATCGN]+)/;
		    		
	    			$just_match_bases =~ s/N//g; #few entries have N's present in string
	    			
	    			if ($match_type eq 'DEL') {    				
	    				if ($just_match_bases ne $indel_bases) {
	    					modules::Exception->throw("ERROR: deletion $indel_entry and match $match_entry have different bases");
	    				}
	    			} else {
	    				#Don't require exact bases for insertion match; 
	    				#next if ($just_match_bases ne $indel_bases);
	    			}
	    			$indel_filter{'filtermatch'} = 1;
					$indel_filter{'filterpass'} = $filterpass == 1?1:0;
					#Set the rs as matched even if there are no allele frequencies reported
					$rs = $local_rs;
					if (!defined $just_match_bases) {
						modules::Exception->throw("ERROR: Can't get bases from match entry $match_entry");
					}
					#If there's allele frequency info
			    	if ($match_bases =~ /\(([0-9\.]+):([0-9\.]+)\)/) {
			    		if ($match_type eq 'INS' && $just_match_bases eq $indel_bases) {
			    			#make sure the inserted bases match
							$matches{$1} = $2;
							$rs{$1} = $local_rs;
			    		} elsif ($match_type eq 'DEL') {
			    			#deletion we already know bases match
			    			$matches{$1} = $2;
			    			$rs{$1} = $local_rs;
			    		} 
			    	} 
			    	
			    	$matched = 1; #At least a coordinate match
					push @matches, $match_entry;
		    	}
		    	
			    if (keys %matches) {
			    	#Get the highest ref_allele value for multiple matches; we use allele freq in determining pass
			    	my ($match_ref_allele) = reverse sort {$a<=>$b} keys %matches; 
			    	$indel_filter{'attribute'} = 'rs='.$rs{$match_ref_allele} . ';ref_allele_freq=' . $match_ref_allele . ';var_allele_freq=' . $matches{$match_ref_allele};
			    } elsif ($matched) {
			    	$indel_filter{'attribute'} = 'rs='.$rs;
			    } 
		    	
		    	if (@matches) {
			    	my $match_entries = join(",",@matches);
					print MATCH join("\t",
										$newindel_chr,
										$newindel_start,
										$newindel_end,
										$indel_entry,
										$match_entries
										) . "\n";
		    	} else {
					print NOMATCH join("\t",
										$newindel_chr,
										$newindel_start,
										$newindel_end,
										$indel_entry,
										'NO_MATCH'
										) . "\n";
		    	}
		    	
		    }
			
			
			#Never found an allele and indel type match so write to nomatch file
			if (!$matched) {
				my ($indel_entry) = keys %{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}};
			} 
		} elsif ($filter_common) {
			#Here any kind of match will count as a pass (ins or deletion)
			for my $indel_entry ( keys %{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}} ) {
				my $match_entry = join(",",@{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}{$indel_entry}});
				$indel_filter{'attribute'} = $match_entry;
				$indel_filter{'filtermatch'} = 1;
				$indel_filter{'filterpass'} = $filterpass == 1?1:0;
				print MATCH join("\t",
										$newindel_chr,
										$newindel_start,
										$newindel_end,
										$indel_entry,
										$match_entry
										) . "\n";
			}
		} elsif ($filter_exac) {
			#Require type to match for filter_exac
			for my $indel_entry ( keys %{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}} ) {
				for my $match_entry (sort @{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}{$indel_entry}}) {
					my ($indel_type_file,$indel_bases,$qual) = split('\^\^\^',$indel_entry);
					$match_entry =~ s/^\d+://; #Get rid of leading coord used to make entry unique
					my $match_type = $match_entry =~ /^\+/?'INS':'DEL';
		    		next unless $match_type eq $indel_type_file; #Don't 'match' insertions and deletions
		    		$indel_filter{'attribute'} = 'exac='.$match_entry;
		    		$indel_filter{'filtermatch'} = 1;
					$indel_filter{'filterpass'} = $filterpass == 1?1:0;
					print MATCH join("\t",
										$newindel_chr,
										$newindel_start,
										$newindel_end,
										$indel_entry,
										$match_entry
										) . "\n";
				}
			}
		} elsif ($filter_gnomad) {
			#Require type to match for filter_gnomad
			for my $indel_entry ( keys %{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}} ) {
				for my $match_entry (sort @{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}{$indel_entry}}) {
					my ($indel_type_file,$indel_bases,$qual) = split('\^\^\^',$indel_entry);
		    		$match_entry =~ s/^\d+://; #Get rid of leading coord used to make entry unique
					my $match_type = $match_entry =~ /^\+/?'INS':'DEL';
		    		next unless $match_type eq $indel_type_file; #Don't 'match' insertions and deletions
		    		$indel_filter{'attribute'} = 'gnomad='.$match_entry;
		    		$indel_filter{'filtermatch'} = 1;
					$indel_filter{'filterpass'} = $filterpass == 1?1:0;
					print MATCH join("\t",
										$newindel_chr,
										$newindel_start,
										$newindel_end,
										$indel_entry,
										$match_entry
										) . "\n";
				}
			}
		} elsif ($filter_clinvar) {
                        for my $indel_entry ( keys %{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}} ) {
                                for my $match_entry (sort @{$overlap_match->{PASS}{$newindel_chr}{$newindel_start}{$newindel_end}{$indel_entry}}) {
                                        my ($indel_type_file,$indel_bases,$qual) = split('\^\^\^',$indel_entry);
                                        my $clinvar_match = $match_entry;
                                        $indel_filter{'attribute'} = 'clinvar='.$clinvar_match;
                                        $indel_filter{'filtermatch'} = 1;
                                        $indel_filter{'filterpass'} = $filterpass == 1?1:0;
                                        print MATCH join("\t",
                                                        $newindel_chr,
                                                        $newindel_start,
                                                        $newindel_end,
                                                        $indel_entry,
                                                        $match_entry
                                                        ) . "\n";
                                }
                        }
		}
	} elsif (exists $overlap_nomatch->{FAIL}{$newindel_chr}{$newindel_start}{$newindel_end}) {
		my ($indel_entry) = keys %{$overlap_nomatch->{FAIL}{$newindel_chr}{$newindel_start}{$newindel_end}};
		print NOMATCH join("\t",
								$newindel_chr,
								$newindel_start,
								$newindel_end,
								$indel_entry,
								'NO_MATCH'
								) ."\n";
		next if $filter_common; #Don't enter non matching entries into db
		next if $filter_exac;
		next if $filter_gnomad;
		next if $filter_clinvar;
		#If the indel didn't match
		$indel_filter{'filtermatch'} = 0;
		$indel_filter{'filterpass'} = $filterpass == 0?1:0;
	} else {
		my $snp_str = $newindel_chr." ".$newindel_start;
		modules::Exception->throw("ERROR: Can't find snp $snp_str in the fail or pass list\n");
	}

	if (!exists $indel_filter{attribute}) {
		$indel_filter{attribute} = 'N/A';
	}
		
	my $mutant_indel_key =  $newindel_chr.":".$newindel_start;
		
	my $end = $newindel_end;
	my $var_base = $newindel->{var_bases};
		
	if ($dbsnp_indel_flag) {
		if ($indel_filter{attribute} eq 'N/A') {
			$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_filter_dbsnp_indel_string} = 'novel';		
		} else {
			$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_filter_dbsnp_indel_string} = $indel_filter{'attribute'};
			$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{pass} = 0;
		}
	} elsif ($filter_common && $indel_filter{'filterpass'} == 1) {
		$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_filter_common_string} = $indel_filter{'attribute'};
	} elsif ($filter_exac && $indel_filter{'filterpass'} == 1) {
		$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_filter_exac_indel_string} = $indel_filter{'attribute'};
	} elsif ($filter_gnomad && $indel_filter{'filterpass'} == 1) {
		$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_filter_gnomad_indel_string} = $indel_filter{'attribute'};
	} elsif ($filter_clinvar && $indel_filter{'filterpass'} == 1) {
    	$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{indel_filter_clinvar_indel_string} = $indel_filter{'attribute'};
  	}

		
	#Only load passed variants except with dbsnp
	if ($indel_filter{'filterpass'} == 1) {
		$mutant_indel_data{$mutant_indel_key}{$end}{$var_base}{pass} = 1;
	} 
		
	
	
	
}


close MATCH;
close NOMATCH;

	
my $file_name = $sample_name . '_' . $runid . '.'.$coord_filter.'.indel.xml';
if ($OPT{chr}) {
	$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_indel_data, -chr=>$OPT{chr});
	if ($last_chr) { #Chr is Y
		my $full_xml_file = $run_dir . '/conf/'.$file_name;
		$var_xml->split_xml(-file_name=>$full_xml_file);
	}
	
}  else {
	$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%mutant_indel_data, -chr=>'all');
	my $full_xml_file = $run_dir . '/conf/'.$file_name;
	$var_xml->split_xml(-file_name=>$full_xml_file);
}

print STDERR "Total match count: $match_count\n";
print STDERR "Total nomatch count: $nomatch_count\n";

