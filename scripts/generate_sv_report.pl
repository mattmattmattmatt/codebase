#! /usr/bin/perl -w

use strict;
use modules::Report;
use modules::Adaptors::Filter;
use modules::Adaptors::BulkInsert;
use modules::Annotation;
use modules::VariantXML;
use Getopt::Long;
use modules::Exception;
use modules::Utils;
use modules::Pipeline;
use modules::SystemCall;
use File::Basename;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

GetOptions(\%OPT, 
		    "help|h",
		    "man|m",
		 	"runid=i",
		    "chr=s",
		 	"start=i",
		 	"end=i",
		 	"ref=s",
		 	"annotation_file=s",
		 	"exon_coord_file=s",
		    "tsv_file=s",
		    "pass_summary_file=s",
		    "pass_file=s",
		    "genefilter_files=s",
			"report_xml=s",
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{tsv_file} || !$OPT{runid} || !$OPT{annotation_file} || !$OPT{exon_coord_file} || !$OPT{pass_summary_file} || !$OPT{pass_file});

	   
=pod

=head1 SYNOPSIS

generate_sv_report.pl -runid runid -ref ref_name(default=hs37d5) -annotation_file annotation_file -exon_coord_file exon_coord_file -debug print_stderr -filter_summary_file filter_summary_output_file -tsv_file output_tsv_file -chr specific_chr(default=all) -start chr_start_coord(default=all) -end chr_end_coord(default=all) -pass_summary_file pass_rows_summary [options]

Required flags: -runid -tsv_file -annotation_file -exon_coord_file -pass_summary_file -pass_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

generate_sv_report.pl -> Script to generate an sv report

=head1 DESCRIPTION

Feb 18, 2011

a script that ...

=head1 AUTHOR

Dan Andrews

=head1 EXAMPLE

generate_sv_report.pl

=cut

my $pipe_conf = modules::Pipeline::get_pipe_conf();

my $sys_call = modules::SystemCall->new();

#Lots of variables to load...
my $report_file = $OPT{report_file};
my $output_file = $OPT{output_file};
my $summary_file = $OPT{pass_summary_file};
my $pass_file = $OPT{pass_file};
my $pipe_config = modules::Pipeline::get_pipe_conf();

my @filter_files = ();
if (defined $OPT{genefilter_files}) {
	@filter_files = split(",",$OPT{genefilter_files});
}

my $runid = $OPT{runid};
my $tsv_file = $OPT{tsv_file};
my $chr = defined $OPT{chr}?$OPT{chr}:'';
my $start = defined $OPT{start}?$OPT{start}:'';
my $end = defined $OPT{end}?$OPT{end}:'';
my $annotation_file = $OPT{annotation_file};
my $ref = defined $OPT{ref}?$OPT{ref}:'hs37d5';
my $source_type = modules::Pipeline::get_source_type(-run_id=>$runid);
my $sample_type = modules::Pipeline::get_sample_type(-run_id=>$runid);
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);

my $filter = my $allele_filter;

my $splice_size = $pipe_conf->read($source_type,'cutoffs','splice_length_cutoff');
my $exon_coord_file = $OPT{exon_coord_file};

if ( !-e $annotation_file ) {
	modules::Exception->throw("File $annotation_file doesn't exist");	
}


if ( !-e $exon_coord_file ) {
	modules::Exception->throw("File $exon_coord_file doesn't exist");	
}

my ($run_obj) = modules::Adaptors::Run->search('id' => $runid);
unless (defined $run_obj) {
    modules::Exception->throw("Did not get all necessary things from database [run_id:$runid]");
}

#Get the gene column name
my $config = defined($OPT{report_xml}) ? modules::Pipeline->get_report_conf(-report_xml=>$OPT{report_xml}) : modules::Pipeline->get_report_conf();
if (!$config->exists('sv','sample_types',$sample_type)) {
   	modules::Exception->throw("ERROR: Cannot get annotation columns for sv $sample_type");
}

my $organism = $pipe_conf->read($source_type,'organism');
my ($gene_col_name) = split(",",$config->read('snv','common','human_file_annotations'));

#generate the gene mapper 
my $annotation = modules::Annotation->new(-annotation_file=>$annotation_file,-exon_coord_file=>$exon_coord_file,-splice_size=>$splice_size,-sample_type=>$sample_type,-mutant_type=>'sv',-organism=>$organism,-config => $config);

my @sv_callers = split(",",$pipe_conf->read('common','sv_callers'));

my @sv_types = split(",",$pipe_conf->read('common','sv_types'));

#Keep track of summary reports for generating final combined report
my @final_reports = ();

for my $sv_caller (@sv_callers) {
	(my $tsv_file_final = $tsv_file) =~ s/SVCALLER/$sv_caller/;
	(my $pass_file_final = $pass_file) =~ s/SVCALLER/$sv_caller/;
	(my $summary_file_final = $summary_file) =~ s/SVCALLER/$sv_caller/;
	
	push @final_reports, $tsv_file_final;
	#Create the snp report params
	my %sv_report_param = ( -run => $run_obj,
		                     -gene_mapper => $annotation,
		                     -sample_type => $sample_type,
		                     -sv=>1,
		                     -sv_caller=>$sv_caller,
		                     -source_type => $source_type,
		                     -gene_col_name => $gene_col_name
		                     );
	
	
	if (@filter_files) {
		$sv_report_param{-filter_files} = \@filter_files
	}
	my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
	
	$sv_report_param{-confdir} = $run_dir . '/conf';
	
	my $sv_report = modules::Report->new(%sv_report_param);
	
	#Load the filters, this is specific to each sample group
	if($OPT{report_xml}) {
		$sv_report->load(-report_xml=>$OPT{report_xml});
	} else{
		$sv_report->load();
	}
		
	if ($chr && $start) {
		$sv_report->generate_pass_fail(-chr=>$chr, -start=>$start, -end=>$end);
	} elsif ($chr) {
		$sv_report->generate_pass_fail(-chr=>$chr);
	} else {
		$sv_report->generate_pass_fail();
	}
	
	
	$sv_report->print_to_files(-tsv_file=>$tsv_file_final,-pass_file=>$pass_file_final);
	$sv_report->summarize(-summary_file=>$summary_file_final);
	
	#First the passed rows
	my $pass_row_data = $sv_report->get_pass_rows();
	my $rare_allele_rows = $sv_report->get_allele_rows();
	my $no_freq_rows = $sv_report->get_no_freq_rows();
	my $fail_row_data = $sv_report->get_fail_rows();
	
	
	my %rows = (
				'pass'=>$pass_row_data,
				'allele'=>$rare_allele_rows,
				'no_freq'=>$no_freq_rows,
				'fail'=>$fail_row_data
				);
	
	my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');
	
	
	for my $row_type (keys %rows) {
		my $row_count = 0;
		my $file_name = $sample_name . '_' . $runid . '.'.$sv_caller .'.sv.xml'; 
		$var_xml->load_xml_keys(-file_name=>$file_name,-sv=>1);
		for my $chr1 (sort keys %{$rows{$row_type}})	{
			for my $chr2 (sort keys %{$rows{$row_type}->{$chr1}})	{
			    for my $start_coord (sort {$a<=>$b} keys %{$rows{$row_type}->{$chr1}{$chr2}}) {
			    	for my $end_coord (sort {$a<=>$b} keys %{$rows{$row_type}->{$chr1}{$chr2}{$start_coord}}) {
			    		for my $sv_caller_tmp (keys %{$rows{$row_type}->{$chr1}{$chr2}{$start_coord}{$end_coord}}) {
							$row_count++;
			    		}
			    	}
				    
				}
			
			
			}
		}	
		
		print STDERR "Reporting $row_count $row_type rows\n";
		
	}
	
	my $final_xml = $sample_name . '_' . $runid .'.'.$sv_caller. '.report_sv.sv.xml';
	$sv_report->final_xml(-xml_file=>$final_xml);
}

#Create tmp dir for overlaps
my $basedir = dirname($final_reports[0]);
my $tmpdir = $basedir .'/tmp';

if (-d $tmpdir) {
	$sys_call->run("rm -Rf $tmpdir"); 
}
$sys_call->run("mkdir $tmpdir");

my $svndir = $ENV{'SVNDIR'};
if (!-d "$svndir/utils") {
	modules::Exception->throw("ERROR: Can't open dir $svndir\n");
}

my $overlap_bin = "$svndir/utils/overlap_files.pl";
my $max_distance = $pipe_config->read($source_type,'cutoffs','sv_max_distance');

my %common_sv_ids = ();

my $combined_event_count = 0;

if (@final_reports == 2) {
	for my $sv_type (@sv_types) {
		#fussy grepping for tab
		my $regex = "'".'\\t'.${sv_type}.'\\t'."'";
		if ($sv_type eq 'tra') {
			#here we need the second chr
			`grep -P $regex $final_reports[0] | cut -d'\t' -f 1,2,3,4,7 | sort | uniq > $tmpdir/file1.$sv_type`;
			`grep -P $regex $final_reports[1] | cut -d'\t' -f 1,2,3,4,7 | sort | uniq > $tmpdir/file2.$sv_type`;
		} else {
			`grep -P $regex $final_reports[0] | cut -d'\t' -f 1,2,4,7 | sort | uniq > $tmpdir/file1.$sv_type`;
			`grep -P $regex $final_reports[1] | cut -d'\t' -f 1,2,4,7 | sort | uniq > $tmpdir/file2.$sv_type`;
		}
		
		if ($sv_type eq 'tra') {
			#Use custom code to overlap both coordinates
			open(FILE1,"$tmpdir/file1.tra") || modules::Exception->throw("Can't open file file1.tra\n");
			open(FILE2,"$tmpdir/file2.tra") || modules::Exception->throw("Can't open file file2.tra\n");
			
			my %trans = ();
			
			while (<FILE1>) {
				chomp;
				my ($chr1,$coord1,$chr2,$coord2,$id) = split("\t");
				my $chr_pair;
				my $key;
				if ($chr1 lt $chr2) {
					$chr_pair = $chr1.':'.$chr2;
					$key = "$coord1:$coord2:$id";
				} else {
					$chr_pair = $chr2.':'.$chr1;
					$key = "$coord2:$coord1:$id";
				}
				$trans{$chr_pair}{first}{$key}++;
			}
			
			while (<FILE2>) {
				chomp;
				my ($chr1,$coord1,$chr2,$coord2,$id) = split("\t");
				my $chr_pair;
				my $key;
				if ($chr1 lt $chr2) {
					$chr_pair = $chr1.':'.$chr2;
					$key = "$coord1:$coord2:$id";
				} else {
					$chr_pair = $chr2.':'.$chr1;
					$key = "$coord2:$coord1:$id";
				}
				$trans{$chr_pair}{second}{$key}++;
			}
			
			
			open(TRA,">$tmpdir/combined.tra") || modules::Exception->throw("Can't open file to write $tmpdir/combined.tra\n");
			
			
			for my $combined_coord (keys %trans) {
				if (keys %{$trans{$combined_coord}} == 2) {
					
					for my $first_str (keys %{$trans{$combined_coord}{first}}) {
						my ($first_coord1,$first_coord2,$first_id) = split(':',$first_str);
					
						for my $second_str (keys %{$trans{$combined_coord}{second}}) {
							my ($second_coord1,$second_coord2,$second_id) = split(':',$second_str);
					
							if (abs($first_coord1-$second_coord1) < $max_distance && abs($first_coord2-$second_coord2) < $max_distance) {
								my ($chr1,$chr2) = split(':',$combined_coord);
								
								$common_sv_ids{tra}{$combined_event_count}{$first_id}++;
								$common_sv_ids{tra}{$combined_event_count}{$second_id}++;
								
								print TRA join("\t",
												$chr1,
												$first_coord1,
												$chr2,
												$first_coord2,
												$first_id,
												$second_id
												) ."\n";
							}
						}
					}
					$combined_event_count++;
				}
			}
			
		} elsif ($sv_type eq 'ins') {
			#Lumpy doesn't call ins
			$sys_call->run("cp $tmpdir/file1.ins $tmpdir/combined.ins");
		} else {
			#Dels/dups/invs get normal handling
			$sys_call->run("$overlap_bin -just_overlap -max_distance $max_distance -ref $tmpdir/file1.$sv_type -coord $tmpdir/file2.$sv_type > $tmpdir/combined.$sv_type"); 
			open(OVERLAP,"$tmpdir/combined.$sv_type") || modules::Exception->throw("Can't open file $tmpdir/combined.$sv_type\n");
			my %repeats = ();
			
			while (<OVERLAP>) {
				chomp;
				my @fields = split("\t");
				my @file2_ids = split('\^\^\^',$fields[4]);
				
				#Group pairs of event which all overlap (i.e. different sides of breakpoints)
				if (!exists $repeats{$fields[4]}) {
					$combined_event_count++;
				} 
				
				$repeats{$fields[4]}++;
				
				
				$common_sv_ids{$sv_type}{$combined_event_count}{$fields[3]}++;
				
				for my $file2_id (@file2_ids) {
					$common_sv_ids{$sv_type}{$combined_event_count}{$file2_id}++;
				}
				
			}	
		}
	}
} else {
	modules::Exception->throw("ERROR: No handling for three SV callers\n");
}

my $max_event_count = $combined_event_count + 1;


my $uniq_id = 'sv_id';
my $status = 'final_status';
my $exon_col_name = 'breakpoint_exon';

my %pass_lines = ();
my %fail_lines = ();
my @headers = ();
my $sv_length_cutoff = 100000;


#Now we combine the reports to generate a final unified prioritised report
for my $file (@final_reports) {
	#Want to order rows using the following algorithm
	#1) SVs found in both and < 100kb
		#For each sv type
		#i) Two callers Novel/rare SVs where breakpoint overlaps an exon (for all sv types) (high)
		#ii) Two callers Novel/rare SVs which overlaps gene but not exon bp (for del/dup) (medium)
		#iii) One caller Novel/rare SVs where breakpoint overlaps an exon (for all sv types) (medium)
		#iv) Novel/rare with no exon but gene (intron) overlap (for inv, tra) (low)
		#v) Novel/rare with no exon or gene overlap (for inv, tra) (low)
		#vi) Lowest priority -> everything left (fails,etc) 
		
	#2) SVs unique to one caller and SVs above > 100kb and delly insertions
		#Same as above
	
	open(FILE,"$file") || modules::Exception->throw("Can't open file $file\n");
	
	my $sv_id_index;
	my $gene_index;
	my $exon_index;
	my $length_index;
	my $status_index;
	
	while (<FILE>) {
		chomp;
		
		next if $_ =~ /^$/;
		my @fields = split("\t");
		next unless $fields[0] =~ /^[0-9XYc]/;
	
		
		
		if ($fields[0] =~ /^chr_1/) {
			@headers = @fields;
			my $header_count = 0;
			for my $header (@headers) {
				if ($header eq $uniq_id) {
					$sv_id_index = $header_count;
				} elsif ($header eq $exon_col_name) {
					$exon_index = $header_count;
				} elsif ($header eq $gene_col_name) {
					$gene_index = $header_count;
				} elsif ($header eq 'length') {
					$length_index = $header_count;
				} elsif ($header eq $status) {
					$status_index = $header_count;
				} 						
				$header_count++;
			}
			next;
		}
		
				
		my ($chr1,$coord1,$chr2,$coord2,$sv_caller,$sv_type,$id) = @fields;
		
		
		
		my $fail;
		if ($fields[$status_index] =~ /PASS/) {
			$fail = 0;
		} else {
			$fail = 1;
		}
		
		if ($fields[$length_index] > $sv_length_cutoff) {
			#only dels and dups fail for length
			if ($sv_type ne 'inv') {
				if ($fail) {
					$fields[$status_index] .= ",length_fail";
				} else {
					$fields[$status_index] = "FAIL_REASONS: length_fail";
				}
				
				$fail = 1;
				
				
			}
		} 
		
		my $sv_id = $fields[$sv_id_index];
		
		
		#Check if event is combined
		my $combined_event = 0;
		
		for my $event_count ( keys %{$common_sv_ids{$sv_type}} ) {
			for my $sv_id_local (keys %{$common_sv_ids{$sv_type}{$event_count}}) {
				if ($sv_id_local eq $sv_id) {
			    	$combined_event = $event_count;
				}
			}
		}
		
		my $exon_event = 1;
		
		if ($fields[$exon_index] eq 'NO_EXON_BP') {
			$exon_event = 0;
		} 
			
		my $gene_event = 1;
		
		if ($fields[$gene_index] eq 'No gene overlap') {
			$gene_event = 0;
		} 
		
		if ($fail || $combined_event == 0) {
			push @{$fail_lines{$sv_type}{$combined_event}{$exon_event}{$gene_event}{lines}}, join("\t",@fields);
		} else {
			push @{$pass_lines{$sv_type}{$combined_event}{$exon_event}{$gene_event}{lines}}, join("\t",@fields);
		}
		
		
	}
	close (FILE);
	
	#$xml_data{$xml_key}{$end_coord}{$sv_caller}{$sv_type}{$column_name} = $column_value;
}

#Clean up cases where one is passed and one is failed -> happens when one bp overlaps dgv and other doesn't
for my $sv_type ( sort keys %pass_lines ) {
	for my $event_id ( sort {$a<=>$b} keys %{$pass_lines{$sv_type}} ) {
		if (exists $fail_lines{$sv_type}{$event_id}) {
			for my $exon_event (keys %{$pass_lines{$sv_type}{$event_id}}) {
				for my $gene_event (keys %{$pass_lines{$sv_type}{$event_id}{$exon_event}}) {
					for my $line (@{$pass_lines{$sv_type}{$event_id}{$exon_event}{$gene_event}{lines}}) {
						push @{$fail_lines{$sv_type}{$event_id}{$exon_event}{$gene_event}{lines}}, $line;
					}
				}
			}
			delete $pass_lines{$sv_type}{$event_id};
		}
	}
}


(my $tsv_combined_final = $tsv_file) =~ s/SVCALLER/combined/;
open(COMBINED,">$tsv_combined_final") || modules::Exception->throw("Can't open file to write $tsv_combined_final/\n");

print COMBINED join("\t",
					"SV Caller Count",
					"Total SV Calls",
					"Exon Overlap",
					"Gene Overlap",
					"Event ID",
					@headers
					) ."\n\n";


#Add 5 fields to lines (one left for text; rest for ---)
my $divider_count = @headers + 4;

#Keep track of lines reported to avoid duplicates
my %reported_svs = ();

#First passed exon SVs
for my $sv_type ( sort keys %pass_lines ) {
	print COMBINED "\nHIGH_PRIORITY $sv_type EXON BP\t" . "---\t" x $divider_count ."\n\n";
    for my $event_id ( sort {$a<=>$b} keys %{$pass_lines{$sv_type}} ) {
		
    	if (exists $pass_lines{$sv_type}{$event_id}{1} && exists $pass_lines{$sv_type}{$event_id}{1}{1}) {
			my @lines = ();
			my %uniq_count = ();	
    		my %sv_callers = ();
    		for my $line (@{$pass_lines{$sv_type}{$event_id}{1}{1}{lines}}) {
    			push @lines, $line;
    			$reported_svs{$line}++;
    			my @tmp = split("\t",$line);
    			$uniq_count{$tmp[6]}++;
    			$sv_callers{$tmp[4]}++;
    		}
    		#Sometimes only one end may overlap exon for trans so need to include both
    		if (exists $pass_lines{$sv_type}{$event_id}{0} && exists $pass_lines{$sv_type}{$event_id}{0}{1}) {
    			for my $line (@{$pass_lines{$sv_type}{$event_id}{0}{1}{lines}}) {
    				push @lines, $line;
    				$reported_svs{$line}++;
    				my @tmp = split("\t",$line);
    				$uniq_count{$tmp[6]}++;
    				$sv_callers{$tmp[4]}++;
    			}
    		}
    		
    		my $uniq_count = keys %uniq_count;
    		my $sv_callers = keys %sv_callers;
    		print COMBINED join("\t",
    						$sv_callers,
    						$uniq_count,
    						'EXON_BP',
    						'GENE_OVERLAP',
    						$event_id,
    						$lines[0]
    						) . "\n";
    		
    		shift @lines;
    		
    		for my $line (@lines) {
	    		print COMBINED join("\t",
	    					'',
	    					'',
	    					'',
	    					'',
	    					$event_id,
	    					$line
	    					) ."\n";
    			
    		} 
    	} 
	}
}


#Now passed gene SVs for dels and dups
for my $sv_type ( qw(del dup) ) {
	print COMBINED "\nMEDIUM_PRIORITY $sv_type GENE_OVERLAP\t" . "---\t" x $divider_count ."\n\n";
    for my $event_id ( sort {$a<=>$b} keys %{$pass_lines{$sv_type}} ) {
		
    	if (exists $pass_lines{$sv_type}{$event_id}{0} && exists $pass_lines{$sv_type}{$event_id}{0}{1}) {
			my @lines = ();
			my %uniq_count = ();	
    		my %sv_callers = ();
    		my $reported = 0;
    		for my $line (@{$pass_lines{$sv_type}{$event_id}{0}{1}{lines}}) {
    			if (exists $reported_svs{$line}) {
    				$reported = 1;
    			}
    			$reported_svs{$line}++;
    			push @lines, $line;
    			my @tmp = split("\t",$line);
    			$uniq_count{$tmp[6]}++;
    			$sv_callers{$tmp[4]}++;
    		}
    		#Skip mixed cases we've already reported
    		next if $reported;
    		my $uniq_count = keys %uniq_count;
    		my $sv_callers = keys %sv_callers;
    		print COMBINED join("\t",
    						$sv_callers,
    						$uniq_count,
    						'NO_EXON_BP',
    						'GENE_OVERLAP',
    						$event_id,
    						$lines[0]
    						) . "\n";
    		
    		shift @lines;
    		
    		for my $line (@lines) {
	    		print COMBINED join("\t",
	    					'',
	    					'',
	    					'',
	    					'',
	    					$event_id,
	    					$line
	    					) ."\n";
    			
    		} 
    	} 
	}
}

#Now exon single sv caller SVs
for my $sv_type ( sort keys %fail_lines ) {
    print COMBINED "\nMEDIUM_PRIORITY $sv_type SINGLE_TOOL EXON_BP_OR_DGV/LENGTH_FAIL\t" . "---\t" x $divider_count ."\n\n";
    for my $event_id ( sort {$a<=>$b} keys %{$fail_lines{$sv_type}} ) {
		
    	if (exists $fail_lines{$sv_type}{$event_id}{1} && exists $fail_lines{$sv_type}{$event_id}{1}{1}) {
			my @lines = ();
			my %uniq_count = ();	
    		my %sv_callers = ();
    		my $reported = 0;
    		
    		for my $line (@{$fail_lines{$sv_type}{$event_id}{1}{1}{lines}}) {
    			if (exists $reported_svs{$line}) {
    				$reported = 1;
    			}
    			push @lines, $line;
    			$reported_svs{$line}++;
    			my @tmp = split("\t",$line);
    			$uniq_count{$tmp[6]}++;
    			$sv_callers{$tmp[4]}++;
    		}
    		next if $reported;
    		
    		my $uniq_count = keys %uniq_count;
    		my $sv_callers = keys %sv_callers;
    		
    		if ($event_id == 0) {
    			for my $line (@lines) {
    				print COMBINED join("\t",
	    						1,
	    						1,
	    						'EXON_BP',
	    						'GENE_OVERLAP',
	    						$max_event_count,
	    						$line
	    						) . "\n";
	    			$max_event_count++;
    			}
    		} else {
	    		print COMBINED join("\t",
	    						$sv_callers,
	    						$uniq_count,
	    						'EXON_BP',
	    						'GENE_OVERLAP',
	    						$event_id,
	    						$lines[0]
	    						) . "\n";
	    		
	    		shift @lines;
	    		
	    		for my $line (@lines) {
		    		print COMBINED join("\t",
		    					'',
		    					'',
		    					'',
		    					'',
		    					$event_id,
		    					$line
		    					) ."\n";
	    			
	    		}
	    			
	    	} 
    		
    	} 
	}
}

#Now non-exon non-gene passed SV (happens for bp SVs like trans and inversion)
for my $sv_type ( qw(tra inv) ) {
	print COMBINED "\nLOW_PRIORITY $sv_type NO_EXON/GENE_OVERLAP\t" . "---\t" x $divider_count ."\n\n";
	for my $event_id ( sort {$a<=>$b} keys %{$pass_lines{$sv_type}} ) {
			
	    if (exists $pass_lines{$sv_type}{$event_id}{0} && exists $pass_lines{$sv_type}{$event_id}{0}{1}) {
			my @lines = ();
			my %uniq_count = ();	
    		my %sv_callers = ();
    		my $reported = 0;
    		for my $line (@{$pass_lines{$sv_type}{$event_id}{0}{1}{lines}}) {
    			if (exists $reported_svs{$line}) {
    				$reported = 1;
    			}
    			$reported_svs{$line}++;
    			push @lines, $line;
    			my @tmp = split("\t",$line);
    			$uniq_count{$tmp[6]}++;
    			$sv_callers{$tmp[4]}++;
    		}
    		#Skip mixed cases we've already reported
    		next if $reported;
    		my $uniq_count = keys %uniq_count;
    		my $sv_callers = keys %sv_callers;
    		print COMBINED join("\t",
    						$sv_callers,
    						$uniq_count,
    						'NO_EXON_BP',
    						'GENE_OVERLAP',
    						$event_id,
    						$lines[0]
    						) . "\n";
    		
    		shift @lines;
    		
    		for my $line (@lines) {
	    		print COMBINED join("\t",
	    					'',
	    					'',
	    					'',
	    					'',
	    					$event_id,
	    					$line
	    					) ."\n";
	    			
	    	} 
		} 
	}
}

#Now non-exon non-gene passed SV (happens for bp SVs like trans and inversion)
for my $sv_type ( qw(tra inv) ) {
	print COMBINED "\nLOW_PRIORITY $sv_type NO_EXON/NO_GENE\t" . "---\t" x $divider_count ."\n\n";
    for my $event_id ( sort {$a<=>$b} keys %{$pass_lines{$sv_type}} ) {
		
    	if (exists $pass_lines{$sv_type}{$event_id}{0} && exists $pass_lines{$sv_type}{$event_id}{0}{0}) {
			my @lines = ();
			my %uniq_count = ();	
    		my %sv_callers = ();
    		my $reported = 0;
    		for my $line (@{$pass_lines{$sv_type}{$event_id}{0}{0}{lines}}) {
    			if (exists $reported_svs{$line}) {
    				$reported = 1;
    			}
    			$reported_svs{$line}++;
    			push @lines, $line;
    			my @tmp = split("\t",$line);
    			$uniq_count{$tmp[6]}++;
    			$sv_callers{$tmp[4]}++;
    		}
    		#Skip mixed cases we've already reported
    		next if $reported;
    		my $uniq_count = keys %uniq_count;
    		my $sv_callers = keys %sv_callers;
    		print COMBINED join("\t",
    						$sv_callers,
    						$uniq_count,
    						'NO_EXON_BP',
    						'NO_GENE',
    						$event_id,
    						$lines[0]
    						) . "\n";
    		
    		shift @lines;
    		
    		for my $line (@lines) {
	    		print COMBINED join("\t",
	    					'',
	    					'',
	    					'',
	    					'',
	    					$event_id,
	    					$line
	    					) ."\n";
    			
    		} 
    	}
    }
}

#Now fail or low priority SVs other than exon cases -> don't bother with dividing lines except for exon overlaps
print COMBINED "\nLOWEST_PRIORITY\t" . "---\t" x $divider_count ."\n\n";
for my $sv_type ( sort keys %fail_lines ) {
    for my $event_id ( sort {$a<=>$b} keys %{$fail_lines{$sv_type}} ) {
		
		for my $gene_overlap (sort {$b<=>$a} keys %{$fail_lines{$sv_type}{$event_id}{0}}) {
			my @lines = ();
			my %uniq_count = ();	
    		my %sv_callers = ();
    		my $reported = 0;
    		my $gene_str = $gene_overlap == 1?'GENE_OVERLAP':'NO_GENE';
    		for my $line (@{$fail_lines{$sv_type}{$event_id}{0}{$gene_overlap}{lines}}) {
    			if (exists $reported_svs{$line}) {
    				$reported = 1;
    			}
    			push @lines, $line;
    			$reported_svs{$line}++;
    			my @tmp = split("\t",$line);
    			$uniq_count{$tmp[6]}++;
    			$sv_callers{$tmp[4]}++;
    		}
    		next if $reported;
    		
    		my $uniq_count = keys %uniq_count;
    		my $sv_callers = keys %sv_callers;
    		
    		if ($event_id == 0) {
    			for my $line (@lines) {
    				print COMBINED join("\t",
	    						1,
	    						1,
	    						'EXON_BP',
								$gene_str,
								$max_event_count,
	    						$line
	    						) . "\n";
	    			$max_event_count++;
    			}
    		} else {
	    		print COMBINED join("\t",
	    						$sv_callers,
	    						$uniq_count,
	    						'EXON_BP',
	    						$gene_str,
	    						$event_id,
	    						$lines[0]
	    						) . "\n";
	    		
	    		shift @lines;
	    		
	    		for my $line (@lines) {
		    		print COMBINED join("\t",
		    					'',
		    					'',
		    					'',
		    					'',
		    					$event_id,
		    					$line
		    					) ."\n";
	    			
	    		}
	    			
	    	} 
    	} 
	}
}


