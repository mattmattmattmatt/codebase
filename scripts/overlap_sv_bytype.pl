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
		   "coord_file=s",
		   "coord_filtername=s",
		   "overlap_outfile=s",
		   "ref_file=s",
		   "sv_type=s",
		   "runid=i"
	   		);
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid} || !$OPT{coord_file} || !$OPT{overlap_outfile} || !$OPT{ref_file});

	   
=pod

=head1 SYNOPSIS

overlap_svs_bytype.pl -sv_type comma_delim_list_of_Var_types(default=del,dup,inv)-runid runid -overlap_outfile output_overlap_file -ref_file input_overlap_file -coord_file formatted_overlap_file [options]

Required flags: -coord_filtername -coord_file -runid -overlap_outfile -ref_file

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

overlap_sv_bytype.pl -> Script to run a overlaps for svs versus dgv

=head1 DESCRIPTION

Mar 30, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE



=cut

my $coord_filter = defined $OPT{coord_filtername}?$OPT{coord_filtername}:'filter_dgv';

my $runid = $OPT{runid};
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $var_xml = modules::VariantXML->new(-outdir=>$run_dir . '/conf');

my $pipe_config = modules::Pipeline::get_pipe_conf();
my $source_name = modules::Pipeline::get_source_name($sample_name);
my $source_type = modules::Pipeline::get_source_type(-source_name=>$source_name);
my $max_distance = $pipe_config->read($source_type,'cutoffs','sv_max_distance');

my $coord_file_base = $OPT{coord_file};
my $chromosome = defined $OPT{chr}?$OPT{chr}:'';

my $last_chr = 0;
if (defined $OPT{chr} && $OPT{chr} eq 'Y') {
	$last_chr = 1;
}

my $sv_type = defined $OPT{sv_type}?$OPT{sv_type}:'del,dup,inv';

my $overlap_args;

my $overlap = modules::Overlap->new();

#Need to create tmp files for overlap mo
my $outfile_match_base = $OPT{overlap_outfile};

my $ref_file_base =  $OPT{ref_file};

#Each file requires regex
if ($coord_file_base !~ /SVTYPE/ || $ref_file_base !~ /SVCALLER/ || $outfile_match_base !~ /SVCALLER/) {
	modules::Exception->throw("ERROR: -coord_file must contain regex SVTYPE, -ref_file file and -overlap_outfile must contain SVCALLER");
}


my %sv_data;
my $pipe_conf = modules::Pipeline::get_pipe_conf();
my @sv_callers = split(",",$pipe_conf->read('common','sv_callers'));
my @sv_types = split(",",$sv_type);

for my $sv_caller ( @sv_callers ) {
    for my $sv_type ( @sv_types ) {
        (my $ref_file = $ref_file_base) =~ s/SVCALLER/$sv_caller/;
        (my $outfile_match = $outfile_match_base) =~ s/SVCALLER/$sv_caller/;
        (my $coord_file = $coord_file_base) =~ s/SVTYPE/$sv_type/;
        
        #Check the regexes worked
        if ( !-e $ref_file ) {
        	modules::Exception->throw("File $ref_file doesn't exist");	
        }
        if ( !-e $coord_file ) {
        	modules::Exception->throw("File $coord_file doesn't exist");	
        }
        
		my ($overlap_match,$match_count,$overlap_nomatch,$nomatch_count);
		#Build up the arguments for the overlap call
		my %args = (
					-ref=>$ref_file,
					-coord=>$coord_file,
					-silent=>1, #Don't need the file
					-max_distance=>$max_distance #Max distance to permit for breakpoints
					);
		
		
		
		#Call for matching
		($overlap_match,$match_count) = $overlap->overlap(%args);
		
		#Now call for nonmatching; need to change some of the arguments
		(my $outfile_nomatch = $outfile_match) =~ s/match/nomatch/;
		$args{-fail} = 1;
		($overlap_nomatch,$nomatch_count) = $overlap->overlap(%args);
		
		my @sorted_svs = ();
		my $file_name = $sample_name . '_' . $runid . '.'.$sv_caller.'.sv.xml';
#		if ( !-e $file_name ) {
#			modules::Exception->throw("File $file_name doesn't exist");	
#		}
		$var_xml->load_xml_keys(-file_name=>$file_name,-sv=>1);
		my $svs = $var_xml->search_xml_keys(-var_type=>'sv');
		my $count = @{$svs};
		if ($svs && @{$svs}) {
			@sorted_svs = sort {$a->{chr} cmp $b->{chr} || $a->{start_coord} <=> $b->{start_coord}} @{$svs};
		}

		#Here we overwrite the default match/nomatch files as coordinate matching isn't enough with indel
		open(MATCH,">$outfile_match.$sv_type") || modules::Exception->throw("Can't open file to write $outfile_match.$sv_type\n");
		open(NOMATCH,">$outfile_nomatch.$sv_type") || modules::Exception->throw("Can't open file to write $outfile_nomatch.$sv_type\n");
		
		for my $sv (@sorted_svs) {
			my $newsv_chr = my $newsv_start = my $newsv_end  = my $newsv_type = my $newsv_caller;
			$newsv_chr = $sv->{chr};
			$newsv_start = $sv->{start_coord};
			$newsv_end = $sv->{end_coord};
			$newsv_type = $sv->{sv_type};
			$newsv_type =~ s/sv//;
			$newsv_caller = $sv->{sv_caller};
			$newsv_caller =~ s/sv//;
			
			next unless $newsv_type eq $sv_type; #Don't 'match' different sv types
			next unless $newsv_caller eq $sv_caller; #Don't match different SV callers
			my $sv_key = $newsv_chr . ':' . $newsv_start;
			
#			print "caller $sv_caller type $sv_type count $count\n";
#			print Dumper $sv;
#			print Dumper $overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end};
#			print Dumper $overlap_match->{FAIL}{$newsv_chr}{$newsv_start}{$newsv_end};
#			print "\n\n";
			
			if (exists $overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end}) {
				my $match_found = 0;
				for my $sv_entry ( keys %{$overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end}} ) {
					my %matches = ();
					my $match_str;
					for my $match_entry (sort @{$overlap_match->{PASS}{$newsv_chr}{$newsv_start}{$newsv_end}{$sv_entry}}) {
						my @fields = split('\^\^\^',$sv_entry);
						my (undef,$match_sv_type) = split('=',$fields[1]);
			    		next unless $newsv_type eq $match_sv_type; #Don't 'match' different sv types
						$match_found = 1;
						$matches{$match_entry}++;
					}
					
					if (keys %matches) {
						my $match_str = join(",",keys %matches);
						print MATCH join("\t",
											$newsv_chr,
											$newsv_start,
											$newsv_end,
											$sv_entry,
											$match_str
											) . "\n";
											
						$sv_data{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{pass} = 0;
						(my $match_formatted = $match_str) =~ s/\^\^\^/;/g;
						$sv_data{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{sv_filter_dgv_string} = $match_formatted;
						
					} else {
						print NOMATCH join("\t",
										$newsv_chr,
										$newsv_start,
										$newsv_end,
										$sv_entry,
										'NO_MATCH'
										) . "\n";
						$sv_data{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{sv_filter_dgv_string} = 'novel';
						$sv_data{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{pass} = 1;
						
					}					
					
				}
				
				
				
			} elsif (exists $overlap_nomatch->{FAIL}{$newsv_chr}{$newsv_start}{$newsv_end}) {
				my ($sv_entry) = keys %{$overlap_nomatch->{FAIL}{$newsv_chr}{$newsv_start}{$newsv_end}};
				print NOMATCH join("\t",
									$newsv_chr,
									$newsv_start,
									$newsv_end,
									$sv_entry,
									'NO_MATCH'
									) . "\n";
				$sv_data{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{sv_filter_dgv_string} = 'novel';
				$sv_data{$sv_key}{$newsv_end}{$sv_caller}{$newsv_type}{pass} = 1;
	    	} else {
	    		my $sv_str = $newsv_chr.":".$newsv_start."-".$newsv_end;
				modules::Exception->throw("ERROR: Can't find snp $sv_str in the fail or pass list\n");
	    	}
			
		}
		
    }
    
}

my $file_name = $sample_name . '_' . $runid . '.'.$coord_filter.'.sv.xml';
$var_xml->create_var_xml(-file_name=>$file_name,-data=>\%sv_data, -chr=>'all',-sv=>1);
my $full_xml_file = $run_dir . '/conf/'.$file_name;
$var_xml->split_xml(-file_name=>$full_xml_file);





