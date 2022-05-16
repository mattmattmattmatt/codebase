#! /usr/bin/perl -w 
use strict;
use Data::Dumper;
use Getopt::Long;
use modules::SystemCall;
use modules::Exception;
use modules::ConfigXML;
use modules::Pipeline;
use modules::Adaptors::Read_File;
use modules::Adaptors::Lane;
use modules::Adaptors::Sample;
use Pod::Usage;
use PDF::API2::Simple; 
use PDF::Table; 
use modules::Utils;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
     	"runid=i"
	   );
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{runid});

	   
=pod

=head1 SYNOPSIS

nata_report.pl -runid runid [options]

Required flags: -runid

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

nata_report.pl -> Generate summary NATA report

=head1 DESCRIPTION

May 03, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

nata_report.pl 

=cut

my $runid = $OPT{runid};
my $sample_name = modules::Pipeline::get_sample_name(-run_id=>$runid);
my $run_dir = modules::Pipeline::get_run_dir(-run_id=>$runid,-sample_name=>$sample_name);
my $genome_size = 2897310462;
my $exon_size = 33363981;
my $summary_dir  = $run_dir .'/summary';
my $log_dir = $run_dir . '/log';
my $sequence_type = modules::Pipeline::get_sequence_type(-sample_name=>$sample_name);
my $clus_config = modules::Pipeline::get_cluster_conf();
my $pipe_config = modules::Pipeline::get_pipe_conf();
my $source_type = modules::Pipeline::get_source_type(-sample_name=>$sample_name);
my ($sample_obj) = modules::Adaptors::Sample->search(sample_name=>$sample_name);

my $date = modules::Utils::GetTime();

my %data = ();
#Load up some basic info
$data{sample}{external_sample_name} = $sample_obj->external_sample_name;
$data{sample}{total_read_pairs} = $sample_obj->total_lanes;
$data{sample}{sequence_type} = $sequence_type;
$data{path}{rundir} = $run_dir;
$data{sample}{sample_type} = $sample_obj->sample_type;



opendir(LOG,$log_dir) || modules::Exception->throw("Can't open directory $log_dir");
my ($log_file) = grep {/.log/} readdir LOG;
open(LOGFILE,"$log_dir/$log_file") || modules::Exception->throw("Can't open file $log_file\n");
my $binary_flag = 0;

while (<LOGFILE>) {
	chomp;
	
	if (/^$/) {
    	$binary_flag = 0;
    }
	
	if (/Svn revision: (\d+)/) {
		$data{software}{svn} = 'revision'.$1;
	}
	
	if (/Starting pipeline run\s+(\S+)\s+at\s+(\S+)\s+on\s+(\S+)/) {
		$data{analysis}{start_time} = $2;
		$data{analysis}{host} = $3;
		($data{sample}{sample_name},$data{analysis}{run_number}) = $1 =~ /(.*)_(\d+)$/;
	}
	
	if (/^SVNDIR:\s+(\S+)/) {
		$data{path}{svn} = $1;
	}
	
	if ($binary_flag) {
		my @fields = split(':',$_);
		$fields[1] =~ s/^ //g;
		$data{software}{$fields[0]} = $fields[1];
	}

	if (/^BINARIES/) {
		$binary_flag = 1;
	}

	
}

#Get steps run
my $xml = $run_dir .'/conf/'. $data{sample}{sample_name} .'_' .$data{analysis}{run_number} .'.conf.xml';

if ( !-e $xml ) {
	modules::Exception->throw("File $xml doesn't exist");	
}

my $xml_config = modules::ConfigXML->new($xml);
my $steps = $xml_config->read('steps_order', 'step');

$data{analysis}{steps_run} = $steps;

#Get annotation versions
open(XML,"$xml") || modules::Exception->throw("Can't open file $xml\n");

while (<XML>) {
	if (/ENTITY ref\s+\"(\S+)\"/) {
		$data{annotation}{ref_genome} = $1;
	}	
	if (/ENTITY dbsnpVersion\s+\"(\S+)\"/) {
		$data{annotation}{dbsnp_version} = $1;
	}	
	if (/ENTITY ExACVersion\s+\"(\S+)\"/) {
		$data{annotation}{exac_version} = $1;
	}
	if (/ENTITY GnomADVersion\s+\"(\S+)\"/) {
		$data{annotation}{gnomad_version} = $1;
	}
	if (/ENTITY ClinvarVersion\s+\"(\S+)\"/) {
		$data{annotation}{clinvar_version} = $1;
	}
	if (/ENTITY exonVersion\s+\"(\S+)\"/) {
		$data{annotation}{exon_version} = $1;
	}
	if (/ENTITY cosmicVersion\s+\"(\S+)\"/) {
		$data{annotation}{cosmic_version} = $1;
	}
	if (/ENTITY geneVersion\s+\"(\S+)\"/) {
		$data{annotation}{gene_version} = $1;
	}
	if (/ENTITY dgvVersion\s+\"(\S+)\"/) {
		$data{annotation}{dgv_version} = $1;
	}
}

#Get few from pipe.xml that aren't in the sample xml
my $annotations = $pipe_config->read($source_type,'annotation_version');


for my $annotation (keys %{$annotations}) {
	if (!exists $data{annotation}{$annotation}) {
		$data{annotation}{$annotation} = $annotations->{$annotation};
	} 
}

#$data{annotation}{ensembl} = $pipe_config->read($source_type,'annotation_version','ensembl_version');



opendir(SUM,$summary_dir) || modules::Exception->throw("Can't open directory $summary_dir");
my @files = grep {/Report/} readdir SUM;

#Variant summaries
for my $file ( @files ) {
    if ($file =~ /read/) {
    	open(READ,"$summary_dir/$file") || modules::Exception->throw("Can't open file $file\n");
    	my $total = my $aligned = my $exon = 0;
    	while (<READ>) {
    		chomp;
    		if (/^total: (\d+)/) {
				$total = $1;
			} elsif (/^aligned: (\d+)/) {
				$aligned = $1;
			} elsif (/^unaligned: (\d+)/) {
				$data{read_pairs}{unaligned} = int($1/2);
			} elsif (/^exon_aligned: (\d+)/) {
				$exon = $1
			} 
    	}
    	
    	my $percent_aligned = sprintf("%.2f",$aligned/$total * 100);
    	my $percent_exon_aligned = sprintf("%.2f",$exon/$total * 100);
    	
    	$data{read_pairs}{total} = int($total/2);
    	$data{read_pairs}{aligned} = int($aligned/2);
    	
    	$data{read_pairs}{percent_ref_aligned} = int($aligned/2) . '/' . int($total/2) . ' ('. $percent_aligned. '%)';
    	$data{read_pairs}{percent_exon_aligned} = int($exon/2) . '/' . int($total/2) . ' ('. $percent_exon_aligned. '%)';

		#Get sample read file to find read length
		my ($lane_obj) = modules::Adaptors::Lane->search(sample_id=>$sample_obj->id);
		my ($read_file_obj) = modules::Adaptors::Read_File->search(lane_id=>$lane_obj->id);
		
		my $read_file_full = $read_file_obj->read_directory .'/' . 	$read_file_obj->file_name;		
    	
    	$data{path}{readdir} = $read_file_obj->read_directory;
    	
    	
    	if (!-e $read_file_full) {
    		modules::Exception->throw("ERROR: Can't open read file $read_file_full");
    	}
    	
    	my $command = "zcat $read_file_full | head -2 | tail -n 1 | wc -c"; #Gives error but doesn't matter 
    	my $read_length = `$command`;
    	$read_length--; #Account for \n
    	$data{read_pairs}{read_length} = $read_length;
    	
    	if ($sequence_type eq 'genome') {
			$data{read_pairs}{coverage} =  sprintf ("%.2f",$read_length * $aligned / $genome_size) . 'X';
		} elsif ($sequence_type eq 'exome') {
			$data{read_pairs}{coverage} =  sprintf ("%.2f",$read_length * $exon / $exon_size) . 'X';
		}
    	
    	
    	
    } elsif ($file =~ /exon/) {
    	my %gene_missing = ();
    	my $total_missing = my $total_exon = 0;
    	open(EXON,"$summary_dir/$file") || modules::Exception->throw("Can't open file $file\n");
    	while (<EXON>) {
    		chomp;
    		my @fields = split("\t");
    		if (/UNCOVERED/) {
    			my ($missing,$total) = $_ =~ /(\d+)\/(\d+)/;
    			$total_missing += $missing;
    			$total_exon += $total;
    		} elsif (/^[0-9XY]/) {
	    		$gene_missing{$fields[4]}++;
    		}
    		
    	}
    	
    	my $gene_cover_str = $total_missing . '/' . $total_exon . ' ('. sprintf("%.2f",($total_missing/$total_exon) *100) . '%)';
    	$data{genes}{Genes_missing_exon_bases} = keys %gene_missing;
		$data{genes}{Total_low_cover_exon_bases} = $gene_cover_str;
		
		
		
    	
    } elsif ($file =~ /snv/) {
    	&Parse_Reports($file,'snv');
    } elsif ($file =~ /indel/) {
    	&Parse_Reports($file,'indel');
    } elsif ($file =~ /delly/) {
    	&Parse_Reports($file,'sv_delly');
    } elsif ($file =~ /lumpy/) {
    	&Parse_Reports($file,'sv_lumpy');
    }
    #TODO: Add combined SV reports
}


#print Dumper \%data;


my $sample_data = [['SAMPLE INFO',''],
					['External Sample Name',$data{sample}{external_sample_name}],
					['Internal Sample Name',$data{sample}{sample_name}],
					['Sample Type',$data{sample}{sample_type}],
					['Sequence Lanes',$data{sample}{total_read_pairs}],
					['Sequence Type',$data{sample}{sequence_type}]
					];

my $read_data = [['READ INFO',''],
					['Read Length',$data{read_pairs}{read_length}],
					['Total Read Pairs',$data{read_pairs}{total}],
					['Aligned Read Pairs',$data{read_pairs}{aligned}],
					['Unaligned Read Pairs',$data{read_pairs}{unaligned}],
					['Percent Genome Aligned',$data{read_pairs}{percent_ref_aligned}],
					['Percent Exon Aligned',$data{read_pairs}{percent_exon_aligned}],
					['Est. Coverage',$data{read_pairs}{coverage}],
					['Genes With Low/No Cover',$data{genes}{Genes_missing_exon_bases}],
					['Total Low/No Cover Exon Bases',$data{genes}{Total_low_cover_exon_bases}],
					];
					

my $steps_str = join("\n",@{$data{analysis}{steps_run}});

my $analysis_data  = [['ANALYSIS INFO',''],
						['Run Number',$data{analysis}{run_number}],
						['Start Time',$data{analysis}{start_time}],
						['Compute Host',$data{analysis}{host}],
						['Read Directory',$data{path}{readdir}],
						['SVN Directory',$data{path}{svn}],
						['Run Directory',$data{path}{rundir}],			
						['Analysis Steps',$steps_str]
						];	
						
						
my $software_data = [['SOFTWARE INFO',''],
						['SVN',$data{software}{svn}],
						['BWA',$data{software}{bwa}],
						['JAVA',$data{software}{javabin}],
						['Picard',$data{software}{picard}],
						['GATK',$data{software}{gatk}],
						['SAMTools',$data{software}{samtools}],
						['BCFTools',$data{software}{bcftools}],
						['Tabix',$data{software}{tabix}],
						['Variant Effect Predictor',$data{software}{variant_predictor}],
						['Lumpy',$data{software}{lumpy}],
						['Delly',$data{software}{delly}]
						];

my $annotation_data = [['ANNOTATION INFO',''],
						['Ensembl',$data{annotation}{ensembl_version}],
						['Reference Genome',$data{annotation}{ref_genome}],
						['DBSNP',$data{annotation}{dbsnp_version}],
						['Exac',$data{annotation}{exac_version}],
						['Gnomad',$data{annotation}{gnomad_version}],
						['Exon version',$data{annotation}{exon_version}],
						['Gene version',$data{annotation}{gene_version}],
						['DGV',$data{annotation}{dgv_version}],
						['Clinvar',$data{annotation}{clinvar_version}],
						['miRNA',$data{annotation}{mirna_version}],
						['aB_FOXP3_regulatory_version',$data{annotation}{aB_FOXP3_regulatory_version}],
						['Cosmic',$data{annotation}{cosmic_version}]
						];


my $snv_data = [['SNV INFO',''],
				['Input Set',join("\n",@{$data{variants}{snv}{input_filters}})],
				['Total SNV Calls',$data{variants}{snv}{total}],
				['Exon SNVs',$data{variants}{snv}{filter_exon}],
				['Splice Site (+-10bp) SNVs',$data{variants}{snv}{splice}],
				['Missense SNV Calls',$data{variants}{snv}{missense}],
				['Synonomous SNV Calls',$data{variants}{snv}{syn}],
				['Nonsense SNV Calls',$data{variants}{snv}{nonsense}],
				['Applied Filters',join("\n",@{$data{variants}{snv}{include_filters}})],
				['Passed SNV Calls',$data{variants}{snv}{passed}{total}],
				['Passed exon SNVs',$data{variants}{snv}{passed}{Exonic}],
				['Passed splice SNVs',$data{variants}{snv}{passed}{'Splice-site'}],
				['Passed median depth SNVs',$data{variants}{snv}{passed}{'Median Depth'}],
				['Passed median qual SNVs',$data{variants}{snv}{passed}{'Median Quality'}]				
				];

my $indel_data = [['Small Indel INFO',''],
				['Input Set',join("\n",@{$data{variants}{snv}{input_filters}})],
				['Total Indel Calls',$data{variants}{indel}{total}],
				['Exon Indels',$data{variants}{indel}{filter_exon}],
				['Splice Site (+-10bp) SNVs',$data{variants}{indel}{splice}],
				['Insertion',$data{variants}{indel}{ins}],
				['Deletions',$data{variants}{indel}{del}],
				['Applied Filters',join("\n",@{$data{variants}{indel}{include_filters}})],
				['Passed SNV Calls',$data{variants}{indel}{passed}{total}],
				['Passed exon SNVs',$data{variants}{indel}{passed}{Exonic}],
				['Passed splice SNVs',$data{variants}{indel}{passed}{'Splice-site'}],
				];

my $sv_data = [['Lumpy SV/CNV INFO',''],
				['Input Set',join("\n",@{$data{variants}{sv_lumpy}{input_filters}})],
				['Total SV Calls',$data{variants}{sv_lumpy}{total}],
				['Deletion Calls',$data{variants}{sv_lumpy}{del}],
				['Duplication Calls',$data{variants}{sv_lumpy}{dup}],
				['Inversion Calls',$data{variants}{sv_lumpy}{inv}],
				['Translocation Calls',$data{variants}{sv_lumpy}{tra}],
				['Applied Filters',join("\n",@{$data{variants}{sv_lumpy}{include_filters}})],
				['Passed SV Calls',$data{variants}{sv_lumpy}{passed}{total}]
				];


my $page_num = 1; 

my $pdf_file = $summary_dir . '/' . $data{sample}{sample_name}.'.pdf';
my $pdf = PDF::API2::Simple->new(  
                                  file => $pdf_file, 
		                          header => \&header, 
		                          footer => \&footer 
		                          ); 
my $pdftable = new PDF::Table;  
$pdf->add_font('VerdanaBold'); 
$pdf->add_font('Verdana'); 

$pdf->add_page();

&PDF_Table($sample_data);
&PDF_Table($read_data);
&PDF_Table($analysis_data);
&PDF_Table($software_data);
&PDF_Table($annotation_data);
&PDF_Table($snv_data);
&PDF_Table($indel_data);
&PDF_Table($sv_data);
      
$pdf->save(); 


sub PDF_Table {
	
	my $data = shift;
	
	# build the table layout 
	$pdftable->table( 
	     # required params 
	     $pdf->pdf, 
	     $pdf->current_page, 
	     $data, 
	     x => $pdf->margin_left, 
	     w => $pdf->effective_width, 
	     start_y => $pdf->margin_bottom + 500, 
	     next_y  => 750, 
	     start_h => 500, 
	     next_h  => 720, 
	     # some optional params 
	     new_page_func  => $pdf->add_page(), 
	     header_props          => {  
	                font       => $pdf->current_font, 
	        font_size  => 9, 
	        font_color => '#000000', 
	        bg_color   => '#FF9046', 
	        repeat     => 1,  
	         }, 
	     font => $pdf->current_font, 
	     font_size      => 9, 
	     padding => 5, 
	     padding_right => 10, 
	     #background_color_odd  => "gray", 
	     #background_color_even => "lightblue", 
	      
	      #cell background color for even rows 
	  );  
	################################################3 
}  
  
  
  
sub header { 
    my $strokecolor = $pdf->strokecolor; 
     my $y = $pdf->height - 20; 
    my $x = ($pdf->width / 2);  
    $pdf->stroke_color( '#555555' ); 
  
    $pdf->next_line; 
    $pdf->set_font( 'VerdanaBold' ); 
    $pdf->text( 'NATA Report For '.$data{sample}{external_sample_name}, 
            x => $x, 
                   y => $y, 
                   font_size => 10, 
            align =>'center' ); 
  
    $pdf->y( $pdf->y - 5 ); 
     
  
     
} 
  
sub footer { 
 
    my $strokecol = $pdf->strokecolor; 
    $pdf->stroke_color( '#000000' ); 
    $pdf->line(     x => $pdf->margin_left, 
            y => 40, 
            to_x => $pdf->effective_width, 
            to_y => 40, 
            stroke => 'on', 
            fill => 'off', 
            width => 1 ); 
    my $fillcolor = $pdf->fill_color; 
    my $font = $pdf->current_font; 
  
    $pdf->fill_color( '#000000' ); 
     
    $pdf->set_font( 'VerdanaBold' ); 
    $pdf->text( $date, 
        x => $pdf->margin_left, 
        y => 30, 
        font_size => 9, 
        align => 'left' ); 

}



sub Parse_Reports () {
	my ($file,$type) = @_;
	
	open(FILE,"$summary_dir/$file") || modules::Exception->throw("Can't open file $file\n");
    my $filter_flag = my $breakdown_flag = 0;
    while (<FILE>) {
    		chomp;
    		
    		if (/^$/) {
    			$filter_flag = $breakdown_flag = 0;
    		}
    		
    		
    		
    		if (/^Total.*called: (\d+)/) {
    			$data{variants}{$type}{total} = $1;
    			my @filters = ($_ =~ /(filter_[a-z]+)/g);
    			$data{variants}{$type}{input_filters} = \@filters; 
    		}

			if (/(\S+)\s+MATCH: (\d+ \([0-9\.\%]+\))/) {
				$data{variants}{$type}{$1} = $2;
			}
    		
    		if ($filter_flag) {
    			$_ =~ s/^ //g;
    			$_ =~ s/^\t//g;
    			push @{$data{variants}{$type}{include_filters}}, $_;
    		}

			if ($breakdown_flag) {
    			$_ =~ s/^ //g;
    			$_ =~ s/^\t//g;
    			my @fields = split(':',$_);
    			next if $fields[0] =~ /Gene/;
    			$fields[1] =~ s/^ //g;
    			
    			$data{variants}{$type}{passed}{$fields[0]} = $fields[1];
    		}

    		if (/^Filters/) {
    			$filter_flag = 1;
    		}
    		
    		if (/^Breakdown/) {
    			$breakdown_flag = 1;
    		}
    		
    		if (/^Total.*variants: (\d+)/) {
    			$data{variants}{$type}{passed}{total} = $1;
    		}
    		
    }
}





