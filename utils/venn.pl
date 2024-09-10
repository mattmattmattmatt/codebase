#! /usr/bin/perl -w

use strict;
use Getopt::Long;
use modules::Exception;
use modules::Utils;
use Pod::Usage;
use Data::Dumper;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "title=s",
	   "files=s",
	   "set_columns=s",
	   "delim=s",
	   "s1=s",
	   "s2=s",
	   "s3=s",
	   "debug"
	   );
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{files});
	   
=pod

=head1 SYNOPSIS

venn.pl -files comma_delim_list_of_files -delim deliminator(default="\t") -set_columns columns_to_use_for_unique_set_id(default=0,1) -s1 s1_name(default=IGL_name) -s2 s2_name(default=IGL_name) -s3 s3_name(default=IGL_name)-title title_for_graph(default=Venn.giff)

Required flags: -files

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

venn.pl -> generate 2 or 3 set venn diagram from overlapping files (requires R library VennDiagram)

=head1 DESCRIPTION

May 02, 2011

a script that ...

=head1 AUTHOR

Matt Field

=head1 EXAMPLE

venn.pl -title "Overlapping SNP Calls of Technical Replicates" -files /home/matt/work/docs/papers/anu/exome_analysis/IGL00133_19_filter_basic_match,/home/matt/work/docs/papers/anu/exome_analysis/IGL00134_20_filter_basic_match,/home/matt/work/docs/papers/anu/exome_analysis/IGL00135_21_filter_basic_match

=cut

my @files = split(",",$OPT{files});

if (@files != 2 && @files != 3 && @files != 4) {
	modules::Exception->throw("Need to pass in 2-4 files");
}

my $title = defined $OPT{title}?$OPT{title}:'Venn Diagram';
my @key_columns = defined $OPT{set_columns}?split(',',$OPT{set_columns}):qw(0 1);
my $delim = defined $OPT{delim}?$OPT{delim}:"\t";	

my $rand = int(rand(10000000));

print Dumper \@key_columns;

my $venn_r_file = "venn$rand.r";

print "R $venn_r_file\n";



open(R,">$venn_r_file") || modules::Exception->throw("Can't open venn file");;
print R "library(VennDiagram)\n";
my $set_count = 1;
my $list_str;
foreach my $file (@files) {
    my $set_file = "s$set_count";
    if ($OPT{s1} && $set_file eq 's1') {
    	$list_str .= "$OPT{s1} = $set_file,"
    } elsif ($OPT{s2} && $set_file eq 's2') {
    	$list_str .= "$OPT{s2} = $set_file,"
    } elsif ($OPT{s3} && $set_file eq 's3') {
    	$list_str .= "$OPT{s3} = $set_file,"
    }
    
    open(FILE,$file) || modules::Exception->throw("Can't open $file");
    my $r_string = "s$set_count<-c(";

    my $count = 0;
    while (<FILE>) {
    	chomp;
		my @cols = split($delim,$_);
		$r_string .= "\"";
		if (@key_columns == 1) {
			$r_string .= $cols[0];
		} else {
			for my $key_entry ( @key_columns ) {
			    $r_string .= $cols[$key_entry] . ':';
			}
			$r_string =~ s/:$//;
			
		}

		$r_string .= "\",";
		
		$count++;
		#last if $count > 25;
    }
    $r_string =~ s/\,$/\)/;
    print R "$r_string\n";
    $set_count++;
}

print R "jpeg(\"$title.jpeg\")\n";
$list_str =~ s/,$//;

if (@files == 4) {
	print R "venn.diagram(list($list_str), \"$title.png\",  main = \"$title\", col=\"transparent\",fill = c(\"cornflowerblue\", \"green\", \"yellow\",\"red\"),scaled = TRUE, imagetype =\"png\")\n";
} elsif (@files == 3) {
	print R "venn.diagram(list($list_str), \"$title.png\",  main = \"$title\", col=\"transparent\",fill = c(\"cornflowerblue\", \"green\", \"yellow\"),scaled = TRUE, imagetype =\"png\")\n";
} else {
	print R "venn.diagram(list($list_str), \"$title.png\",  main = \"$title\", col=\"transparent\",fill = c(\"cornflowerblue\", \"green\"), imagetype =\"png\")\n";
}


system("R --save < $venn_r_file 1>/dev/null");
#system("rm $venn_r_file");
#system("mv RPlots.pdf ")


