#!/usr/bin/perl -wall

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
	   "match_nums=s",
	   "nonmatch_nums=s",
	   "key=s",
	   "entries=s",
	   "debug"
	   );
	   
pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{match_nums} || !$OPT{nonmatch_nums} || !$OPT{entries});
	   
=pod

=head1 SYNOPSIS

histogram.pl -title title_for_graph(default=histogram.tiff) -entries entry_name_for_X_axis_labels -match_nums comma_delim_match_nums -nonmatch_nums comma_delim_nonmatch_nums -key comma_delim_key_for_match_and_nonmatch(default=match/nonmatch)

Required flags: -match_num -nonmatch_num -entries

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

histogram.pl -> generate a match/nonmatch histogram

=head1 DESCRIPTION

May 02, 2011

a script that ...

=head1 AUTHOR

Matt Field
Simple hist code:
> histdata <- read.table("quality_scores_pileup")
> histcol = histdata[,1]
> histcol <- as.numeric(histcol)
> hist(histcol,main="Quality scores from pileup",xlab="Quality score",ylab="Count")


=head1 EXAMPLE

./histogram.pl -match_num 1765542758,1736203305,1434516888,1391979506 -nonmatch_num 83244831,88402858,56512129,64476680 -entries RAM_A15,RAM_A15_LCL,AGRF_A15,AGRF_A15_LCL -title AGRF_and_RAMACIOTTI_READ_DEPTH -key Aligned,Non-Aligned 

=cut

my $title = defined $OPT{title}?$OPT{title}:'Match and Non-match Histogram';
my $key = defined $OPT{key}?$OPT{key}:"Match,Non-match";
my $match_str = $OPT{match_nums};
my $non_match_str = $OPT{nonmatch_nums};
my $entry = $OPT{entries};



unless ($match_str =~ /^[0-9,]+$/ && $non_match_str =~ /^[0-9,]+/) {
	modules::Exception->throw("ERROR: Match and non-match string must be comma delim integers");
}

my @match_nums = split(",",$match_str);
my @non_match_nums = split(",",$non_match_str);
my @entries = split(",",$entry);
my @quoted_keys;
for my $key (split(",",$key)) {
	push @quoted_keys, '"'. $key . '"';
}
my $quoted_keys = join(",",@quoted_keys);

unless (@match_nums == @non_match_nums && @match_nums == @entries) {
	modules::Exception->throw("ERROR: Require the same number of match_nums, non_match_nums, and entries");
}

my $frame_data;
my $entry_labels;

for (my $i=0; $i < @match_nums; $i++) {
	$frame_data .= 'c("'.$match_nums[$i].'","'.$non_match_nums[$i].'"),';
	$entry_labels .= '"'.$entries[$i].'",';
}

$frame_data =~ s/,$//;
$entry_labels =~ s/,$//;



my $hist_r_file = 'histogram.r';

open(RCODE,">$hist_r_file") || modules::Exception->throw("Can't open hist file");
print RCODE "hist_data = data.frame($frame_data)\n";
print RCODE "par(xpd=T, mar=par()\$mar+c(0,0,0,4))\n";
print RCODE "barplot(as.matrix(hist_data), main=\"$title\", ylab=\"Total Count\", col=c(\"greenyellow\",\"red3\"), space=0.1, cex.axis=0.7, las=1, names.arg=c($entry_labels), cex=0.8)\n";
print RCODE "legend('topleft', c($quoted_keys), cex=0.7, fill=c(\"greenyellow\",\"red3\"))\n";
close RCODE;

(my $title_filename = $title) =~ s/ /_/g;

system("R --save < $hist_r_file");
system("mv Rplots.pdf $title_filename.pdf");
#system("rm histogram.r");





