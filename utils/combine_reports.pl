#!/usr/bin/perl
use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "report1=s",
	   "report2=s",
	   "combine_file=s"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{report1} || !$OPT{report2});



=pod

=head1 SYNOPSIS

script_name [options]

Required flags: NONE

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

script_name.pl -> One line description

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

sample -f file

=cut


my %data = ();

my $report1 = $OPT{report1};
my $report2 = $OPT{report2};
my @reports = ($report1,$report2);
my %headers = ();
my %key_count = ();


my $combine = defined $OPT{combine_file}?$OPT{combine_file}:"combined.tsv"; 

if ( !-e $report1 || !-e $report2) {
	modules::Exception->throw("File $report1 or $report2 doesn't exist");
}

open(COMBINE,">$combine") || modules::Exception->throw("Can't open file $combine\n");

for my $file (@reports) {
	open(FILE,"$file") || modules::Exception->throw("Can't open file $file\n");
	
	while (<FILE>) {
		chomp;
		
		my @fields = split("\t");
		my $count = 0;
		if (/^chr/) {
			
			for my $field (@fields) {
				if (exists $headers{$count} && $headers{$count} != $field) {
					modules::Exception->throw("ERROR: Headers must be the same");
				}
				$headers{$count} = $field;
				$count++;
			}
			next;
		}
		next if /^$/;
		my $key = "$fields[0]:$fields[1]:$fields[2]:$fields[13]";
		
		$key_count{$key}++;
		
		for my $field (@fields) {
			$data{$key}{$count}{$file} = $field;
			$count++;
		}
		

		#$data{$key}{$file} = "$fields[7]:$fields[8]:$fields[9]";
	}
	
}

my %multi_entries = ();
for my $key ( keys %key_count ) {
	$multi_entries{$key_count{$key}}{$key} = $data{$key};
}

open(COMBINE,">$combine") || modules::Exception->throw("Can't open file $combine\n");

my @headers = ();
for my $count (sort {$a<=>$b} keys %headers) {
	push @headers, $headers{$count};
}

print COMBINE join ("\t",
					"Algorithm(s)",
					@headers
					) . "\n";


my %total_counts = ();

for my $count ( sort {$b<=>$a} keys %multi_entries ) {
	if ($count > 1) {
		print COMBINE "Entries in both reports $report1 and $report2\n";
	} else {
		print COMBINE "\n\nEntries unique to either $report1 or $report2\n";
	}
	
	
	
    for my $key ( sort {my ($a_chr,$a_coord) = $a =~ /([0-9XYM]+):(\d+)/; my ($b_chr,$b_coord) = $b =~ /([0-9XYM]+):(\d+)/; $a_chr cmp $b_chr || $a_coord <=> $b_coord} keys %{$multi_entries{$count}} ) {
		my @col_data = ();
    	if ($count > 1) {
    		push @col_data, "Strelka_Mutect";
    		$total_counts{'Strelka_Mutect'}++;
    	} 
    	
		for my $col_count (sort {$a<=>$b} keys %{$multi_entries{$key_count{$key}}{$key}}) {
	    	my @val = ();
	    	
	    	
			for my $file (sort keys %{$multi_entries{$key_count{$key}}{$key}{$col_count}}) {
		    	if ($col_count == 0 && $count == 1)  {
		    		my @algs = split("_",$file);
    				push @col_data, $algs[1];
    				$total_counts{$algs[1]}++;
	    		}
				push @val, $multi_entries{$key_count{$key}}{$key}{$col_count}{$file};
				
			}		
			
			if ($val[0] eq $val[1]) {
				push @col_data,$val[0];
			} else {
				push @col_data, join(":",@val);
				
			}
		}
		
		print COMBINE join("\t",
							@col_data
							) . "\n";
		
		
	}
}

print COMBINE "\nSUMMARY:\n";

for my $alg (sort keys %total_counts) {
	print COMBINE $alg.":\t".$total_counts{$alg}."\n";
}

#print Dumper \%headers;

#print Dumper \%multi_entries;



