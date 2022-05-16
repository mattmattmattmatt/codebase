#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use vars qw(%OPT);
use modules::Exception;

GetOptions(\%OPT,
	   "samplenamesFile=s");

my @genofiles;

my $sample_names_dir="./sample_names";
my $outdir = "./human_snp_validation";
my $genodir = "./human_snp_validation";

open(SAMPLEFILE,$sample_names_dir."\/".$OPT{samplenamesFile}) || modules::Exception->throw("Can't open file $OPT{samplenamesFile}\n");
while(<SAMPLEFILE>){
  chomp;
  my $sample_file = "$_.genotype.tsv";
  push @genofiles, $sample_file;
}
close(SAMPLEFILE);

my %geno_data;
my $outfile = $outdir."\/".$OPT{samplenamesFile};
$outfile =~ s/sample_names.txt/all.genotype.tsv/;
my @header = ("chr","coord");
my %snp_loc;
my $count_file=0;	#get chr and coord from the first file

for my $geno_file (@genofiles){
  my $sample_name;
  $count_file++;
  open(GENOFILE,$genodir."\/".$geno_file) || modules::Exception->throw("Can't open file $geno_file");
  while(<GENOFILE>){
      chomp;
      if($.==1){
	  my @header_tmp = split("\t");
	  $sample_name = $header_tmp[2];
	  push @header, $sample_name;
      }else{
	  my ($chr,$coord,$geno) = split("\t");
	  $geno_data{$sample_name}{$chr}{$coord} = $geno;
	  if($count_file==1){
	    $snp_loc{$chr}{$coord}=1;
	  }
      }
  }
  close(GENOFILE);
}

open(OUTFILE,">$outfile") || modules::Exception->throw("Can't open file to write $outfile\n");

print OUTFILE join("\t", @header)."\n";		#header

for my $chr (sort {$a<=>$b} keys %snp_loc){
	for my $coord (sort {$a<=>$b} keys %{$snp_loc{$chr}}){
	  
		my $geno_line = join("\t",$chr,$coord);
		for my $sample (@header[2..$#header]){
		    $geno_line = $geno_line."\t".$geno_data{$sample}{$chr}{$coord};
		}
			print OUTFILE "$geno_line\n";
	}
	
}

close OUTFILE;

#print Dumper \@genofiles;
#print Dumper \%geno_data;
#print Dumper \%snp_loc;
#print Dumper \$outfile;
