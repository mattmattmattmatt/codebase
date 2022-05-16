#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd qw(abs_path getcwd);
use DBI;
use POSIX qw(strftime);

GetOptions(
	'input=s'      => \my $inpName,       # name of the input tsv file (sample <-> file_name)
	'output=s'     => \my $outName,       # name prefix for the two output files - <prefix>.qsub (file coping) and <prefix>.csv (sample_info.csv entries)
  'project=s'    => \my $prid,          # CACPI CPIEEU APOSLE CPINJ
  'facility=s'   => \my $seqc,          # Annoroad Macrogen
  'dir=s'        => \my $dir,           # directory with source fastq files
  'reads=s'      => \my $reads,         # directory where to create dirs with read files
  'checkfq'      => \my $fastqcheck,    # check existance of fastq files
  'overwrite'    => \my $overwrite,     # overwrite reads existsing dirs and files
  'familyfile=s' => \my $familyFile,    # read family relations from a file not CPI database, "~/annoroad_family_batch2.vbsv"
  'verbose'      => \my $verbose        # verbose output
) || die "ERROR: Invalid options passed to $0\n";

pod2usage(1) if (!defined $inpName || !defined $outName || !defined $prid || !defined $seqc);

=pod

=head1 SYNOPSIS

prepare_samples.pl 
  --input <sample_seqfile.tsv> required (file name with sample_id <TAB> fastq_name)
  --output <out_prefix> required (file name prefix for .csv and .qsub files)
  --project <project id> required (project id)
  --facility <seq. facility> required (facility id)
  --familyfile <family_file.txt> optional (if family file to be used instead of APF database) 
  --dir <source fastq dir> optional (default=.) 
  --reads <reads dir> optional (default=/short/u86/reads) 
  --checkfq optional (check if source fastq files are accessable)
  --exec optional (create read files directory structure) 
  --overwrite optional (overwrite read files directories if present)
  --verbose optional (be more vocal)
  [options]

=head1 OPTIONS

=head1 NAME

prepare_samples.pl -> Script to prepare samples for variantdb processesing

=head1 DESCRIPTION

Nov 21, 2017

a script that ...

=head1 AUTHOR

Marcin Adamski

=head1 EXAMPLE

prepare_samples.pl --input test.tsv --output test --project CCG --facility BRF

Example of a familyfile:
RequestId|PatientDatabaseId|IndividualName|Proband|Gender|Affected|FatherName|MotherName|FamilyId|Relationship
|SH301|SH301|Yes|Female|Yes|||7|proband
|SH302|SH302|No|Male|No|||7|father
|SH303|SH303|No|Female|No|||7|mother
|SH304|SH304|Yes|Female|Yes|||8|proband
|SH305|SH305|No|Male|No|||8|father
|SH306|SH306|No|Female|No|||8|mother
|SH307|SH307|Yes|Female|Yes|||9|proband
|SH308|SH308|No|Male|No|||9|father
|SH309|SH309|No|Female|No|||9|mother
|SH310|SH310|Yes|Female|Yes|||10|proband
|SH311|SH311|No|Male|No|||10|father
|SH312|SH312|No|Female|No|||10|mother
|SH313|SH313|Yes|Male|Yes|||11|proband
|SH314|SH314|No|Male|No|||11|father
|SH315|SH315|No|Female|No|||11|mother

=cut


$dir     = '.'                 if(!defined $dir);
$reads   = '/short/u86/reads'  if(!defined $reads);

#let's stop any ralative path dramas before they even start
$dir   = abs_path($dir);
$reads = abs_path($reads);

die "ERROR: can't access source directory $dir\n"  if(!-d $dir);
die "ERROR: can't access reads directory $reads\n" if(!-d $reads);

die "ERROR: need argument --input\n"    if(!defined $inpName);
die "ERROR: need argument --output\n"   if(!defined $outName);
die "ERROR: need argument --project\n"  if(!defined $prid);
die "ERROR: need argument --facility\n" if(!defined $seqc);

#vardb database details:
my $host = "130.56.244.153";
my $user = "ro_user";
my $pass = "QPwoEI92";  # read-only user and password
my $db   = "vardb_genome_2_3";
my $dbh  = DBI->connect ("DBI:mysql:database=$db:host=$host", $user, $pass) || die "ERROR: Couldn't connect to database: $DBI::errstr\n";

#globals
my %smpl;
my %family;

####################################################################################################
#available family relations
#my %relations = (proband => 1, mother => 1, father => 1, parent => 1, brother => 1, sister => 1, sibling => 1, grandfather => 1, grandmother => 1, son => 1, daughter => 1, child => 1, other => 1, aunt => 1, uncle => 1, niece => 1, nephew => 1, cousin => 1, other => 1, spouse => 1);
my %relations;
{
	my $sth = $dbh->prepare("show columns from human_related_samples where Field='relation'") || die "ERROR: Couldn't prepare statement: " . $dbh->errstr;
	$sth->execute() || die "ERROR: Couldn't execute statement: " . $sth->errstr;
	my $type = $sth->fetchrow_hashref->{Type};
	$sth->finish;
	$type =~ s/^enum\(//;
	$type =~ s/\)$//;
	$type =~ s/\'//g;
	$type =~ s/\s//g;
	my @types = split ',', $type;
	foreach(@types)
	{
		$relations{$_} = 1;
	}
}

####################################################################################################
#starting numbers for related and single for the project specified
my $ssid1;      #starting index for single
my $srid1;      #starting index for related
{
	my $sth = $dbh->prepare("select s.external_source_name name, s.source_type type from sources s inner join projects p on p.id=s.project_id where p.project_name = ?") || die "Couldn't prepare statement: " . $dbh->errstr;
	$sth->execute($prid) || die "ERROR: Couldn't execute statement: " . $sth->errstr;
	my(@sng, @rel);
	while(my $h = $sth->fetchrow_hashref)
	{
		die "ERROR: Couldn't decode serial number suffix from sample name '$h->{name}'\n" if($h->{name} !~ /(\d+)$/);
		push @rel, $1 if($h->{type} eq 'human_related_gatk');
		push @sng, $1 if($h->{type} eq 'human_single_gatk');
	}
	$sth->finish;
  @rel = sort{$b <=> $a} @rel;
  @sng = sort{$b <=> $a} @sng;
	$srid1 = scalar @rel? $rel[0]: 0;
	$ssid1 = scalar @sng? $sng[0]: 0;
	$srid1++;
	$ssid1++;
}
print STDERR "starting numbers for project $prid are: related: $srid1, single: $ssid1\n";

####################################################################################################
#won't need mysql db anymore
$dbh->disconnect; 

####################################################################################################
#reading sampleId <=> filename tsv file
print STDERR "reading tsv file with sampleId <=> filename\n";
open F, "$inpName" || die "ERROR: couldn't open file $inpName for reading\n";
while(<F>)
{
	chomp;
	next if(/^#/);
	
	my @a = split "\t";
	next if(!scalar @a);
	
	#splitting name of the fastq file into parts and updating content of the @a array
	die "ERROR: couldn't decode file name '$a[1]' for sampleId '$a[0]'\n" if($a[1] !~ /(.+)(_[Rr]?[12])(\.f(ast)?q(\.gz)?)$/ && $a[1] !~ /(.+)(_[Rr]?[12])(_\d+\.f(ast)?q(\.gz)?)$/); #the second regex is for files names like <name>_R1_001.fastq.gz
	#print STDERR "INFO: sample's fname split: $1\t$2\t$3\n";
	@a = ($a[0], $1, $2, $3);
		
	die "ERROR: sampleId '$a[0]' already seen\n" if(defined $smpl{$a[0]});
	#$smpl{$a[0]} = $a[1];
	push @{$smpl{shift @a}}, @a;
}
close F;
print STDERR "loaded ".(scalar keys %smpl)." sample ids\n";

####################################################################################################
#checking if fastq files are present
if($fastqcheck)
{
	print STDERR "checking if fastq files are present\n";
	foreach(sort keys %smpl)
	{
		my $r1 = $smpl{$_}[1];
		my $r2 = $r1; $r2 =~ s/1$/2/;
		my $f1 = $smpl{$_}[0].$r1.$smpl{$_}[2];
		my $f2 = $smpl{$_}[0].$r2.$smpl{$_}[2];
		if(! -e "$dir/$f1" || ! -e "$dir/$f2")
		{
			print STDERR "  $_\tWARNING: NO fatsq '$dir/$f1' '$dir/$f2'!; sample removerd.\n";
			delete $smpl{$_};
		}
		else
		{
			#print STDERR "  $_\t$f1, $f2\tok\n";
		}
	}
	print STDERR "left with ".(scalar keys %smpl)." samples having fastq data\n";
}
else
{
	print STDERR "skipping ckecking for presence of fastq files\n";
}

####################################################################################################
#loading family trees from databases.apf.edu.au or provided file
print STDERR "loading family trees from databases.apf.edu.au\n";
foreach(sort keys %smpl)
{
	print STDERR "  $_...\n" if($verbose);
	my @fam;
	if(!defined $familyFile)
	{
		@fam = `wget https://databases.apf.edu.au/ApfBioinformaticsRequest/webservice/requestFamilyTree?externalId='$_'\\&key=UJl6BOLPO0f72cOGGK63M14H4fD3yXYz -qO -`;
	}
	else
	{
		print STDERR "  using $familyFile\n" if($verbose);
		my $fid = `grep -F "|$_|" $familyFile |cut -d "|" -f 9`;
		#print STDERR "$_ => fid=:$fid:\n" if($verbose);
		if(defined $fid && $fid ne '')
		{
			@fam = `head -n 1 $familyFile`;
			push @fam, `awk -F "|" '\$9 == $fid' $familyFile`;
			#print STDERR "@fam\n"; exit;
		}
	}
	shift @fam if(@fam); #get rid of the header
	if(!@fam)
	{
		#print STDERR "  ERROR: no record for sample '$_'. '$_' will be treated as single. It will not be possible to import this sample into https://databases.apf.edu.au/.\n";
		print STDERR "  ERROR: no record for sample '$_'. '$_' will be excluded from the analysis. Add '$_' into https://databases.apf.edu.au/individuals/ if you want it included.\n";
		#exit; #needs to be treated accordingly
	}
	else
	{
		#shift @fam; #get rid of the header
	}
	print STDERR "    has family with ".(scalar @fam)." members\n" if($verbose);
	my $isFirst    = 1;
	foreach(@fam)
	{
		chomp;
		s/\n//g;s/\r//g;
		my @a = split '\|';
		$a[8] = "smpl_$a[1]" if($a[8] eq '');
		if($isFirst)
		{
			if(defined $family{$a[8]})
			{
				print STDERR "    this family was already recorded\n" if($verbose);
				last;				
			}
		}
		$isFirst = 0;

		my @b = split ',', $a[2];
		undef $a[0];
		foreach(@b)
		{
			if(defined $smpl{$_})
			{
				$a[0] = $_;
				last;
			}
		}
		if(defined $a[0])
		{
			$a[4] = lc($a[4]);
			$a[5] = lc($a[5]);
			$a[9] = lc($a[9]);
			$a[9] =~ s/\d+$// if(defined $a[9]); #clean cases with sister1, nice2, etc.
			#if(defined $a[4] && $a[5] && $a[9] && defined $relations{$a[9]})
			if(defined $a[4] && $a[5] && $a[9])
			{
				print STDERR "  $a[0]:$a[9] $a[5]\n" if($verbose);
				print STDERR "  WARNING: unknow gender '$a[4]'\n" if($a[4] ne 'male' && $a[4] ne 'female');
				print STDERR "  WARNING: unknow affected yes/no '$a[5]'\n" if($a[5] ne 'yes' && $a[5] ne 'no');
				print STDERR "  WARNING: unknow family relationship '$a[9]'\n" if(!defined $relations{$a[9]});
				$a[5] = $a[5] eq 'yes'? 'affected': 'unaffected';
				$family{$a[8]}{$a[0]}{sex}      = lc($a[4]);
				$family{$a[8]}{$a[0]}{affected} = lc($a[5]);
				$family{$a[8]}{$a[0]}{relation} = lc($a[9]);
			}
			else
			{
				print STDERR "  ERROR: sample '$a[2]' doesn't have all the necessary data. It will be excluded from the analysis. Update '$a[2]' record at https://databases.apf.edu.au/individuals/ if you want it included.\n";
			}			
		}
		else
		{
			print STDERR "  WARNING: sample '$a[2]':$a[9] is not present in the input file '$inpName'. It will be excluded from the analysis. Add it to '$inpName' if you want it included.\n";
		}
	}#foreach(@fam)
}

####################################################################################################
#checking if all families have a proband
foreach my $fid(keys %family)
{
	my $hasProband = 0;
	foreach my $sid(keys %{$family{$fid}})
	{
		$hasProband = 1 if($family{$fid}{$sid}{relation} eq 'proband');
	}
	if(!$hasProband)
	{
		print STDERR "ERROR: family $fid (samples ".join(',', keys %{$family{$fid}}).") has no proband. Affected samples will not be used in the analysis. Correct this family record at https://databases.apf.edu.au/individuals/.\n";
		delete $family{$fid};
	}
}

####################################################################################################
#preparing dirs - checking relations within affected/unaffected are uniq
foreach my $fid(keys %family)
{
	my %a;
	my %u;
	foreach my $sid(keys %{$family{$fid}})
	{
		$family{$fid}{$sid}{dir} = $family{$fid}{$sid}{relation};
		push @{$a{$family{$fid}{$sid}{relation}}}, $sid if($family{$fid}{$sid}{affected} eq 'affected');
		push @{$u{$family{$fid}{$sid}{relation}}}, $sid if($family{$fid}{$sid}{affected} eq 'unaffected');
	}
	foreach(keys %a)
	{
		if(scalar @{$a{$_}} > 1)
		{
			for(my $i = 0; $i < scalar @{$a{$_}}; $i++)
			{
				$family{$fid}{$a{$_}[$i]}{dir} .= $i;
			}
		}
	}
	foreach(keys %u)
	{
		if(scalar @{$u{$_}} > 1)
		{
			print STDERR "INFO: unafected $_ present ".(scalar @{$u{$_}})." times\n";
			for(my $i = 0; $i < scalar @{$u{$_}}; $i++)
			{
				print STDERR "  INFO: renaming unafected $_ => $_".($i + 1)."\n";
				$family{$fid}{$u{$_}[$i]}{dir} .= $i + 1;
			}
		}
	}
}


####################################################################################################
#summary of the cohorts
print STDERR "summary of the cohorts - families and singles (samples excluded from the analysis are not shown below)\n";
foreach my $fid(keys %family)
{
	print STDERR "family $fid\n";
	foreach my $sid(keys %{$family{$fid}})
	{
		print STDERR "  $sid $family{$fid}{$sid}{relation} $family{$fid}{$sid}{affected} $family{$fid}{$sid}{sex}\n";
	}
}

####################################################################################################
#making directory structure and printing lines for sample_info.csv
#needs to detect single-member families and treat them as singles

print STDERR "generating output files\n";

my $cwd = getcwd();

open QSUB, ">$outName.qsub" || die "ERROR: couldn't open file '$outName.qsub' for writing\n";
open CSV,  ">$outName.csv"  || die "ERROR: couldn't open file '$outName.csv' for writing\n";

print QSUB qq ~#!/bin/bash
#PBS -P u86
#PBS -q copyq
#PBS -l walltime=10:00:00,mem=1GB,ncpus=1
#PBS -W umask=0007
#PBS -W group_list=u86
#PBS -o $outName.qsub.out
#PBS -e $outName.qsub.err

~;
print QSUB "#relative source path used, changing directory to $cwd\ncd $cwd\n\n" if($dir !~ /^\//); #there should not be relative paths, but just in case...

my $username = `getent passwd \$(whoami) | cut -d : -f 5`;
chomp $username;
print CSV "#runs entered by $username on ".(strftime "%d-%m-%Y", localtime)."\n"; #a standard comment
foreach my $fid(keys %family)
{
	if(scalar keys %{$family{$fid}} == 1)
	{
		##################################################################################################
		#single
		my $sid = (keys %{$family{$fid}})[0];
		print STDERR "single $prid\_single$ssid1 - $sid (family $fid)\n";
		#print STDERR "creating directories and copying files\n";
		my $d = "$reads/$prid\_single$ssid1";
		if(-d $d)
		{
			print STDERR "WARNING: directory $d already exists; ";
			if($overwrite)
			{
				print STDERR "deleting\n";
				print QSUB "echo \"delete $d\"\nrm -rf $d\n\n";
			}
			else
			{
				print STDERR "skipping\n";
			}
		}
		print QSUB "echo \"create $d\"\nmkdir -p $d\n\n";

		my $r1 = $smpl{$sid}[1];
		my $r2 = $r1; $r2 =~ s/1$/2/;
		my $f1 = $smpl{$sid}[0].$r1.$smpl{$sid}[2];
		my $f2 = $smpl{$sid}[0].$r2.$smpl{$sid}[2];

		print QSUB "echo \"$f1, $f2 => $d\"\ncp -n $dir/$f1 $dir/$f2 $d\n\n";
		
		#print STDERR "printing sample_info.csv record\n";
		print CSV "human_single_gatk,$prid,$prid\_single$ssid1,$r1,$seqc,genome,$sid,$family{$fid}{$sid}{affected}\n";
		$ssid1++;
	}
	else
	{
		##################################################################################################
		#related
		#print STDERR "family $fid will be processed as 'related' $prid\_cohort$srid1\n";
		print STDERR "related $prid\_cohort$srid1 (family $fid)\n";
		my @sida;   #affected sample ids
		my @sidu;   #unaffected sample ids
		my @snamea; #affected sample names
		my @snameu; #unaffected sample names
		my @rela;   #affected family relation
		my @relu;   #unaffected family relation
		my $sex;    #proband's sex
		my $affct;  #in the family samples are affeced, unaffeced, or both
		foreach my $sid(sort{$family{$fid}{$a}{relation} eq 'proband'? -1: $family{$fid}{$a}{relation} eq $family{$fid}{$b}{relation}? 0: 1} keys %{$family{$fid}})
		{
			print STDERR "  $sid $family{$fid}{$sid}{sex} $family{$fid}{$sid}{affected} $family{$fid}{$sid}{relation}\n";
			$sex = $family{$fid}{$sid}{sex} if($family{$fid}{$sid}{relation} eq 'proband');
			
			if($family{$fid}{$sid}{affected} eq 'affected')
			{
				$affct |= 0b10; #set second bit to mark affected
				push @sida, $sid;
				push @snamea, "cohort$srid1\_affected_$family{$fid}{$sid}{dir}";
				push @rela, $family{$fid}{$sid}{relation};
			}
			else
			{
				$affct |= 0b01; #set first bit to mark unaffected
				push @sidu, $sid;
				push @snameu, "cohort$srid1\_unaffected_$family{$fid}{$sid}{dir}";
				push @relu, $family{$fid}{$sid}{relation};
			}
		}

		##################################################################################################
		#print STDERR "creating directories and copying files\n";
		my $d = "$reads/$prid\_cohort$srid1";
		if(-d $d)
		{
			print STDERR "WARNING: directory $d already exists; ";
			if($overwrite)
			{
				print STDERR "deleting\n";
				print QSUB "echo \"delete $d\"\nrm -rf $d\n\n";
			}
			else
			{
				print STDERR "skipping\n";
			}
		}
		for(my $i = 0; $i < scalar @sida; $i++)
		{
			#print STDERR "  $snamea[$i]...\n";
			print QSUB "echo \"create $reads/$prid\_cohort$srid1/$snamea[$i]\"\nmkdir -p $reads/$prid\_cohort$srid1/$snamea[$i]\n\n";
			
			my $r1 = $smpl{$sida[$i]}[1];
			my $r2 = $r1; $r2 =~ s/1$/2/;
			my $f1 = $smpl{$sida[$i]}[0].$r1.$smpl{$sida[$i]}[2];
			my $f2 = $smpl{$sida[$i]}[0].$r2.$smpl{$sida[$i]}[2];

			print QSUB "echo \"$f1, $f2 => $reads/$prid\_cohort$srid1/$snamea[$i]\"\ncp -n $dir/$f1 $dir/$f2 $reads/$prid\_cohort$srid1/$snamea[$i]\n\n";
		}
		for(my $i = 0; $i < scalar @sidu; $i++)
		{
			#print STDERR "  $snameu[$i]...\n";
			print QSUB "echo \"create $reads/$prid\_cohort$srid1/$snameu[$i]\"\nmkdir -p $reads/$prid\_cohort$srid1/$snameu[$i]\n\n";

			my $r1 = $smpl{$sidu[$i]}[1];
			my $r2 = $r1; $r2 =~ s/1$/2/;
			my $f1 = $smpl{$sidu[$i]}[0].$r1.$smpl{$sidu[$i]}[2];
			my $f2 = $smpl{$sidu[$i]}[0].$r2.$smpl{$sidu[$i]}[2];

			print QSUB "echo \"$f1, $f2 => $reads/$prid\_cohort$srid1/$snameu[$i]\"\ncp -n $dir/$f1 $dir/$f2 $reads/$prid\_cohort$srid1/$snameu[$i]\n\n";
		}
		
		##################################################################################################
		#print STDERR "printing sample_info.csv record\n";
		my $r1 = $smpl{$sida[0]}[1]; #a nasty hack - assumes all fastqs in a cohort follow the same pattern
		print CSV "human_related_gatk,$prid,$prid\_cohort$srid1,$r1,$seqc,genome,";
		print CSV join(",", join(":", @sida), join(":", @sidu), join(":", @snamea), , join(":", @snameu), join(":", @rela), join(":", @relu));
		print CSV ",$sex,".($affct == 0b11? 'both': $affct == 0b10? 'affected': 'unaffected')."\n";

		$srid1++;
	}
}
print QSUB "echo \"DONE.\"\n";
close QSUB;
close CSV;
print STDERR "all done\n";
