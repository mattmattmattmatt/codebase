#! /usr/bin/perl -w

use strict;
use modules::Exception;
use modules::Vcf;
use modules::SystemCall;
use Getopt::Long;
use File::Basename;
use Pod::Usage;
use Data::Dumper;
use Cwd 'abs_path';
use vars qw(%OPT);

# Command line arguments
GetOptions(\%OPT, 
	   "help|h",
	   "man|m",
	   "vcf_in=s",
	   "ref=s",
	   "outdir=s",
	   "outbase=s",
	   "gene_anno_file=s",
	   "gene_anno_coord=s",
	   "no_run",
	   "joint"
	   );

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help} || !$OPT{vcf_in});

=pod

=head1 SYNOPSIS

quick_manta.pl -vcf_in output_dir -outdir outdir(default=cwd) -joint joint_calls(default=somatic) -ref ref_genome(default=GRCh38) -no_run list_commands_and_quit -gene_anno_file gene_annotation_file -outbase outbase_filename(default=vcf_manta) 

Required flags: -vcf_in 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 USAGE

quick_manta.pl -> annotate a vcf from manta

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

=cut

my $svndir;
if (!exists $ENV{'SVNDIR'}) {
	modules::Exception->throw("ERROR:  You need to set system paths and svn directory by running 'source ../conf/export_env.txt ..' from the current working svn/scripts directory\n");
} else {
	$svndir = $ENV{'SVNDIR'};
}

my $joint_flag = defined $OPT{joint}?'-joint':''; #default is somatic



if (!-d $svndir) {
	modules::Exception->throw("ERROR: $svndir doesn't exist\n");
}

my $ref = defined $OPT{ref}?$OPT{ref}:"GRCh38";

my $gene_anno_file = my $gene_coord_file;

my $conf_dir = "$svndir/conf/human/".$ref;

my $gene_dir = &GetLatest('gene');

if (defined $OPT{gene_anno_file}) {
	$gene_anno_file = $OPT{gene_anno_file};
} else {
	my $file = $ref . '.gene.all';
	$gene_anno_file = $gene_dir . '/'.$file;
} 

if ($OPT{gene_coord_file}) {
	$gene_coord_file = $OPT{gene_coord_file};
} else {
	my $file = $ref . '.gene.overlap.all';
	$gene_coord_file = $gene_dir . '/'.$file;
}
	
if ( !-e $gene_coord_file ) {
	modules::Exception->throw("File $gene_coord_file doesn't exist");	
}
if ( !-e $gene_anno_file ) {
	modules::Exception->throw("File $gene_anno_file doesn't exist");	
}

my $vcf = $OPT{vcf_in};
$vcf = abs_path($vcf);

if ( !-e $vcf ) {
	modules::Exception->throw("File $vcf doesn't exist");	
}

my $vcf_base;

if ($OPT{outbase}) {
	$vcf_base = $OPT{outbase};
} else {
	my ($vcf_short) = basename($vcf);
	($vcf_base = $vcf_short) =~ s/.vcf/_manta/;
}



my $outdir = defined $OPT{outdir}?$OPT{outdir}:`pwd`;
chomp $outdir;
$outdir = abs_path($outdir);

if (!-d $outdir) {
        modules::Exception->throw("Dir $outdir doesn't exist");
}


my $parse_manta = "$svndir/utils/parse_manta.pl";
my $summarise_manta = "$svndir/utils/manta_summarise.pl";
my $overlap = "$svndir/utils/overlap_files.pl";


#Build up command list
my @commands = ();

push @commands, "$parse_manta -manta_vcf $vcf -outdir $outdir -outfile_base $vcf_base $joint_flag";
push @commands, "grep -v inv $outdir/${vcf_base}.txt > $outdir/${vcf_base}.nobp";
push @commands, "cp $outdir/${vcf_base}.tra $outdir/${vcf_base}.bp";
push @commands, "grep inv $outdir/${vcf_base}.txt >> $outdir/${vcf_base}.bp";
push @commands, "grep bnd $outdir/${vcf_base}.bp".' | awk \''.'{print $1,$2,$2,$5"\n"$3,$4,$4,$5}\'' ." > $outdir/${vcf_base}.bp.split";
push @commands, "grep inv $outdir/${vcf_base}.bp".' | awk \''.'{print $1,$2,$2,$4"\n"$1,$3,$3,$4}\'' ." >> $outdir/${vcf_base}.bp.split";


push @commands, "$overlap -ref $outdir/${vcf_base}.nobp -just_overlap -coord  $gene_coord_file -all | sed -e \"s/\\t1\$/\\tNO_GENE/\" > $outdir/${vcf_base}.nobp.gene";
push @commands, "$overlap -ref $outdir/${vcf_base}.bp.split -just_overlap -coord  $gene_coord_file -all | sed -e \"s/\\t1\$/\\tNO_GENE/\" > $outdir/${vcf_base}.bp.gene";
push @commands, "rm -f $outdir/${vcf_base}.all";
push @commands, "rm -f $outdir/${vcf_base}.all.gene";
push @commands, "cat $outdir/${vcf_base}.txt  $outdir/${vcf_base}.tra >> $outdir/${vcf_base}.all";
push @commands, "cat $outdir/${vcf_base}.nobp.gene $outdir/${vcf_base}.bp.gene >> $outdir/${vcf_base}.all.gene";
push @commands, "$summarise_manta -outdir $outdir -gene_file $gene_anno_file -outfile_base $vcf_base $joint_flag > $outdir/${vcf_base}_annotated.tsv";


my $sys_call = modules::SystemCall->new();

#Run the commands
for my $command (@commands) {
	print "$command\n";
	`$command` unless $OPT{no_run};
}


#Get the latest directory based on formats DDMMYY
sub GetLatest {
	my ($name) = @_;
	opendir(DIR,"$conf_dir/$name/") || modules::Exception->throw("ERROR: Cannot open directory $conf_dir/$name/");
	my @files = grep {/^\d/} readdir DIR;
	closedir DIR;
	my ($dir_tmp) = reverse(sort {my ($aday,$amonth,$ayear) = $a =~ /(\d\d)(\d\d)(\d\d)/; my ($bday,$bmonth,$byear) = $b =~ /(\d\d)(\d\d)(\d\d)/; $ayear<=>$byear||$amonth<=>$bmonth||$aday<=>$bday} @files);
	return "$conf_dir/$name/$dir_tmp";
}



