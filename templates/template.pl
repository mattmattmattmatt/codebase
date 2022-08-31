#! /usr/bin/perl -w

use strict;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use modules::Exception;
use modules::Utils;
use vars qw(%OPT);

GetOptions(\%OPT, 
	   	"help|h",
	   	"man|m",
	   	);

pod2usage(-verbose => 2) if $OPT{man};
pod2usage(1) if ($OPT{help});



=pod

=head1 SYNOPSIS

FILE.pl -help -man [options]

Required flags: 

=head1 OPTIONS

    -help  brief help message

    -man   full documentation

=head1 NAME

FILE.pl 

=head1 DESCRIPTION

date

a script that ...

=head1 AUTHOR

Matthew Field

=head1 EXAMPLE

./FILE.pl    

=cut

