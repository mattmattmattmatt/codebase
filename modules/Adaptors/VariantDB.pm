package modules::Adaptors::VariantDB;

use strict;
use warnings;
use modules::Adaptors::Pass;
use modules::ConfigXML;
use modules::Exception;
use base 'Class::DBI';
use FindBin;

my $proj_conf;
if (exists $ENV{'SVNDIR'}) {
	$proj_conf = $ENV{'SVNDIR'}.'/conf/pipe.xml' ;
} else {
	$proj_conf = '/drive2/variantdb/trunk/conf/pipe.xml';
  #$proj_conf = "/Users/mattfield/work/variantdb_trunk/conf/pipe.xml";
}

#my $proj_conf = $ENV{'SVNDIR'}.'/conf/pipe.xml' ;
#my $proj_conf = "/drive2/variantdb/trunk/conf/pipe.xml"; #Hack for epic
#my $proj_conf = "/Users/mattfield/work/variantdb_trunk/conf/pipe.xml"; #Hack for epic

my ($dsn, $username, $password) = getConfig();

modules::Adaptors::VariantDB->set_db('Main',
				      $dsn,
				      $username,
				      $password,
				      {AutoCommit=>1},
    );

sub getConfig {
	#Get the cluster specific db info
	my $config = modules::ConfigXML->new($proj_conf);
	my $db = $config->read('common','database','db_name');
	my $host = $config->read('common','database','host');
	my $user = $config->read('common','database','user');

	#Need to change this if db different
	my $pw = modules::Adaptors::Pass->decode();
    return ("dbi:mysql:$db:$host",$user, $pw);
}

return 1;
