package modules::Adaptors::MelDB;

use strict;
use warnings;
use modules::Adaptors::Pass;
use modules::ConfigXML;
use modules::Exception;
use base 'Class::DBI';
use FindBin;

#my $proj_conf = $ENV{'SVNDIR'}.'/conf/pipe.xml' ;
my $proj_conf = "/drive2/melanoma/v1.1/conf/pipe.xml"; #Hack for epic
my ($dsn, $username, $password) = getConfig();

modules::Adaptors::MelDB->set_db('Main',
				      $dsn,
				      $username,
				      $password,
				      {AutoCommit=>1},
    );

sub getConfig {
	#Get the cluster specific db info
	my $config = modules::ConfigXML->new($proj_conf);
	my $db = $config->read('database','db_name');
	my $host = $config->read('database','host');
	my $user = $config->read('database','user');

	#Need to change this if db different
	my $pw = modules::Adaptors::Pass->decode();
	#return ("dbi:mysql:melanoma_production_1_0:igl1.anu.edu.au",'melanoma_user', $pw);
    return ("dbi:mysql:$db:$host",$user, $pw);
}

return 1;
