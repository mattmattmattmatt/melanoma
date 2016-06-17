package modules::Adaptors::Sample;

use strict;
use base 'modules::Adaptors::MelDB';
use modules::ConfigXML;

my $proj_conf = "/drive2/melanoma/v1.1/conf/pipe.xml"; #Hack for epic
#my $proj_conf = $ENV{'SVNDIR'}.'/conf/pipe.xml' ;

my $config = modules::ConfigXML->new($proj_conf);
my $db_name = $config->read('database','db_name');

modules::Adaptors::Sample->table('samples');
modules::Adaptors::Sample->columns(All => qw/id sample_number tissue tumour cell_line tumour_type total_lanes sample_type sample_name bpa_id sample_group_id/);
modules::Adaptors::Sample->has_a(sample_group_id => 'modules::Adaptors::Sample_Group');
modules::Adaptors::Sample->set_sql('all' => 'Select * from samples');
modules::Adaptors::Sample->set_sql('sample_type' => "SELECT column_type FROM information_schema.columns WHERE table_name = \'samples\' AND column_name = \'sample_type\' AND table_schema=\'$db_name\'");

1;