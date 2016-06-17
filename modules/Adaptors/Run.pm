package modules::Adaptors::Run;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Run->table('runs');
modules::Adaptors::Run->columns(All => qw/id run_directory status commenced production sample_id/);
modules::Adaptors::Run->has_a(sample_id => 'modules::Adaptors::Sample');
modules::Adaptors::Run->set_sql('latest_run_date' => ' Select * from runs where sample_id= ? order by id DESC limit 1 ');

1;

