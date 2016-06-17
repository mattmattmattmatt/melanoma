package modules::Adaptors::Pipeline_Step;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Pipeline_Step->table('pipeline_steps');
modules::Adaptors::Pipeline_Step->columns(All => qw/id name description/);
modules::Adaptors::Pipeline_Step->set_sql('all' => 'Select * from pipeline_steps order by id');

1;