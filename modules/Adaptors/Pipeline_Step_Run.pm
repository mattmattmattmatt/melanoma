package modules::Adaptors::Pipeline_Step_Run;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Pipeline_Step_Run->table('pipeline_steps_runs');
modules::Adaptors::Pipeline_Step_Run->columns(All => qw/id pipeline_step_id run_id/);
modules::Adaptors::Pipeline_Step_Run->has_a(run_id => 'modules::Adaptors::Run');
modules::Adaptors::Pipeline_Step_Run->has_a(pipeline_step_id => 'modules::Adaptors::Pipeline_Step');

1;