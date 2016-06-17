package modules::Adaptors::Syscall;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Syscall->table('syscalls');
modules::Adaptors::Syscall->columns(All => qw/id command level analysis_step_name chr align_lane_id pipeline_steps_run_id/);
modules::Adaptors::Syscall->has_a(align_lane_id => 'modules::Adaptors::Align_Lane');
modules::Adaptors::Syscall->has_a(pipeline_steps_run_id => 'modules::Adaptors::Pipeline_Step_Run');

1;