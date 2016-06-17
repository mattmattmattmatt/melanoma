package modules::Adaptors::Lane_Run;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Lane_Run->table('lanes_runs');
modules::Adaptors::Lane_Run->columns(All => qw/id lane_id run_id/);
modules::Adaptors::Lane_Run->has_a(lane_id => 'modules::Adaptors::Lane');
modules::Adaptors::Lane_Run->has_a(run_id => 'modules::Adaptors::Run');

1;