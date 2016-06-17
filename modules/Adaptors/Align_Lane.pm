package modules::Adaptors::Align_Lane;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Align_Lane->table('align_lanes');
modules::Adaptors::Align_Lane->columns(All => qw/id quality_encoding lane_id/);
modules::Adaptors::Align_Lane->has_a(lane_id => 'modules::Adaptors::Lane');

1;