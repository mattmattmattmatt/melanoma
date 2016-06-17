package modules::Adaptors::Lane;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Lane->table('lanes');
modules::Adaptors::Lane->columns(All => qw/id lane_number lane_name sequencing_centre_id sample_id/);
modules::Adaptors::Lane->has_a(sample_id => 'modules::Adaptors::Sample');
modules::Adaptors::Lane->has_a(sequencing_centre_id => 'modules::Adaptors::Sequencing_Centre');

1;