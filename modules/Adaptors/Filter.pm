package modules::Adaptors::Filter;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Filter->table('filters');
modules::Adaptors::Filter->columns(All => qw/id name description pipeline_step_id/);
modules::Adaptors::Filter->has_a(pipeline_step_id => 'modules::Adaptors::Pipeline_Step');

1;
