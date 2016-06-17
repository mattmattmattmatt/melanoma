package modules::Adaptors::SNV_Filter;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::SNV_Filter->table('snvs_filters');
modules::Adaptors::SNV_Filter->columns(All => qw/id filtermatch filterpass attribute filter_id snv_id/);
modules::Adaptors::SNV_Filter->has_a(snv_id => 'modules::Adaptors::SNV');
modules::Adaptors::SNV_Filter->has_a(filter_id => 'modules::Adaptors::Filter');


1;
