package modules::Adaptors::Variant_Filter;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Variant_Filter->table('variants_filters');
modules::Adaptors::Variant_Filter->columns(All => qw/id filtermatch filterpass attribute filter_id variant_id/);
modules::Adaptors::Variant_Filter->has_a(variant_id => 'modules::Adaptors::Variant');
modules::Adaptors::Variant_Filter->has_a(filter_id => 'modules::Adaptors::Filter');


1;
