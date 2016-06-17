package modules::Adaptors::Variant_Annotation;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Variant_Annotation->table('variantrow_values');
modules::Adaptors::Variant_Annotation->columns(All => qw/id column_number column_name column_value variant_row_id/);
modules::Adaptors::Variant_Annotation->has_a(variant_row_id => 'modules::Adaptors::Variant_Row');

1;
