package modules::Adaptors::SNV_Annotation;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::SNV_Annotation->table('snvrow_values');
modules::Adaptors::SNV_Annotation->columns(All => qw/id column_number column_name column_value snv_row_id/);
modules::Adaptors::SNV_Annotation->has_a(snv_row_id => 'modules::Adaptors::SNV_Row');
1;
