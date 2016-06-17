package modules::Adaptors::Variant_Row;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Variant_Row->table('variant_rows');
modules::Adaptors::Variant_Row->columns('All' => qw/id variant_id run_id/);
modules::Adaptors::Variant_Row->has_a(variant_id => 'modules::Adaptors::Variant');
modules::Adaptors::Variant_Row->has_a(run_id => 'modules::Adaptors::Run');


return 1;
