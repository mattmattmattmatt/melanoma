package modules::Adaptors::SNV_Row;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::SNV_Row->table('snv_rows');
modules::Adaptors::SNV_Row->columns('All' => qw/id snv_id run_id/);
modules::Adaptors::SNV_Row->has_a(snv_id => 'modules::Adaptors::SNV');
modules::Adaptors::SNV_Row->has_a(run_id => 'modules::Adaptors::Run');

return 1;
