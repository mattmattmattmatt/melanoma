package modules::Adaptors::Patient;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Patient->table('patients');
modules::Adaptors::Patient->columns(All => qw/id patient_external_name/);
modules::Adaptors::Patient->set_sql('all' => 'Select * from patients');

1;