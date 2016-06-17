package modules::Adaptors::Sample_Group;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Sample_Group->table('sample_groups');
modules::Adaptors::Sample_Group->columns(All => qw/id total_samples sample_group_number sample_group_name patient_id/);
modules::Adaptors::Sample_Group->has_a(patient_id => 'modules::Adaptors::Patient');
modules::Adaptors::Sample_Group->set_sql('all' => 'Select * from sample_groups');


1;