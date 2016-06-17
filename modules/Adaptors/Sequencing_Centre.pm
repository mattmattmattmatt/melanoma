package modules::Adaptors::Sequencing_Centre;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Sequencing_Centre->table('sequencing_centres');
modules::Adaptors::Sequencing_Centre->columns(All => qw/id sequencing_centre_name/);
modules::Adaptors::Sequencing_Centre->set_sql(sequencing_centre => 'Select * from sequencing_centres,lanes where lanes.id=? AND lanes.sequencing_centre_id=sequencing_centres.id');


1;