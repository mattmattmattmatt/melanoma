package modules::Adaptors::Read_File;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Read_File->table('read_files');
modules::Adaptors::Read_File->columns(All => qw/id file_name is_compressed compression_suffix read_file_number read_directory read_file_name lane_id/);
modules::Adaptors::Read_File->has_a(lane_id => 'modules::Adaptors::Lane');

1;