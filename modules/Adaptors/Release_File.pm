package modules::Adaptors::Release_File;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Release_File->table('release_files');
modules::Adaptors::Release_File->columns(All => qw/id file_name file_type total_lanes pipeline_steps_run_id/);
modules::Adaptors::Release_File->has_a(pipeline_steps_run_id => 'modules::Adaptors::Pipeline_Step_Run');
modules::Adaptors::Release_File->set_sql('release_files_run'=>'Select file_name from release_files,pipeline_steps_runs where pipeline_steps_runs.id=release_files.pipeline_steps_run_id AND pipeline_steps_runs.run_id=?');
modules::Adaptors::Release_File->set_sql('release_files_sample'=>'Select file_name from release_files,pipeline_steps_runs,runs where pipeline_steps_runs.id=release_files.pipeline_steps_run_id AND pipeline_steps_runs.run_id=runs.id AND runs.sample_id=?');
modules::Adaptors::Release_File->set_sql(total_lanes => 'Select total_lanes from release_files,pipeline_steps_runs where pipeline_steps_runs.run_id=? AND pipeline_steps_runs.pipeline_step_id=? AND pipeline_steps_runs.id = release_files.pipeline_steps_run_id');


1;