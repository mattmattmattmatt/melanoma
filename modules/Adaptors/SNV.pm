package modules::Adaptors::SNV;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::SNV->table('snvs');
modules::Adaptors::SNV->columns('All' => qw/id chr coord ref_base var_base ref_base_freq var_base_freq median_quality_score snv_score clr_score read_depth normal_base_string tumour_base_string snv_class run_id/);
modules::Adaptors::SNV->has_a(run_id => 'modules::Adaptors::Run');
modules::Adaptors::SNV->set_sql('by_experiment_and_passed_filter' => 'SELECT snvs.id FROM snvs, snvs_filters WHERE snvs.run_id = ? AND snvs.id = snvs_filters.snv_id AND snvs_filters.filter_id = ? AND snvs_filters.filterpass = 1 ORDER BY chr, coord');
modules::Adaptors::SNV->set_sql('chr_by_experiment_and_passed_filter' => 'SELECT snvs.id FROM snvs, snvs_filters WHERE snvs.run_id = ? AND snvs.id = snvs_filters.snv_id AND snvs_filters.filter_id = ? AND snvs_filters.filterpass = 1 AND chr = ? ORDER BY chr, coord');
modules::Adaptors::SNV->set_sql('region_by_experiment_and_passed_filter' => 'SELECT snvs.id FROM snvs, snvs_filters WHERE snvs.run_id = ? AND snvs.id = snvs_filters.snv_id AND snvs_filters.filter_id = ? AND snvs_filters.filterpass = 1 AND chr = ? AND coord >= ? AND coord <= ? ORDER BY chr, coord');
modules::Adaptors::SNV->set_sql('snv_ids' => 'SELECT snvs.id FROM snvs WHERE run_id = ?');


return 1;
