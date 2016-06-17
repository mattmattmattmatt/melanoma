package modules::Adaptors::Variant;

use strict;
use base 'modules::Adaptors::MelDB';

modules::Adaptors::Variant->table('variants');
modules::Adaptors::Variant->columns('All' => qw/id chr start_coord end_coord var_type var_score clr_score inserted_bases var_caller read_depth ref_base_freq var_base_freq normal_base_string tumour_base_string median_quality_score variant_class run_id/);
modules::Adaptors::Variant->has_a(run_id => 'modules::Adaptors::Run');
modules::Adaptors::Variant->set_sql('by_experiment_and_passed_filter' => 'SELECT variants.id FROM variants, variants_filters WHERE variants.run_id = ? AND variants.id = variants_filters.variant_id AND variants_filters.filter_id = ? AND variants_filters.filterpass = 1 ORDER BY chr, start_coord');
modules::Adaptors::Variant->set_sql('by_experiment_and_passed_filter_and_type' => 'SELECT variants.id FROM variants, variants_filters WHERE variants.run_id = ? AND variants.id = variants_filters.variant_id AND variants_filters.filter_id = ? AND variants_filters.filterpass = 1 AND variants.var_type = ? ORDER BY chr, start_coord');
modules::Adaptors::Variant->set_sql('chr_by_experiment_and_passed_filter' => 'SELECT variants.id FROM variants, variants_filters WHERE variants.run_id = ? AND variants.id = variants_filters.variant_id AND variants_filters.filter_id = ? AND variants_filters.filterpass = 1 AND chr = ? ORDER BY chr, start_coord');
modules::Adaptors::Variant->set_sql('region_by_experiment_and_passed_filter' => 'SELECT variants.id FROM variants, variants_filters WHERE variants.run_id = ? AND variants.id = variants_filters.variant_id AND variants_filters.filter_id = ? AND variants_filters.filterpass = 1 AND chr = ? AND start_coord >= ? AND start_coord <= ? ORDER BY chr, start_coord');
modules::Adaptors::Variant->set_sql('variant_ids' => 'SELECT variants.id FROM variants WHERE run_id = ?');
modules::Adaptors::Variant->set_sql('variant_ids_type' => 'SELECT variants.id FROM variants WHERE run_id = ? AND variants.var_type = ?');

return 1;
