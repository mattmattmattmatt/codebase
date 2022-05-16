package modules::Adaptors::SV;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::SV->table('structural_variants');
modules::Adaptors::SV->columns('All' => qw/id chr1 chr2 start_coord1 end_coord1 start_coord2 end_coord2 sv_type sv_score supporting_read_pairs sv_class run_id/);
modules::Adaptors::SV->has_a(run_id => 'modules::Adaptors::Run');
modules::Adaptors::SV->set_sql('by_experiment_and_passed_filter' => 'SELECT structural_variants.id FROM structural_variants, structural_variants_filters WHERE structural_variants.run_id = ? AND structural_variants.id = structural_variants_filters.snv_id AND structural_variants_filters.filter_id = ? AND structural_variants_filters.filterpass = 1 ORDER BY chr1, start_coord1');
modules::Adaptors::SV->set_sql('chr_by_experiment_and_passed_filter' => 'SELECT structural_variants.id FROM structural_variants, structural_variants_filters WHERE structural_variants.run_id = ? AND structural_variants.id = structural_variants_filters.snv_id AND structural_variants_filters.filter_id = ? AND structural_variants_filters.filterpass = 1 AND chr = ? ORDER BY chr1, start_coord1');
modules::Adaptors::SV->set_sql('region_by_experiment_and_passed_filter' => 'SELECT structural_variants.id FROM structural_variants, structural_variants_filters WHERE structural_variants.run_id = ? AND structural_variants.id = structural_variants_filters.snv_id AND structural_variants_filters.filter_id = ? AND structural_variants_filters.filterpass = 1 AND chr1 = ? AND start_coord1 >= ? AND start_coord1 <= ? ORDER BY chr1, start_coord1');
modules::Adaptors::SV->set_sql('sv_ids' => 'SELECT structural_variants.id FROM structural_variants WHERE run_id = ?');


return 1;
