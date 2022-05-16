package modules::Adaptors::Sample;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Sample->table('samples');
modules::Adaptors::Sample->columns(All => qw/id sample_number total_lanes sample_name sample_description external_sample_name sample_type sequence_type source_group_id apf_request_id/);
modules::Adaptors::Sample->has_a(source_group_id => 'modules::Adaptors::Source_Group');
modules::Adaptors::Sample->set_sql('all' => 'Select * from samples');
1;