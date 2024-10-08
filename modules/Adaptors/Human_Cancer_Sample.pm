package modules::Adaptors::Human_Cancer_Sample;

use strict;
use base 'modules::Adaptors::VariantDB';
use modules::ConfigXML;

modules::Adaptors::Human_Cancer_Sample->table('human_cancer_samples');
modules::Adaptors::Human_Cancer_Sample->columns(All => qw/id tissue tumour cell_line tumour_type cancer_type sample_id/);
modules::Adaptors::Human_Cancer_Sample->has_a(sample_id => 'modules::Adaptors::Sample');
modules::Adaptors::Human_Cancer_Sample->set_sql('all' => 'Select * from human_cancer_samples');

1;