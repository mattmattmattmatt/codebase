package modules::Adaptors::Human_Single_Sample;

use strict;
use base 'modules::Adaptors::VariantDB';
use modules::ConfigXML;

modules::Adaptors::Human_Single_Sample->table('human_single_samples');
modules::Adaptors::Human_Single_Sample->columns(All => qw/id affected sample_id/);
modules::Adaptors::Human_Single_Sample->has_a(sample_id => 'modules::Adaptors::Sample');
modules::Adaptors::Human_Single_Sample->set_sql('all' => 'Select * from human_single_samples');

1;