package modules::Adaptors::Human_Related_Sample;

use strict;
use base 'modules::Adaptors::VariantDB';
use modules::ConfigXML;

modules::Adaptors::Human_Related_Sample->table('human_related_samples');
modules::Adaptors::Human_Related_Sample->columns(All => qw/id affected relation sex sample_id/);
modules::Adaptors::Human_Related_Sample->has_a(sample_id => 'modules::Adaptors::Sample');
modules::Adaptors::Human_Related_Sample->set_sql('all' => 'Select * from human_related_samples');

1;