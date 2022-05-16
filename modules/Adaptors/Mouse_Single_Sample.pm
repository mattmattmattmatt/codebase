package modules::Adaptors::Mouse_Single_Sample;

use strict;
use base 'modules::Adaptors::VariantDB';
use modules::ConfigXML;

modules::Adaptors::Mouse_Single_Sample->table('mouse_single_samples');
modules::Adaptors::Mouse_Single_Sample->columns(All => qw/id strain internal sample_id/);
modules::Adaptors::Mouse_Single_Sample->has_a(sample_id => 'modules::Adaptors::Sample');
modules::Adaptors::Mouse_Single_Sample->set_sql('all' => 'Select * from mouse_single_samples');

1;