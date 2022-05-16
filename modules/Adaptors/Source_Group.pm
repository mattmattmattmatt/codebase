package modules::Adaptors::Source_Group;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Source_Group->table('source_groups');
modules::Adaptors::Source_Group->columns(All => qw/id total_samples source_group_number source_group_name source_id/);
modules::Adaptors::Source_Group->has_a(source_id => 'modules::Adaptors::Source');
modules::Adaptors::Source_Group->set_sql('all' => 'Select * from source_groups');


1;