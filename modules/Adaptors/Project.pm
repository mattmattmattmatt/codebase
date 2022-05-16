package modules::Adaptors::Project;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Project->table('projects');
modules::Adaptors::Project->columns(All => qw/id project_name project_description/);
modules::Adaptors::Project->set_sql('all' => 'Select * from projects');

1;