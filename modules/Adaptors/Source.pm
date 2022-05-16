package modules::Adaptors::Source;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Source->table('sources');
modules::Adaptors::Source->columns(All => qw/id source_type external_source_name organism project_id/);
modules::Adaptors::Source->has_a(project_id => 'modules::Adaptors::Project');
modules::Adaptors::Source->set_sql('all' => 'Select * from sources');
modules::Adaptors::Source->set_sql('latest_source_name' => 'SELECT * FROM sources WHERE project_id = ? AND source_type = ? ORDER by id desc limit 1');

1;
