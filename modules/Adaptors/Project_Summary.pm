package modules::Adaptors::Project_Summary;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Project_Summary->table('project_summaries');
modules::Adaptors::Project_Summary->columns(All => qw/id file_name total_samples summary_date project_id/);
modules::Adaptors::Project_Summary->has_a(project_id => 'modules::Adaptors::Project');


1;