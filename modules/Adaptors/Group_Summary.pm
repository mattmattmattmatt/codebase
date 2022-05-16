package modules::Adaptors::Group_Summary;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Group_Summary->table('group_summaries');
modules::Adaptors::Group_Summary->columns(All => qw/id file_name total_samples summary_date source_group_id/);
modules::Adaptors::Group_Summary->has_a(source_group_id => 'modules::Adaptors::Source_Group');


1;