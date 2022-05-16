package modules::Adaptors::Lane;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Lane->table('lanes');
modules::Adaptors::Lane->columns(All => qw/id lane_number lane_name sequencing_centre sample_id/);
modules::Adaptors::Lane->has_a(sample_id => 'modules::Adaptors::Sample');

1;