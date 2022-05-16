package modules::Adaptors::SV_Row;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::SV_Row->table('structural_variant_rows');
modules::Adaptors::SV_Row->columns('All' => qw/id structural_variant_id run_id/);
modules::Adaptors::SV_Row->has_a(structural_variant_id => 'modules::Adaptors::SV');
modules::Adaptors::SV_Row->has_a(run_id => 'modules::Adaptors::Run');

return 1;
