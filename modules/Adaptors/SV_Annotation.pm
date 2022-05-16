package modules::Adaptors::SV_Annotation;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::SV_Annotation->table('structural_variantrow_values');
modules::Adaptors::SV_Annotation->columns(All => qw/id column_number column_name column_value structural_variant_row_id/);
modules::Adaptors::SV_Annotation->has_a(structural_variant_row_id => 'modules::Adaptors::SV_Row');
1;
