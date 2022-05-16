package modules::Adaptors::Variant_Filter;

use strict;
use base 'modules::Adaptors::VariantDB';

modules::Adaptors::Variant_Filter->table('structural_variants_filters');
modules::Adaptors::Variant_Filter->columns(All => qw/id filtermatch filterpass attribute filter_id variant_id/);
modules::Adaptors::Variant_Filter->has_a(structural_variant_id => 'modules::Adaptors::SV');
modules::Adaptors::Variant_Filter->has_a(filter_id => 'modules::Adaptors::Filter');


1;
