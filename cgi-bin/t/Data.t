#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Data');
}
require_ok('miRkwood::Data');

## get_mirbase_file ##
ok( my $mirbase_file = miRkwood::Data::get_mirbase_file(),
    'Can call get_mirbase_file()'
);

## get_matrix_file ##
ok( my $matrix_file = miRkwood::Data::get_matrix_file(),
    'Can call get_matrix_file()'
);

## get_mirdup_data_path ##
ok( my $mirdup_data_path = miRkwood::Data::get_mirdup_data_path(),
    'Can call get_mirdup_data_path()'
);

## get_mirdup_model_name ##
ok( my $model_name = miRkwood::Data::get_mirdup_model_name(),
    'Can call get_mirdup_model_name()'
);
my $model_name_expected = 'plant.model';
is( $model_name, $model_name_expected,
    'get_mirdup_model_name() returns the expected value');
