#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Parsers');
}
require_ok('PipelineMiRNA::Parsers');

# Vienna line

my $input1    = '(((((.(((.(((((.))....))).)))...))))) (-16.93)';
my @expected1 = [ '(((((.(((.(((((.))....))).)))...)))))', '-16.93' ];
my @result1   = PipelineMiRNA::Parsers::parse_Vienna_line($input1);
is_deeply( \@result1, @expected1, 'Standard line correctly parsed' );

my $input2    = '(((((.(((.(((((.))....))).)))...))))) ( -8.58)';
my @expected2 = [ '(((((.(((.(((((.))....))).)))...)))))', '-8.58' ];
my @result2   = PipelineMiRNA::Parsers::parse_Vienna_line($input2);
is_deeply( \@result2, @expected2,
    'Line with whitespace in energy correctly parsed' );

# P_value

my $p_value_file = input_file('pvalue.txt');
file_exists_ok($p_value_file);
my $expected3 = '0.125000';
my $result3   = PipelineMiRNA::Parsers::parse_pvalue($p_value_file);
is( $result3, $expected3, 'p_value file is correctly parsed' );

# MFEI

my $mfei_file = input_file('outMFEI.txt');
file_exists_ok($mfei_file);
my @expected4 = [ '-1.00166666666667', '-60.1', '-37.0987654320988' ];
my @result4 = PipelineMiRNA::Parsers::parse_mfei($mfei_file);
is_deeply( \@result4, @expected4, 'MFEI file is correctly parsed' );

