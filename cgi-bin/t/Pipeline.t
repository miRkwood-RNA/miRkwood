#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;
use Test::Exception;
use Data::Dumper;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::Pipeline');
}
require_ok('miRkwood::Pipeline');

## is_included ##

my @fail_value = ( 1, 7, 2, 9 );
dies_ok { miRkwood::Pipeline::is_included(@fail_value); }
'is_included dies if positions are not ordered';

dies_ok { miRkwood::Pipeline::is_included( 2, 7, 1 ); }
'is_included dies if not enough values are provided';

my @is_included_values = (
    [ [ 2, 7,  1, 9  ], 1 ],
    [ [ 2, 12, 1, 9  ], '' ],
);

foreach my $couple (@is_included_values) {
    my @input           = @{ @{$couple}[0] };
    my $expected        = @{$couple}[1];
    my $is_to_merge_res = miRkwood::Pipeline::is_included(@input);
    is( $is_to_merge_res, $expected, "is_included (@input) --> $expected ok" );
}

## is_overlapping ##

dies_ok { miRkwood::Pipeline::is_overlapping(@fail_value); }
'is_overlapping dies if positions are not ordered';

dies_ok { miRkwood::Pipeline::is_overlapping( 2, 7, 1 ); }
'is_overlapping dies if not enough values are provided';

my @is_overlapping_values = (
    [ [ 2, 7,  1, 9 ], 1 ],
    [ [ 1, 12, 1, 9 ], 1 ],
    [ [ 4, 12, 1, 9 ], 1 ],
    [ [ 5, 12, 1, 9 ], '' ],
);

foreach my $couple (@is_overlapping_values) {
    my @input           = @{ @{$couple}[0] };
    my $expected        = @{$couple}[1];
    my $is_to_merge_res = miRkwood::Pipeline::is_overlapping(@input);
    is( $is_to_merge_res, $expected,
        "is_overlapping (@input) --> $expected ok" );
}
