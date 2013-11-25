#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::WebTemplate');
}
require_ok('PipelineMiRNA::WebTemplate');

ok( my $result1 = PipelineMiRNA::WebTemplate::get_link_back_to_results(1),
    'can call get_link_back_to_results()');
is( $result1, './resultsWithID.pl?run_id=1',
    'get_link_back_to_results returns expected value');

ok( my $result2 = PipelineMiRNA::WebTemplate::get_css_file(),
    'can call get_css_file()');
is( $result2, '/arn/style/script.css',
    'get_css_file returns expected value');

ok( my $result3 = PipelineMiRNA::WebTemplate::get_js_file(),
    'can call get_js_file()');
is( $result3, '/arn/js/miARN.js',
    'get_js_file returns expected value');

ok( my $result4 = PipelineMiRNA::WebTemplate::get_error_page("Error"),
    'can call get_error_page()');
my $expected4 = <<"HTML";
Content-type: text/html

<html>
    <head>
        <LINK rel="stylesheet" type="text/css" href="$result2" />
        <script src="$result3" type="text/javascript" LANGUAGE="JavaScript"></script>
        <title>MicroRNA identification</title>
    </head>
    <body>
        Error
    </body>
HTML
is( $result4, $expected4,
    'get_error_page returns expected value');

ok( my $result5 = PipelineMiRNA::WebTemplate::make_mirbase_link('ABCDE'),
    'can call make_mirbase_link()');
is( $result5, 'http://mirbase.org/cgi-bin/mirna_entry.pl?acc=ABCDE',
    'make_mirbase_link returns expected value');

ok( $ENV{SERVER_NAME} = 'toto',
    'Can set SERVER_NAME variable');
ok( my $result6 = PipelineMiRNA::WebTemplate::make_url('ABCDE'),
    'can call make_url()');
is( $result6, 'http://toto/cgi-bin/ABCDE',
    'make_url returns expected value');

