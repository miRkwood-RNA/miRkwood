#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use File::Temp;
use File::Spec;
use File::Basename;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('miRkwood::FileUtils');
}
require_ok('miRkwood::FileUtils');

## Setup ##

my $tmp_dir = File::Temp::tempdir(CLEANUP => 1);
my $sub_dir1 = File::Spec->catdir($tmp_dir, '1');
my $sub_dir2 = File::Spec->catdir($tmp_dir, '2');
my $sub_dir3 = File::Spec->catdir($tmp_dir, '10');
mkdir $sub_dir1;
mkdir $sub_dir2;
mkdir $sub_dir3;
my ($FH, $tmp_file) = File::Temp::tempfile( DIR => $tmp_dir );
print {$FH} 'test';
close $FH or die('Error closing');
my $tmp_file_name = File::Basename::basename($tmp_file);


## is_directory() ##

my $res1 = miRkwood::FileUtils::is_directory('1', $tmp_dir);
ok ($res1, 'is_directory returns True on a directory');
my $res2 = miRkwood::FileUtils::is_directory('3', $tmp_dir);
ok (!$res2, 'is_directory returns False on an unexisting item');
my $res3 = miRkwood::FileUtils::is_directory($tmp_file_name, $tmp_dir);
ok (!$res3, 'is_directory returns False on a file');


## get_dirs_from_directory() ##

ok(
    my @contents =
      miRkwood::FileUtils::get_dirs_from_directory($tmp_dir),
    'Can call get_dirs_from_directory()'
);
my @expected = qw(  1 2 10);
is_deeply(\@contents, \@expected,
          'get_dirs_from_directory returns the correct values');
