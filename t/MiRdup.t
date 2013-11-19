
#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::MiRdup');
}
require_ok('PipelineMiRNA::MiRdup');

my $name_base = "dummy_seq";
my $mature_seq = "ATGCTAGCTGATGCTGATCGTAGTCGATGCTAGCTG";
my $structure  = "..((((((....((((....))))....))))))..";
my @alignments;
push @alignments, '10-19';
push @alignments, '20-28';

my $validation_source_expected = <<END;
${name_base}__$alignments[0]\tATGCTGATC\t$mature_seq\t$structure
${name_base}__$alignments[1]\tTAGTCGAT\t$mature_seq\t$structure
END

ok( my $validation_source = PipelineMiRNA::MiRdup->make_validation_source_file($name_base, $mature_seq, $structure, @alignments),
    'Can call make_validation_source_file' );

is( $validation_source . "\n", $validation_source_expected,
    "make_validation_source_file returns a correct file") ;


my $validation_file = input_file('sequencesToValidate.txt.plant.model.miRdup.tab.txt');
file_exists_ok($validation_file);
ok( my %result = PipelineMiRNA::MiRdup->parse_validation_output($validation_file),
    'Can call parse_validation_output');
my %expected = ('chr9_14772_c' => 1,
                'chr22_8102_c' => 1,
                'chr11_1369_c' => 1,
                'chr22_8101_c' => 0,
);
is_deeply(\%result, \%expected,
          'parse_validation_output correctly parsed MiRdup validation output');
