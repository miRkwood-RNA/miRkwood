#!/usr/bin/perl

use strict;
use warnings;

use Test::More qw/no_plan/;
use Test::File;

use FindBin;
require File::Spec->catfile( $FindBin::Bin, 'Funcs.pl' );

BEGIN {
    use_ok('PipelineMiRNA::Utils');
}
require_ok('PipelineMiRNA::Utils');

my $fastaFile1 = input_file('fasta1.fa');
file_exists_ok($fastaFile1);

open( my $INPUT_FH, '<', $fastaFile1 ) or die "Fail to open: $!";
ok( my %tab = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH),
    'Can call parse_multi_fasta()' );
close $INPUT_FH;

my %expected = (
    '>fasta11' => 'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC',
    '>fasta12' =>
'ACTGAGATCGCAACTAATTTATTTATTCGCTCGTATAATGTATACATTAGATAGAGGCCTAGCCTCTTAGTCGAAAAGCCC',
);
is_deeply( \%tab, \%expected, 'FASTA parsing with parse_multi_fasta is ok' );

my $fastaFile2 = input_file('fasta_long_header.fa');
file_exists_ok($fastaFile2);
open( my $INPUT_FH2, '<', $fastaFile2 ) or die "Fail to open: $!";
ok( my %tab2 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH2),
    'Can call parse_multi_fasta()' );
close $INPUT_FH2;
my %expected2 =
  ( '>fasta1' => 'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC', );
is_deeply( \%tab2, \%expected2,
           'Parsing FASTA with long header with parse_multi_fasta ok' );

my $fastaFile3 = input_file('fasta_with_pipes.fa');
file_exists_ok($fastaFile3);
open( my $INPUT_FH3, '<', $fastaFile3 ) or die "Fail to open: $!";
ok( my %tab3 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH3),
    'Can call parse_multi_fasta()' );
close $INPUT_FH3;
my %expected3 = ( '>gi-425626932-gb-JX648278.1-' =>
                  'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC', );

is_deeply( \%tab3, \%expected3,
           'Parsing FASTA with pipes using parse_multi_fasta ok' );

my $fastaFile4 = input_file('fasta2.fa');
file_exists_ok($fastaFile4);
open( my $INPUT_FH4, '<', $fastaFile4 ) or die "Fail to open: $!";
ok( my %tab4 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH4),
    'Can call parse_multi_fasta()' );
close $INPUT_FH4;
my %expected4 = ( '>contig15750' =>
                  'aatgagtaagataaattgctaattaaatgcgacgagaggttcatacatgaagagaagagtgctcttattatgtagccaaggatgaattgcctaatgacagctcaagtcgtttaaaaaacgactctttgttggtttattaggcgttcatttcttgactgacttaatcggctttttttcatcatgttagatcttctcaacttgttacgagcatatcgttcaatattttcatagtcttcttgtaatatgactttgtcaagtcatttcatatagctacttatgtgtagctattattgtcataattattatatagattatatacttaaagagagacttgtaagggatttaagatgtttagataatcatgtaacattcttgtcaagttatgatcaagcattat',
                  '>contig15916' =>
                  'aaaaaacctcacatacagcccccgtatctctctctctctataattgataggctattttcttctctctctagaaatgagcttacatggcatgcagatccattgcttatttataggtatagatacagcagatatatattatttattcatatatgtgtatcgaggtatcggaagaagaaattttcattgttacggcggttttctgattcgcttggtgcaggtcgggaacggcttggccgacggtttcatatttgtctccactgtgtgaaacctcgtagcttgagtactgtcctgccttgcatcaactgaatctgaaccgatgtaaatgatctgtgaccggtgtaggagaattggatgaatattgttggagat'
                  );
is_deeply( \%tab4, \%expected4,
           'Parsing FASTA with pipes using parse_multi_fasta ok' );

##################################################################
# Testing find_matching_count()
is( PipelineMiRNA::Utils::find_matching_count('()..'), 1, 'one - matching left' );
is( PipelineMiRNA::Utils::find_matching_count('..()'), 1, 'one - matching right' );
is( PipelineMiRNA::Utils::find_matching_count('.().'), 1, 'one - matching middle' );
is( PipelineMiRNA::Utils::find_matching_count('(..)'), 1, 'one - matching across' );
is( PipelineMiRNA::Utils::find_matching_count('(.().)'), 2, 'two - matching across' );

##################################################################
# Testing make_loop()

ok( my @res5 = PipelineMiRNA::Utils::make_loop('123'),
    'Can call make_loop() with 3 elements');

is_deeply( \@res5, [ [], ['1'], [' ', '2'], ['3'], [] ], 'make_loop() with 3 elements ok' );

ok( my @res6 = PipelineMiRNA::Utils::make_loop('1234'),
    'Can call make_loop() with 4 elements');

is_deeply( \@res6, [ ['1'], [' ', '2'], [' '], [' ', '3'], ['4'] ], 'make_loop() with 4 elements ok' );

ok( my @res7 = PipelineMiRNA::Utils::make_loop('12345'),
    'Can call make_loop() with 5 elements');
is_deeply( \@res7, [ ['1'], [' ', '2'], [' ', '3'], [' ', '4'], ['5'] ], 'make_loop() with 5 elements ok' );

ok( my @res8 = PipelineMiRNA::Utils::make_loop('123456'),
    'Can call make_loop() with 5 elements');
my @expected8 = [ ['12'], [' ', ' ', '3'], [' ', ' '], [' ', ' ', '4'], ['65'] ];
is_deeply( \@res8, @expected8, 'make_loop() with 6 elements ok' );

my $sequence1 = "cuuauuauguagccaaggaugaauugccuaaugacagcucaagucguuuaaaaaacgacucuuuguugguuuauuaggcguucauuucuugacugacuuaaucgg";
my $vienna1   = "((.((((.((((.(((((((((((.(((((((((((((...((((((((...))))))))....))))..)))))))))))))).)))))).)).)).)))).))";

ok( my $result9 = PipelineMiRNA::Utils::make_ASCII_viz($sequence1, $vienna1),
    'Can call make_ASCII_viz()');
my $expected9 = "  u    u  -  c      -     u         --    uca-        
cu auua gu ag caagga ugaau gccuaauga  cagc    agucguuua
|| |||| || || |||||| ||||| |||||||||  ||||    |||||||| a
gg uaau ca uc guucuu acuug cggauuauu  guug    ucagcaaaa
  c    u  g  a      u     -         ug    uuuc        
";
is( $result9, $expected9, 'make_ASCII_viz returns a correct hairpin');


my $sequence2 = "gucuccacugugugaaaccucguagcuugaguacuguccugccuugcaucaacugaaucugaaccgauguaaaugaucugugaccgguguaggagaauuggaugaauauuguuggagau";
my $vienna2   = "((((((((.((((..(.((...(..((((..((((((((....(((((((...............))))))).......).)).)))))))))..)...)).)..)))).).)))))))";
ok( my $result10 = PipelineMiRNA::Utils::make_ASCII_viz($sequence2, $vienna2),
    'Can call make_ASCII_viz()');
my $expected10 = "       - u    ga a  ucg ag    ag     -  - ugcc---       aacuga
gucucca c gugu  a cc   u  cuug  uacug uc c       uugcauc      a
||||||| | ||||  | ||   |  ||||  ||||| || |       |||||||      u
uagaggu g uaua  u gg   a  ggau  guggc ag g       aauguag      c
       u u    ag a  uua ga    --     c  u ucuagua       ccaagu
";
is( $result10, $expected10, 'make_ASCII_viz returns a correct hairpin with big loop');

my $sequence3 = "ccgacgguuucauauuugucuccacugugugaaaccucguagcuugaguacuguccugccuugcaucaacugaaucugaaccgauguaaaugaucugugaccgguguaggagaauuggaugaauauuguugg";
my $vienna3   = "(((((((((((((.....(((((((((.((.(....((((.................................................))))....).))))))...))))).....))))).))))))))";
ok( my $result11  = PipelineMiRNA::Utils::make_ASCII_viz($sequence3, $vienna3),
    'Can call make_ASCII_viz()');
    my $expected11 = "        -     auuug     ---    u  g aacc    agcuugaguacuguccugccuug
ccgacggu uucau     ucucc   acug gu a    ucgu                       c
|||||||| |||||     |||||   |||| || |    ||||                       a
gguuguua aagua     agagg   uggc ca u    agua                       u
        u     gguua     aug    -  g gucu    aauguagccaagucuaagucaac
";
is( $result11, $expected11, 'make_ASCII_viz returns a correct hairpin with giant loop');

my $seq_with_T = "ATGCATGC";
my $seq_with_U = 'AUGCAUGC';
my $result12 = "";
open my ($input_seq_fh), '<', \$seq_with_T;
open my ($result12_fh),  '>', \$result12;
PipelineMiRNA::Utils::rewrite_fasta_with_TU('U', $input_seq_fh, $result12_fh);
close $input_seq_fh;
close $result12_fh;
is ( $result12, $seq_with_U, 'rewrite_fasta_with_TU correctly replace T with U');

open $input_seq_fh, '<', \$seq_with_U;
open $result12_fh,  '>', \$result12;
PipelineMiRNA::Utils::rewrite_fasta_with_TU('T', $input_seq_fh, $result12_fh);
is ( $result12, $seq_with_T, 'rewrite_fasta_with_TU correctly replace U with T');


my $top1 = "   g  auauu----       a                 c ccuc          -     u  a      c       caccuuucuagcagaucaacaaugaauuuuguggaauagauguugga";
my $left1 = 53;
my $right1 = 21;

ok( my ($true_left1, $size1) = PipelineMiRNA::Utils::compute_mature_boundaries($left1, $right1, $top1),
    'Can call compute_mature_boundaries');
is_deeply([$true_left1, $size1], [57, 22], 'compute_mature_boundaries ok with gaps before start');



my $hairpin1 = <<"END";
   a--      ca  cac-      g u       a-  g   aaaaaa     g
ggu   gagacu  uc    ccggca c cuguaau  gg acu      gugau a
|||   ||||||  ||    |||||| | |||||||  || |||      ||||| 
cca   uucugg  ag    ggcugu g gguauug  cc uga      uacua u
   gca      cg  acua      g u       cg  g   g-----     a
END
ok (my $hairpin_with_mature = PipelineMiRNA::Utils::make_hairpin_with_mature($hairpin1, 56, 75, 105),
    'Can call make_hairpin_with_mature on edge case');

my $expected13 = <<"END";
   a--      ca  cac-      g u       a-  g   aaaaaa     g
ggu   gagacu  uc    ccggca c cuguaau  gg acu      gugau a
|||   ||||||  ||    |||||| | |||||||  || |||      ||||| 
cca   uucugg  ag    ggcugu g g<span class="mature">guauug  cc uga      uacua</span> u
   gca      cg  acua      g u <span class="mature">      cg  g   g-----     </span>a
END

is ($hairpin_with_mature, $expected13,
    'make_hairpin_with_mature ok on edge case');

# Bug 2014-03-13
my $hairpin2 = <<"END";
    gg   c  --------      uaaaa   aac    c    u
agug  ugg gc        gggagc     uca   ucua gcug u
||||  ||| ||        ||||||     |||   |||| |||| 
ucac  acc cg        cccucg     agu   aggu uggc u
    aa   -  uuuaacuu      uggga   c--    -    a
';
END
ok (my $hairpin_with_mature2 = PipelineMiRNA::Utils::make_hairpin_with_mature($hairpin2, 45, 66, 84),
    'Can call make_hairpin_with_mature on other edge case');

my $expected2 = <<"END";
    gg   c  --------      uaaaa   aac    c    u
agug  ugg gc        gggagc     uca   ucua gcug u
||||  ||| ||        ||||||     |||   |||| |||| 
ucac  acc cg       <span class="mature"> cccucg     agu   aggu u</span>ggc u
    aa   -  uuuaacu<span class="mature">u      uggga   c--    - </span>   a
END

is ($hairpin_with_mature2, $expected2,
    'make_hairpin_with_mature ok on other edge case');

my @input_fasta_ok = (
    '>Seq
aatgagtaagataaa
',
    '>Seq1
aatgagtaagataaa
>Seq2
agctagctgatgcga
',
    '>zma-mir167i_mi0001823
acuucgcuggugugagagcuugaagcugcccggccucccaagugcucgaucgguggcgcuucaccagaucauguugcagcuucacucucucgcaaccagcgaa
'
);

foreach my $input_fasta (@input_fasta_ok) {
    ok( my $fasta_ok_res = PipelineMiRNA::Utils::is_fasta($input_fasta),
        'Can call is_fasta on a correct FASTA' );
    is( $fasta_ok_res, 1, 'is_fasta detects a correct FASTA' );
}

my @input_fasta_wrong = (
    '>Seq
aatgagtnaagataaa
',
    '>Seq1
aatgagtaagataaa
>Seq2
agctangctgatgcga
',
    'Toto'
);
foreach my $input_fasta (@input_fasta_wrong) {
    is( PipelineMiRNA::Utils::is_fasta($input_fasta),
        0, 'is_fasta detects a wrong FASTA' );
}

my @reverse_complement_values = (
    [ "AA", "UU" ],
    [ "UU", "AA" ],
    [ "CC", "GG" ],
    [ "GG", "CC" ],
    [ "AC", "GU" ],
    [ "AU", "AU" ],
    [ "AG", "CU" ],
    [ $sequence1,  'ccgauuaagucagucaagaaaugaacgccuaauaaaccaacaaagagucguuuuuuaaacgacuugagcugucauuaggcaauucauccuuggcuacauaauaag'],
);

foreach my $couple (@reverse_complement_values) {
    ok(
        my $reverse_complement_res =
          PipelineMiRNA::Utils::reverse_complement( @{$couple}[0] ),
        'Can call reverse_complement'
    );
    is( $reverse_complement_res, @{$couple}[1],
        "reverse_complement correctly works @{$couple}[0] --> @{$couple}[1]" );
    is(
        PipelineMiRNA::Utils::reverse_complement($reverse_complement_res),
        @{$couple}[0],
"reverse_complement correctly works back @{$couple}[0] --> @{$couple}[1] --> @{$couple}[0]"
    );
}

ok(
    my $compute_mfei_res =
      PipelineMiRNA::Utils::compute_mfei( $sequence1, -34.30 ),
    'Can call compute_mfei'
);
is( $compute_mfei_res, -0.879487179487179,
    'compute_mfei returns the expected value' );

ok(
    my $compute_amfe_res =
      PipelineMiRNA::Utils::compute_amfe( $sequence1, -34.30 ),
    'Can call compute_amfe'
);

is(
    $compute_amfe_res,
    ( -34.3 / 105 ) * 100,
    'compute_amfe returns the expected value (exact)'
);

is( $compute_amfe_res, -32.6666666666667,
    'compute_amfe returns the expected value (rounded)' );

ok( my $compute_gc_res = PipelineMiRNA::Utils::compute_gc_content($sequence1),
    'Can call compute_gc_content' );
is(
    $compute_gc_res,
    ( 39 / 105 * 100 ),
    'compute_gc_content returns the expected value (exact)'
);

is( $compute_gc_res, 37.142857142857146,
    'compute_gc_content returns the expected value (rounded)' );

ok(
    my $make_mirbase_link_output =
      PipelineMiRNA::Utils::make_mirbase_link('ABCDE'),
    'can call make_mirbase_link()'
);
is(
    $make_mirbase_link_output,
    'http://mirbase.org/cgi-bin/mirna_entry.pl?acc=ABCDE',
    'make_mirbase_link returns expected value'
);

my @opposite_strand_values = (
    [ [ 3,   7,   10 ],  [ 2,  6 ] ],
    [ [ 3,   7,   20 ],  [ 12, 16 ] ],
    [ [ 2,   4,   9 ],   [ 4,  6 ] ],
    [ [ 143, 366, 369 ], [ 2,  225 ] ],
);

foreach my $couple (@opposite_strand_values) {
    my @arrary = @{$couple};
    ok(
        my $opposite_strand_res =
          PipelineMiRNA::Utils::get_position_from_opposite_strand(
            @{ @{$couple}[0] }
          ),
        'Can call get_position_from_opposite_strand'
    );
    is_deeply( $opposite_strand_res, @{$couple}[1],
"get_position_from_opposite_strand (@{@{$couple}[0]}) --> (@{@{$couple}[1]}) ok"
    );
}
