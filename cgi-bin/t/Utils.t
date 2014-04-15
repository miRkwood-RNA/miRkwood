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

my $fastaFile1 = input_file('Utils.fasta1.fa');
file_exists_ok($fastaFile1);

open( my $INPUT_FH, '<', $fastaFile1 ) or die "Fail to open: $!";
ok( my @fasta_array = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH),
    'Can call parse_multi_fasta()' );
close $INPUT_FH;

my @expected11 = ('>fasta11', 'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC');
my @expected12 = ('>fasta12', 'ACTGAGATCGCAACTAATTTATTTATTCGCTCGTATAATGTATACATTAGATAGAGGCCTAGCCTCTTAGTCGAAAAGCCC');
my @expected1 = (\@expected11, \@expected12);
is_deeply( \@fasta_array, \@expected1, 'FASTA parsing with parse_multi_fasta is ok' );

my $fastaFile2 = input_file('Utils.fasta_long_header.fa');
file_exists_ok($fastaFile2);
open( my $INPUT_FH2, '<', $fastaFile2 ) or die "Fail to open: $!";
ok( my @fasta_array2 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH2),
    'Can call parse_multi_fasta()' );
close $INPUT_FH2;

my @expected21 =
  ( '>fasta1_Drosophila_simulans_strain_Eden32', 'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC' );
my @expected2 = (\@expected21);
is_deeply( \@fasta_array2, \@expected2,
           'Parsing FASTA with long header with parse_multi_fasta ok' );

my $fastaFile3 = input_file('Utils.fasta_with_pipes.fa');
file_exists_ok($fastaFile3);
open( my $INPUT_FH3, '<', $fastaFile3 ) or die "Fail to open: $!";
ok( my @fasta_array3 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH3),
    'Can call parse_multi_fasta()' );
close $INPUT_FH3;
my @expected31 = ( '>gi|425626932|gb|JX648278.1|_Drosophila_simulans_strain_Eden32',
                   'AATGTGCCAATCCCAATGTTAACCAAAAACTAAAAAAGTGAAACGAACATTGTC');
my @expected3 = (\@expected31);
is_deeply( \@fasta_array3, \@expected3,
           'Parsing FASTA with pipes using parse_multi_fasta ok' );

my $fastaFile4 = input_file('Utils.fasta2.fa');
file_exists_ok($fastaFile4);
open( my $INPUT_FH4, '<', $fastaFile4 ) or die "Fail to open: $!";
ok( my @fasta_array4 = PipelineMiRNA::Utils::parse_multi_fasta($INPUT_FH4),
    'Can call parse_multi_fasta()' );
close $INPUT_FH4;
my @expected41 = ('>contig15750',
                  'aatgagtaagataaattgctaattaaatgcgacgagaggttcatacatgaagagaagagtgctcttattatgtagccaaggatgaattgcctaatgacagctcaagtcgtttaaaaaacgactctttgttggtttattaggcgttcatttcttgactgacttaatcggctttttttcatcatgttagatcttctcaacttgttacgagcatatcgttcaatattttcatagtcttcttgtaatatgactttgtcaagtcatttcatatagctacttatgtgtagctattattgtcataattattatatagattatatacttaaagagagacttgtaagggatttaagatgtttagataatcatgtaacattcttgtcaagttatgatcaagcattat'
                  );
my @expected42 = ('>contig15916',
                  'aaaaaacctcacatacagcccccgtatctctctctctctataattgataggctattttcttctctctctagaaatgagcttacatggcatgcagatccattgcttatttataggtatagatacagcagatatatattatttattcatatatgtgtatcgaggtatcggaagaagaaattttcattgttacggcggttttctgattcgcttggtgcaggtcgggaacggcttggccgacggtttcatatttgtctccactgtgtgaaacctcgtagcttgagtactgtcctgccttgcatcaactgaatctgaaccgatgtaaatgatctgtgaccggtgtaggagaattggatgaatattgttggagat'
                  );

my @expected4 = (\@expected41, \@expected42);
is_deeply( \@fasta_array4, \@expected4,
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
ok (my $make_hairpin_with_mature_html_1 = PipelineMiRNA::Utils::make_hairpin_with_mature($hairpin1, 56, 75, 105, 'html'),
    'Can call make_hairpin_with_mature in HTML on edge case');

my $make_hairpin_with_mature_html_1_out = slurp_file( input_file('Utils.make_hairpin_with_mature.1.html.txt') );

is ($make_hairpin_with_mature_html_1, $make_hairpin_with_mature_html_1_out,
    'make_hairpin_with_mature in HTML ok on edge case');

ok (my $make_hairpin_with_mature_ascii_1 = PipelineMiRNA::Utils::make_hairpin_with_mature($hairpin1, 56, 75, 105, 'ascii'),
    'Can call make_hairpin_with_mature in ASCII on edge case');

my $make_hairpin_with_mature_ascii_1_out = slurp_file( input_file('Utils.make_hairpin_with_mature.1.ascii.txt') );

is ($make_hairpin_with_mature_ascii_1, $make_hairpin_with_mature_ascii_1_out,
    'make_hairpin_with_mature in ASCII ok on edge case');


# Bug 2014-03-13
my $hairpin2 = <<"END";
    gg   c  --------      uaaaa   aac    c    u
agug  ugg gc        gggagc     uca   ucua gcug u
||||  ||| ||        ||||||     |||   |||| |||| 
ucac  acc cg        cccucg     agu   aggu uggc u
    aa   -  uuuaacuu      uggga   c--    -    a
';
END
ok (my $make_hairpin_with_mature_html_2 = PipelineMiRNA::Utils::make_hairpin_with_mature($hairpin2, 45, 66, 84, 'html'),
    'Can call make_hairpin_with_mature in HTML mode on other edge case');
my $make_hairpin_with_mature_html_2_out = slurp_file( input_file('Utils.make_hairpin_with_mature.2.html.txt') );
is ($make_hairpin_with_mature_html_2, $make_hairpin_with_mature_html_2_out,
    'make_hairpin_with_mature in HTML ok on other edge case');


ok (my $make_hairpin_with_mature_ascii_2 = PipelineMiRNA::Utils::make_hairpin_with_mature($hairpin2, 45, 66, 84, 'ascii'),
    'Can call make_hairpin_with_mature in HTML mode on other edge case');
my $make_hairpin_with_mature_ascii_2_out = slurp_file( input_file('Utils.make_hairpin_with_mature.2.ascii.txt') );
is ($make_hairpin_with_mature_ascii_2, $make_hairpin_with_mature_ascii_2_out,
    'make_hairpin_with_mature in ASCII ok on edge case');


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
    my ($compute_mfei_res, $compute_amfe_res) =
      PipelineMiRNA::Utils::compute_mfei_and_amfe( $sequence1, -34.30 ),
    'Can call compute_mfei_and_amfe'
);
is( $compute_mfei_res, -0.879487179487179,
    'compute_mfei_and_amfe returns the expected MFEI value' );

my $expected_amfe = ( -34.3 / 105 ) * 100;

is( $compute_amfe_res, $expected_amfe,
    'compute_mfei_and_amfe returns the expected AMFE value' );

ok(
    my $compute_amfe_res_bis =
      PipelineMiRNA::Utils::compute_amfe( $sequence1, -34.30 ),
    'Can call compute_amfe'
);

is(
    $compute_amfe_res, $expected_amfe,
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


my $no_gc_sequence = 'UAUAUAUAUAUAUAUAUAUAUAUAUAUAUA';
my $compute_gc_res_no_gc = PipelineMiRNA::Utils::compute_gc_content($no_gc_sequence);
is( $compute_gc_res_no_gc, 0, 'Can call compute_gc_content on no GC sequence' );

ok(
    my ($compute_mfei_res_no_gc, $compute_amfe_res_no_gc) =
      PipelineMiRNA::Utils::compute_mfei_and_amfe( $no_gc_sequence, 60 ),
    'Can call compute_mfei_and_amfe on no GC sequence'
);
is( $compute_mfei_res_no_gc, 0, 'compute_mfei_and_amfe on no GC sequence returns a MFEI of 0' );
is( $compute_amfe_res_no_gc, 200, 'compute_mfei_and_amfe on no GC sequence returns the correct AMFE' );


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

my $cleanup_fasta_sequence = ">sequence1\nTAGCTGATGCATCGAGCGAT\r";
my $cleanup_fasta_sequence_expected = ">sequence1
tagctgatgcatcgagcgat
";
ok( my $cleanup_fasta_sequence_out = PipelineMiRNA::Utils::cleanup_fasta_sequence($cleanup_fasta_sequence),
    'Can call cleanup_fasta_sequence' );
is( $cleanup_fasta_sequence_out, $cleanup_fasta_sequence_expected,
    "cleanup_fasta_sequence correctly cleans the sequence");

##### sanitize_sequence_name ####

my @sanitize_sequence_name_values = (
    [ '>contig15916', '>contig15916' ],
    [ '> contig15916', '>_contig15916' ],
    [ '>random_seq_from_cds__no_72421',
      '>random_seq_from_cds__no_72421' ],
    [ '>random 54%gc, 100000nt',
      '>random_54%gc,_100000nt'],
    [ '>aly-MIR169a MI0014545',
      '>aly-MIR169a_MI0014545'],
    [ '>gi|425626932|gb|JX648278.1| Drosophila simulans strain Eden3',
      '>gi-425626932-gb-JX648278.1-_D'],
    [ '>gi|425626932|gb|JX648278.11| Drosophila simulans strain Eden3',
      '>gi-425626932-gb-JX648278.11-'],
);
foreach my $couple (@sanitize_sequence_name_values) {
    my ($input, $expected) = @{$couple};
    ok(
        my $sanitize_sequence_name_res =
          PipelineMiRNA::Utils::sanitize_sequence_name($input),
        'Can call sanitize_sequence_name'
    );
    is_deeply( $sanitize_sequence_name_res, $expected,
"sanitize_sequence_name ($input) --> ($expected) ok"
    );
}