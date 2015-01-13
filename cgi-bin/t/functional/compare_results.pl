#!/usr/bin/perl -w
use strict;
use YAML::XS;
use Data::Dumper;

######################################################################
# Script to compare the contents of 2 mirkwood output directories.
# Since the order in which data are serialized can vary, this script
# deserializes the given files and compares them in depth.
######################################################################


my $results_folder_1 = $ARGV[0];
my $results_folder_2 = $ARGV[1];


########## Compare the yml files

my $yml_folder_1 = "$results_folder_1/candidates";
my $yml_folder_2 = "$results_folder_2/candidates";

my @dir1;
my @dir2;

##### Get yml folders contents
opendir(D, $yml_folder_1) or die "ERROR : cannot opendir $yml_folder_1 : $!";
while (readdir D) {
    if ( $_ ne "." && $_ ne ".." && $_ ne "basic_candidates.yml" && /.yml$/ ){
        push @dir1, $_;
    }
}
closedir D;

opendir(D, $yml_folder_2) or die "ERROR : cannot opendir $yml_folder_2 : $!";
while (readdir D) {
    if ( $_ ne "." && $_ ne ".." && $_ ne "basic_candidates.yml" && /.yml$/ ){
        push @dir2, $_;
    }
}
closedir D;

##### Compare yml folders contents
my $content1 = join(" ", sort(@dir1));
my $content2 = join(" ", sort(@dir2));


if ( $content1 ne $content2 ){
    print "Lists of files in both folders are different.\n";
    print "---$content1\n";
    print "---$content2\n";
    exit 1;
}

foreach ( sort(@dir1) ){
    
    my $yml_file_1 = "$yml_folder_1/$_";
    my $yml_file_2 = "$yml_folder_2/$_"; 
    
    my %attributes_1 = YAML::XS::LoadFile($yml_file_1);
    my %attributes_2 = YAML::XS::LoadFile($yml_file_2);
    
    compare_2_hash( \%attributes_1, \%attributes_2, $_); 
}


########## Compare the other files
my $excludes = "--exclude=.svn --exclude=pvalue.txt --exclude=outBlast.txt --exclude=*miRdupOutput.txt --exclude=*.log --exclude=*.cfg --exclude=*.png --exclude=*.yml --exclude=*.html";
open CMD, "diff $excludes $results_folder_1 $results_folder_2 | grep -v 'Binary' | grep -v 'Les sous-répertoires.*sont identiques' | grep -v 'Common subdirectories:' |";
while ( <CMD> ){
    if ( ! /^\n$/ ){
        print $_;
    }
}
close CMD;



########## Functions
=method compare_2_hash
Compares 2 hash tables given by reference.
If the values are hash, recursive call to compare_2_hash
=cut
sub compare_2_hash {
    my ($hash1, $hash2) = @_;
    
    my $keys_1 = join(" ", sort (keys%$hash1) );
    my $keys_2 = join(" ", sort (keys%$hash2) );
    
    if ( $keys_1 ne $keys_2 ){
        print "Hash tables don't have the same features.\n";
        print "--- $keys_1\n";
        print "--- $keys_2\n";        
    }
    else{
        foreach ( sort (keys%$hash1) ) {
            if ( $_ ne "image" ){
                if ( eval { scalar(keys%{$hash1->{$_}}) >= 0 } ){   # value is a hash
                    compare_2_hash( $hash1->{$_}, $hash2->{$_} );
                }
                elsif ( eval { scalar(keys@{$hash1->{$_}}) >= 0 } ){    # value is an array
                    my $value_1 = join(" ", sort( keys@{$hash1->{$_}} ) );
                    my $value_2 = join(" ", sort( keys@{$hash2->{$_}} ) );
                    
                    if ( $value_1 ne $value_2 ){
                        print "Arrays don't have the same contents.\n";
                        print "--- $value_1\n";
                        print "--- $value_2\n";        
                    }                   
                }
                elsif ( $$hash1{$_} ne $$hash2{$_} ){   # value is a scalar
                    print "Different values for key $_.\n";
                }
            }
        }
    }  
}
