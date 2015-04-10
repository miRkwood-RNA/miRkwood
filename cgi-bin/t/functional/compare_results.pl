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

compare_folders( $results_folder_1, $results_folder_2 );



########## Functions
=method compare_folders

  Compare recursively folders and files.
  If files are yaml, deserialize them and use 
  function compare_2_hash.

=cut
sub compare_folders {
    my ($folder_1, $folder_2) = @_;
    
    ### List all sub-folders and files in each folder
    my ($list_files_1, $list_subfolders_1) = list_files_and_folders( $folder_1 );
    my ($list_files_2, $list_subfolders_2) = list_files_and_folders( $folder_2 );

    ### Compare files
    my $files_1 = join(' ', sort( @{$list_files_1} ) );
    my $files_2 = join(' ', sort( @{$list_files_2} ) );
    if ( $files_1 ne $files_2 ){
        print "Lists of files in $folder_1 and $folder_2 are different.\n";
        print "---$files_1\n";
        print "---$files_2\n";
    }
    else{
        foreach ( @{$list_files_1} ){
            if ( $_ eq 'outBlast.txt' or $_ eq 'pvalue.txt' or $_ eq 'basic_candidates.yml' or $_ eq 'basic_known_candidates.yml' ){
                # don't compare these files
            }
            elsif ( /[.]yml/ ){
                my $yml_file_1 = "$folder_1/$_";
                my $yml_file_2 = "$folder_2/$_";

                my %attributes_1 = YAML::XS::LoadFile($yml_file_1);
                my %attributes_2 = YAML::XS::LoadFile($yml_file_2);

                compare_2_hash( \%attributes_1, \%attributes_2, $_);
            }
            elsif ( $_ !~/[.]svn|[.]log|[.]cfg|[.]png|[.]html|miRdupOutput[.]txt$/ ) {
                open CMD, "diff $folder_1/$_ $folder_2/$_ |";
                while ( my $line = <CMD> ){
                    if ( $line ne /^\n$/ ){
                        print $line;
                    }
                }
                close CMD;
            }
        }
    }

    ### Compare subfolders
    my $subfolders_1 = join(' ', sort( @{$list_subfolders_1} ) );
    my $subfolders_2 = join(' ', sort( @{$list_subfolders_2} ) );
    if ( $subfolders_1 ne $subfolders_2 ){
        print "Lists of subfolders in $folder_1 and $folder_2 are different.\n";
        print "---$subfolders_1\n";
        print "---$subfolders_2\n";
    }
    else{
        foreach ( @{$list_subfolders_1} ){
            compare_folders( "$folder_1/$_", "$folder_2/$_" );
        }
    }

}

=method list_files_and_folders

  List files and subfolders in the folder
  given in parameter.
  Output : 2 reference arrays, one for files
  and one for folders.
  
=cut
sub list_files_and_folders {
    my ($folder) = @_;
    my $list_files = [];
    my $list_folders = [];

    opendir (my $DIR, $folder) or die "ERROR : cannot open directory $folder : $!";
    while ( readdir $DIR ){
        if ( $_ ne "." && $_ ne ".." ){
            if ( -f "$folder/$_" ){
                push @{$list_files}, $_;
            }
            elsif ( -d "$folder/$_" ){
                push @{$list_folders}, $_;
            }
        }
    }
    closedir $DIR;

    return ($list_files, $list_folders);

}


=method compare_2_hash

  Compares 2 hash tables given by reference.
  If the values are hash, recursive call to compare_2_hash

=cut
sub compare_2_hash {
    my ($hash1, $hash2, $file) = @_;

    my $keys_1 = join(" ", sort (keys%{$hash1}) );
    my $keys_2 = join(" ", sort (keys%{$hash2}) );

    if ( $keys_1 ne $keys_2 ){
        print "(File $file) Hash tables don't have the same features.\n";
        print "--- $keys_1\n";
        print "--- $keys_2\n";
    }
    else{
        foreach ( sort (keys%$hash1) ) {
            if ( $_ ne "image" ){
                if ( eval { scalar(keys%{$hash1->{$_}}) >= 0 } ){   # value is a hash
                    compare_2_hash( $hash1->{$_}, $hash2->{$_}, $file );
                }
                elsif ( eval { scalar(keys@{$hash1->{$_}}) >= 0 } ){    # value is an array
                    my $value_1 = join(" ", sort( keys@{$hash1->{$_}} ) );
                    my $value_2 = join(" ", sort( keys@{$hash2->{$_}} ) );
                    if ( $value_1 ne $value_2 ){
                        print "(File $file) Arrays don't have the same contents.\n";
                        print "--- $value_1\n";
                        print "--- $value_2\n";
                    }
                }
                elsif ( $$hash1{$_} ne $$hash2{$_} ){   # value is a scalar
                    print "(File $file) Different values for key $_.\n";
                }
            }
        }
    }
}
