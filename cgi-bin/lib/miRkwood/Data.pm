package miRkwood::Data;

# ABSTRACT: Managing external data

use strict;
use warnings;

use miRkwood::Paths;

=method get_mirbase_file

Return the path to the Mirbase file

=cut

sub get_mirbase_file {
    my @args = @_;
    return File::Spec->catfile( miRkwood::Paths->get_data_path(), 'MirbaseFile.txt' );
}


=method get_matrix_file

Return the path to the matrix file

=cut

sub get_matrix_file {
    my @args = @_;
    return File::Spec->catfile( miRkwood::Paths->get_data_path(), 'matrix' );
}

=method get_mirdup_data_path

Return the path to the mirdup data directory

=cut

sub get_mirdup_data_path {
    my @args = @_;
    return File::Spec->catdir( miRkwood::Paths->get_data_path(), 'mirdup');
}

=method get_mirdup_model_name

Return the name of the mirdup model

=cut

sub get_mirdup_model_name {
    my @args = @_;
    return 'plant.model';
}

1;