#!/usr/bin/perl -w

# PODNAME: export_csv.pl
# ABSTRACT: Export the CSV for a given job

use warnings;
use strict;

use Pod::Usage;
use FindBin;

BEGIN {
    use lib File::Spec->catdir( $FindBin::Bin, '..', 'lib' );
    use PipelineMiRNA;
    use PipelineMiRNA::Results;
}

pod2usage("$0: No job identifier given.") if ( @ARGV == 0 );

my ($jobId) = @ARGV;

my %results = PipelineMiRNA::Results->get_structure_for_jobID($jobId);
my $csv_file = "$jobId.csv";
open( my $CSV_FILE, '>', $csv_file )
  or die("Cannot open file $csv_file: $!");
print {$CSV_FILE} PipelineMiRNA::Results->resultstruct2csv(\%results)
  or die("Cannot write in file $csv_file: $!");
close($CSV_FILE)
  or die("Cannot close file $csv_file: $!");


__END__

=head1 SYNOPSIS

export_csv.pl [JOB ID]

=head1 DESCRIPTION

Export the CSV for a given miRkwood job

=cut]