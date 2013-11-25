use PipelineMiRNA::Components;

my ( $idirData, $idirJob, $iplant ) = @ARGV;

PipelineMiRNA::Components::filter_CDS( $idirData, $idirJob, $iplant );
