use PipelineMiRNA::Components;
use PipelineMiRNA::Programs;

my ( $idirData, $idirJob, $iplant ) = @ARGV;

PipelineMiRNA::Programs::init_programs();
PipelineMiRNA::Components::filter_CDS( $idirData, $idirJob, $iplant );
