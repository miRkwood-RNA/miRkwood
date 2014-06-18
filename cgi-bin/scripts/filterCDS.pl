use PipelineMiRNA::Maskers;
use PipelineMiRNA::Programs;

my ( $idirData, $idirJob, $iplant ) = @ARGV;

PipelineMiRNA::Programs::init_programs();
PipelineMiRNA::Maskers::get_coding_region_masking_information( $idirData, $idirJob, $iplant );
