use miRkwood::Maskers;
use miRkwood::Programs;

my ( $idirData, $idirJob, $iplant ) = @ARGV;

miRkwood::Programs::init_programs();
miRkwood::Maskers::get_coding_region_masking_information( $idirData, $idirJob, $iplant );
