package miRkwood::FilterBED;

# ABSTRACT : methods to filter a BED file

use strict;
use warnings;

use File::Spec;
use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::Paths;

=method filterBEDfile

Method to filter a BED from the reads corresponding to known 
miRNAS and from other features (CDS, rRNA, tRNA, snoRNA).
Filtered reads are written in auxiliary files.
Input : - a BED file
        - a species name (to find the annotations files)

=cut
sub filterBEDfile {
    my @args = @_;
    my $class   = shift @args;
    my $bed_file = shift @args;
    my $species = shift @args;
    my $basename = "";
    
    my $data_path = miRkwood::Paths->get_data_path();
    my $mirbase_file = File::Spec->catfile( $data_path, "miRBase/${species}_miRBase.gff3" );
    my $annotation_file = File::Spec->catfile( $data_path, "annotations/${species}_annotations.gff" );   
    
    if ( $bed_file =~ /(.*)\.bed/ ){
        $basename = $1;
    }
    
    my $mirna_reads = "${basename}_miRNAs.bed";
    my $features_reads = "${basename}_other.bed"; 
    my $filtered_bed = "${basename}_filtered.bed";   
    my $tmp_bed = "${basename}_tmp.bed";  
    
    if ( -r $mirbase_file ) {
        # Create a file with known miRNAs
        store_overlapping_reads( $bed_file, $mirbase_file, $mirna_reads, "-f 1"); 
        
        # Delete reads corresponding to known miRNAs from the BED
        store_non_overlapping_reads( $bed_file, $mirbase_file, $tmp_bed);
        
        if ( -r $annotation_file ){
            # Create a file with other features
            store_overlapping_reads( $tmp_bed, $annotation_file, $features_reads, "");
              
            # Delete reads corresponding to other featuress from the BED
            store_non_overlapping_reads( $tmp_bed, $annotation_file, $filtered_bed); 
        } 
        else{
            debug( "WARNING : no annotation file found.", miRkwood->DEBUG() );    
        }  
        
    }
    else{
        debug( "WARNING : no miRNA file found.", miRkwood->DEBUG() ); 
        
        if ( -r $annotation_file ){
            # Create a file with other features
            store_overlapping_reads( $bed_file, $annotation_file, $features_reads, "");
              
            # Delete reads corresponding to other featuress from the BED
            store_non_overlapping_reads( $bed_file, $annotation_file, $filtered_bed); 
        } 
        else{
            debug( "WARNING : no annotation file found.", miRkwood->DEBUG() );    
        }
                                      
    }
    

    
    unlink $tmp_bed;
    
}

=method store_overlapping_reads

Method to write reads overlapping with the given reference file
in an output file.
Input : - a BED file
        - a file with the features to search for
        - an output file

=cut
sub store_overlapping_reads {
    my @args = @_;
    my $bed_file = shift @args;
    my $referenceFile = shift @args;
    my $outputFile = shift @args;
    my $additionalArgs = shift @args;
    
    my $job = "intersectBed -a $bed_file -b $referenceFile -s -wa -wb $additionalArgs > $outputFile";
    system($job);
    debug( "$outputFile created", miRkwood->DEBUG() );
    
}

=method store_non_overlapping_reads

Method to write reads non overlapping with the given reference file
in an output file.
Input : - a BED file
        - a file with the features to search for
        - an output file

=cut
sub store_non_overlapping_reads {
    my @args = @_;
    my $bed_file = shift @args;
    my $referenceFile = shift @args;
    my $outputFile = shift @args;
    
    my $job = "intersectBed -a $bed_file -b $referenceFile -s -v > $outputFile";
    system($job);    
    debug( "$outputFile created", miRkwood->DEBUG() );
    
}

1;
