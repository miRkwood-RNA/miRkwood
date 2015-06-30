package miRkwood::BEDHandler;

# ABSTRACT : methods to handle a BED or BED-like file

use strict;
use warnings;

use File::Spec;
use Log::Message::Simple qw[msg error debug];

use miRkwood;
use miRkwood::Paths;
use miRkwood::Candidate;
use miRkwood::CandidateHandler;


=method filterBEDfile_for_model_organism

Method to filter a BED from the reads corresponding to known 
miRNAS and from other features (CDS, rRNA, tRNA, snoRNA) for a 
model organism.
Filtered reads are written in auxiliary files.
Input : - a BED file
        - a species name (to find the annotations files)

=cut
sub filterBEDfile_for_model_organism {
    my @args = @_;
    my $class              = shift @args;
    my $basename           = shift @args;
    my $bed_file           = shift @args;
    my $species            = shift @args;
    my $filter_CDS         = shift @args;
    my $filter_tRNA_rRNA   = shift @args;
    my $filter_multimapped = shift @args;

    my $data_path = miRkwood::Paths->get_data_path();
    my $mirbase_file = File::Spec->catfile( $data_path, "miRBase/${species}_miRBase.gff3" );
    my $CDS_file = File::Spec->catfile( $data_path, "annotations/${species}_CDS.gff" );
    my $otherRNA_file = File::Spec->catfile( $data_path, "annotations/${species}_otherRNA.gff" );

    if ( $bed_file =~ /[\/\\]([^\/\\]+)[.]bed/ ){
        $basename .= '/'.$1;
    }

    my $mirna_reads       = "${basename}_miRNAs.bed";
    my $CDS_reads         = "${basename}_CDS.bed";
    my $otherRNA_reads    = "${basename}_otherRNA.bed";
    my $multimapped_reads = "${basename}_multimapped.bed";
    my $filtered_bed      = "${basename}_filtered.bed";
    my $tmp_1 = "${basename}_tmp_1.bed";
    my $tmp_2 = "${basename}_tmp_2.bed";
    my $tmp_3 = "${basename}_tmp_3.bed";

    ### Filter out known miRNAs
    if ( -r $mirbase_file ) {
        # Create a file with known miRNAs
        store_overlapping_reads( $bed_file, $mirbase_file, $mirna_reads, '-f 1');

        # Delete reads corresponding to known miRNAs from the BED
        #~ store_non_overlapping_reads( $bed_file, $mirbase_file, $tmp_1);  # un-comment when tests are done
        system("cp $bed_file $tmp_1");  # comment when tests are done
        debug( 'Known miRNAS have been filtered out from BED.', 1 );
    }
    else{
        debug( "WARNING : no miRNA file found for $species.", 1 );
        system("cp $bed_file $tmp_1");
    }

    ### Filter out CDS
    if ( $filter_CDS ){
        if ( -r $CDS_file ){
            # Create a file with CDS
            store_overlapping_reads( $tmp_1, $CDS_file, $CDS_reads, '');

            # Delete reads corresponding to CDS from the BED
            store_non_overlapping_reads( $tmp_1, $CDS_file, $tmp_2);

            debug( 'CDS have been filtered out from BED.', 1 );
        }
        else{
            debug( "WARNING : no CDS file found for $species.", 1 );
            system("cp $tmp_1 $tmp_2");
        }
    }
    else{
        #~ $tmp_2 = $tmp_1;
        system("cp $tmp_1 $tmp_2");
    }

    ### Filter out rRNA, tRNA, snoRNA
    if ( $filter_tRNA_rRNA ){
        if ( -r $otherRNA_file ){
            # Create a file with CDS
            store_overlapping_reads( $tmp_2, $otherRNA_file, $otherRNA_reads, '');

            # Delete reads corresponding to CDS from the BED
            store_non_overlapping_reads( $tmp_2, $otherRNA_file, $tmp_3);

            debug( 'Other RNA have been filtered out from BED.', 1 );
        }
        else{
            debug( "WARNING : no rRNA/tRNA/snoRNA file found for $species.", 1 );
            system("cp $tmp_2 $tmp_3");
        }
    }
    else{
        system("cp $tmp_2 $tmp_3");
    }

    ### Filter out multimapped reads
    if ( $filter_multimapped ){
        filter_multimapped_reads( $tmp_3, $filtered_bed, $multimapped_reads );
        debug( 'Multi_mapped reads have been filtered out from BED.', 1 );
    }
    else{
        rename $tmp_3, $filtered_bed;
    }

    unlink $tmp_1;
    unlink $tmp_2;
    unlink $tmp_3;

    return ($filtered_bed, $mirna_reads);

}

=method filterBEDfile_for_user_sequence

Method to filter a BED from the reads corresponding to known 
miRNAS and from other features (CDS, rRNA, tRNA, snoRNA) for a
sequence given by the user.

=cut
sub filterBEDfile_for_user_sequence {
    my @args = @_;
    my $class              = shift @args;
    my $out_dir            = shift @args;
    my $bed_file           = shift @args;
    my $filter_CDS         = shift @args;
    my $filter_tRNA_rRNA   = shift @args;
    my $filter_multimapped = shift @args;

    ##### Not yet implemented. Call to BlastX, rnammer, tRNAscan-SE... ?

    return;
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

    return;
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

    return;
}


=method filter_multimapped_reads

Method to discard reads with more than 5 mapping positions

=cut
sub filter_multimapped_reads {
    my @args = @_;
    my $bed_file = shift @args;
    my $output_file = shift @args;
    my $discarded_file = shift @args;

    my $max_positions = 5;
    my $counts;
    my @line;

    open (my $BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!\n";

    ### Read the BED file a first time to count how many times
    ### each id is present
    while ( <$BED> ){

        chomp;
        @line = split( /\t/xms );
        if ( ! exists($counts->{$line[3]}) ){
            $counts->{$line[3]} = 0;
        }
        $counts->{$line[3]}++;

    }

    close $BED;

    ### Read the BED file a second time to write on output only the 
    ### non multi-mapped reads
    open ($BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!\n";
    open (my $KEEP, '>', $output_file) or die "ERROR while creating $output_file : $!\n";
    open (my $OUT, '>', $discarded_file) or die "ERROR while creating $discarded_file : $!\n";

    while ( <$BED> ){

        @line = split( /\t/xms );
        if ( $counts->{$line[3]} <= $max_positions ){
            print $KEEP $_;
        }
        else{
            print $OUT $_;
        }

    }

    close $KEEP;
    close $OUT;
    close $BED;

    return;
}

=method count_reads_in_bed_file

Method to count the nb of total reads and unique reads in a BED(-like) file.
Returns first the nb total and then the nb of unique reads.

=cut
sub count_reads_in_bed_file {
    my @args = @_;
    my $bed_file = shift @args;
    my $start = shift @args;
    my $end = shift @args;    

    my @tab;
    my $reads = {};
    my $nb_total_reads = 0;

    if ( -e $bed_file ){
        open(my $BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!";
        while ( <$BED> ){
            @tab = split ( /\t/ );
            if ( ( $start == -1 && $end == -1 ) || ( $tab[1] >= $start  && $tab[2] <= $end ) ) {
                $reads->{"$tab[1]-$tab[2]"} = $tab[4];
            }
        }
        close $BED;

        foreach ( keys%{$reads} ){
            $nb_total_reads += $reads->{$_};
        }

        return ($nb_total_reads, scalar keys%{$reads});
    }
    else {
        return (0, 0);
    }

}

=method make_reads_length_diagramm

  Method to draw a diagramm of reads length in raw text.
  Takes a BED(-like) in parameter.

=cut
sub make_reads_length_diagramm {
    my @args = @_;
    my $bed_file = shift @args;
    my $max_width = 80;

    my $diagramm = '<pre style="font-size:11px;">';
    my %reads_length;

    open(my $BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!";
    while ( <$BED> ){
        my @tab = split ( /\t/ );
        my $length = $tab[2] - $tab[1];
        if ( !defined( $reads_length{ $length } ) ){
            $reads_length{ $length } = 0;
        }
        $reads_length{ $length }++;
    }
    close $BED;

    my $max = 0;
    foreach my $key ( keys%reads_length ){
        if ( $reads_length{$key} > $max ){
            $max = $reads_length{$key};
        }
    }

    my %reads_width;
    foreach my $key ( sort( keys%reads_length ) ){
        my $width = int( $reads_length{$key} * $max_width / $max + 0.5 );
        $reads_width{$key} = $width;
        $diagramm .= "$key nt | ";
        my $i = 0;
        while ( $i < $width ){
            $diagramm .= '*';
            $i++;
        }
        $diagramm .= "\n";
    }
    $diagramm .= '</pre>';

    return $diagramm;
}


1;
