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


=method filterBEDfile

Method to filter a BED from the reads corresponding to known 
miRNAS and from other features (given in gff) for a 
model organism.
Filtered reads are written in auxiliary files.
Input : - a BED file
        - a species name (to find the annotations files)

=cut
sub filterBEDfile {
    my @args = @_;
    my $class    = shift @args;
    my $bed_file = shift @args;
    my $cfg      = miRkwood->CONFIG();
    my $species            = $cfg->param('job.plant');
    my $filter_multimapped = $cfg->param('options.filter_multimapped');
    my $annotation_gff     = $cfg->param( 'options.annotation_gff' );
    my @annotation_gff     = split( /\&/, $annotation_gff );
    my $job_dir            = $cfg->param('job.directory');
    my $basename           = $cfg->param('job.bed');
    my $mirbase_file       = $cfg->param( 'options.mirbase_gff' );
    my $multimapped_interval = $cfg->param('options.multimapped_interval');
    my $min_nb_positions = 0;
    my $max_nb_positions = 5;
    if ( $multimapped_interval =~ /\[(.*);(.*)\]/ ){
        $min_nb_positions = $1;
        $max_nb_positions = $2;
    }

    my $mirna_reads       = File::Spec->catfile( $job_dir , "${basename}_miRNAs.bed");
    my $multimapped_reads = File::Spec->catfile( $job_dir , "${basename}_multimapped.bed");
    my $filtered_bed      = File::Spec->catfile( $job_dir , "${basename}_filtered.bed");

    my $output_bed = File::Spec->catfile( $job_dir , "${basename}_tmp_1.bed");

    ### Filter out known miRNAs
    if ( -r $mirbase_file ) {
        debug( 'Filter out known miRNAS' . ' [' . localtime() . ']', miRkwood->DEBUG() );

        # Create a file with known miRNAs
        store_overlapping_reads( $bed_file, $mirbase_file, $mirna_reads, '-f 1');

        # Delete reads corresponding to known miRNAs from the BED
        store_non_overlapping_reads( $bed_file, $mirbase_file, $output_bed);  # un-comment when tests are done
        #~ system("cp $bed_file $output_bed");  # comment when tests are done
    }
    else{
        debug( "WARNING : no miRNA file found for $species.", 1 );
        system("cp $bed_file $output_bed");
    }

    ### Filter out according to given annotation gff (other than known miRNAs)
    my $i = 1;
    if ( @annotation_gff ){
        foreach my $gff ( @annotation_gff ){
            if ( -r $gff ){
                debug( "Filter out according to $gff" . ' [' . localtime() . ']', miRkwood->DEBUG() );
                $i++;

                my $basename_gff = '';
                if ( $gff =~ /[\/\\]([^\/\\]+)[.](gff3?|dat)/ ){
                    $basename_gff .= $1;
                }
                my $input_gff_bed = File::Spec->catfile( $job_dir , "${basename}_tmp_".($i-1).'.bed');
                my $output_gff_bed = File::Spec->catfile( $job_dir , "${basename}_tmp_$i.bed");
                my $discarded_reads = File::Spec->catfile( $job_dir , "${basename}_$basename_gff.bed");

                filter_according_given_gff( $gff, $input_gff_bed, $output_gff_bed, $discarded_reads );
            }
            else {
                debug( "$gff is not a readable file.", miRkwood->DEBUG() );
            }
        }
    }
    else{
        debug( 'No annotation gff file were provided.', miRkwood->DEBUG() );
    }

    ### Filter out multimapped reads
    if ( $max_nb_positions != 0 ){
        debug( "Filter out multimapped reads (only keep reads between $min_nb_positions and $max_nb_positions positions)" . ' [' . localtime() . ']', miRkwood->DEBUG() );
        filter_multimapped_reads( 
             File::Spec->catfile( $job_dir , "${basename}_tmp_$i.bed"),
             $filtered_bed,
             $multimapped_reads,
             $min_nb_positions,
             $max_nb_positions );
    }
    else{
        debug( 'Don\'t filter out multimapped reads from BED.' . ' [' . localtime() . ']', miRkwood->DEBUG() );
        rename File::Spec->catfile( $job_dir , "${basename}_tmp_$i.bed"), $filtered_bed;
    }

    for (my $j = 1; $j <= $i; $j++){
        unlink File::Spec->catfile( $job_dir , "${basename}_tmp_$j.bed");
    }

    return ($filtered_bed, $mirna_reads);

}


sub filter_according_given_gff {
    my @args = @_;
    my $annotation_gff  = shift @args;
    my $input_bed       = shift @args;
    my $output_bed      = shift @args;
    my $discarded_reads = shift @args;

    # Create a file with reads corresponding to the features in the gff file
    store_overlapping_reads( $input_bed, $annotation_gff, $discarded_reads, '');

    # Delete reads corresponding to the features in the gff file
    store_non_overlapping_reads( $input_bed, $annotation_gff, $output_bed);

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
    #~ debug( "   $job", 1);
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
    #~ debug( "   $job", 1);
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
    my $min_nb_positions = shift @args;
    my $max_nb_positions = shift @args;      

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
        if ( $counts->{$line[3]} >= $min_nb_positions && $counts->{$line[3]} <= $max_nb_positions ){
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
    my $nb_total_reads = 0;
    my $nb_read_unq = 0;

    if ( -e $bed_file ){
        open(my $BED, '<', $bed_file) or die "ERROR while opening $bed_file : $!";
        while ( <$BED> ){
            @tab = split ( /\t/ );
            if ( ( $start == -1 && $end == -1 ) || ( $tab[1] >= $start  && $tab[2] <= $end ) ) {
                $nb_total_reads += $tab[4];
            }
        }
        close $BED;
        $nb_read_unq = `wc -l $bed_file`;
        if ( $nb_read_unq =~ /(\d+) .*/ ){
            $nb_read_unq = $1;
        }

        return ($nb_total_reads, $nb_read_unq);
    }
    else {
        return (0, 0);
    }

}

=method store_reads_nb_in_BED_file

  Count the number of reads and unique
  reads in a BED file and write it in
  the given txt file.

=cut
sub store_reads_nb_in_BED_file {
    my @args = @_;
    my $BED_file = shift @args;
    my $log_file = shift @args;

    my $basename_bed = '';
    if ( $BED_file =~ /.*[\/\\]([^\/\\]+)/ ){
        $basename_bed = $1;
    }
    my ($nb_reads, $nb_unique_reads) = count_reads_in_bed_file( $BED_file, -1, -1 );
    #~ $nb_reads = miRkwood::Utils::make_numbers_more_readable( $nb_reads );
    #~ $nb_unique_reads = miRkwood::Utils::make_numbers_more_readable( $nb_unique_reads );
    open (my $FH, '>>', $log_file) or die "ERROR while opening $log_file : $!";
    print $FH "$basename_bed\t$nb_reads\t$nb_unique_reads\n";
    close $FH;
    return;
}

=method zipBEDfile

=cut
sub zipBEDfile{
    my @args = @_;
    my $BED_file = shift @args;
    my $path = '';
    my $basename = '';
    if ( -e $BED_file ){
        if ( $BED_file =~ /(.*)\/([^\/]+)\.bed/ ){
            $path = $1;
            $basename = $2;
        }
        system("tar zcf $path/$basename.tar.gz -C $path $basename.bed");
        unlink "$BED_file";
    }
    return;
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
