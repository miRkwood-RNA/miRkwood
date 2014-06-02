package PipelineMiRNA::Clusters;

# ABSTRACT: Handling the cluster making from a BAM file.

use strict;
use warnings;

use PipelineMiRNA::Utils;

=method get_clustered_sequences_from_bam

Given a BAM file

=cut

sub get_clustered_sequences_from_bam {
    my ( $self, @args ) = @_;
    my ( $bamfile, $genome, $mindepth, $pad ) = @args;
    my @clusters = $self->get_clusters( $bamfile, $genome, $mindepth, $pad );
    my @sequences = $self->get_sequences_from_clusters( $genome, \@clusters );
    return @sequences;
}

=method get_sequences_from_clusters

Returns the sequences extracted from the given genome
based on the given clusters.

=cut

sub get_sequences_from_clusters {
    my ( $self,   @args )     = @_;
    my ( $genome, $clusters ) = @args;
    my @clusters = @{$clusters};
    my %sequences_hash = PipelineMiRNA::Utils::multifasta_to_hash($genome);
    my @result;
    foreach my $cluster (@clusters) {
        my ( $chr, $start, $stop ) = @{$cluster};
        #$self->extend_cluster($cluster);
        my $original_seq = $sequences_hash{$chr};
        my $sequence     = substr $original_seq, $start, $stop;
        my $new_name     = $chr . "__$start-$stop";
        my @res          = ( $new_name, $sequence );
        push @result, \@res;
    }
    return @result;
}

=method get_clusters

Retrieve the clusters from the given BAM, genome
and clustering options.

Serves essentially as a wrapper around
- get_faidx_file
- get_islands
- merge_clusters

 Usage : my @clusters = get_clusters( $bamfile, $genome_file, $mindepth, $pad)

=cut

sub get_clusters {
    my ( $self, @args ) = @_;
    my ( $bamfile, $genome, $mindepth, $pad ) = @args;
    my $expected_faidx = $self->get_faidx_file($genome);
    my @islands = $self->get_islands( $bamfile, $mindepth, $expected_faidx );
    my @clusters = $self->merge_clusters( \@islands, $pad, $genome );
    return @clusters;
}

=method get_faidx_file


=cut

sub get_faidx_file {
    my ( $self, @args ) = @_;
    my $genome_file    = shift @args;
    my $expected_faidx = $genome_file . '.fai';
    if ( !-e $expected_faidx ) {
        my $samtools_cmd = "samtools faidx $genome_file";
        system $samtools_cmd;
    }
    return $expected_faidx;
}

# The following methods are taken from ShortStack by Michael J. Axtell
# Licensed under the GNU GPL v3

=method get_islands


my @islands = get_islands($bamfile,$mindepth,$expected_faidx,$read_group);

=cut

sub get_islands {
    my ( $self, @args ) = @_;
    my ( $bamfile, $mindepth, $expected_faidx, $read_group ) = @args;

    # go chr by chr, using the .fai index file to get the chr names
    my @chrs       = ();
    my @fields     = ();
    my $genome_sum = 0;
    my %fai_hash   = ();

    open( my $FAI, '<', "$expected_faidx" )
      or die "Error when opening $expected_faidx: $!";
    while (<$FAI>) {
        chomp;
        @fields = split( "\t", $_ );
        push( @chrs, $fields[0] );
    }
    close $FAI;

    my @islands = ();
    foreach my $chr (@chrs) {
        my $samtools_cmd = '';
        if ($read_group) {
            $samtools_cmd =
"samtools view -F 0x4 -r $read_group -b -u $bamfile $chr | samtools depth /dev/stdin 2> /dev/null |";
        }
        else {
            $samtools_cmd =
"samtools view -F 0x4 -b -u $bamfile $chr | samtools depth /dev/stdin 2> /dev/null |";
        }
        open( my $DEPTH, $samtools_cmd );
        my @chr_islands =
          $self->process_samtools_depth( $DEPTH, $chr, $mindepth );
        push @islands, @chr_islands;
        close $DEPTH;
    }
    return @islands;
}

=method process_samtools_depth


=cut

sub process_samtools_depth {
    my ( $self, @args ) = @_;
    my ( $DEPTH, $chr, $mindepth ) = @args;
    my @fields;
    my $last_start = -1;
    my $last_ok    = -1;
    my @islands    = ();
    while (<$DEPTH>) {
        chomp;
        @fields = split( "\t", $_ );
        if ( $fields[2] >= $mindepth ) {
            if ( $last_ok == -1 ) {
                ## first start on the chr
                $last_start = $fields[1];
                $last_ok    = $fields[1];
            }
            elsif ( $fields[1] == $last_ok + 1 ) {
                ## continuing to extend an island
                $last_ok = $fields[1];
            }
            else {
                ## close the previous island, open a new one
                my @island = ( $chr, $last_start, $last_ok );
                push( @islands, \@island );
                $last_start = $fields[1];
                $last_ok    = $fields[1];
            }
        }
    }
    ## close the last one, if there is one
    if ( $last_ok != -1 ) {
        my @island = ( $chr, $last_start, $last_ok );
        push( @islands, \@island );
    }
    return @islands;
}

=method merge_clusters

my @clusters = merge_clusters(\@islands,\$pad,\$genome);

=cut

sub merge_clusters {
    my ( $self, @args ) = @_;
    my ( $input, $pad, $genome ) = @args;
    my @output = ();

    my $this_start;
    my $this_stop;
    my $last_start;
    my $last_stop;
    my $last_padded_start;
    my $last_padded_stop;
    my $this_padded_start;
    my $this_padded_stop;
    my $last_chr  = 'null';
    my %chr_sizes = ();
    ## grab the chrom sizes, which you need to ensure that you don't pad off the end of the chroms
    # chrom sizes in column 1 from the fai file
    my $fai_file = "$genome" . '.fai';
    if ( !-e $fai_file ) {
        warn(
"\nFatal in sub-routine merge_clusters : expected fai file $fai_file does not exist\n"
        );
        exit;
    }
    my @fai_fields = ();
    open( my $FAI, '<', "$fai_file" )
      or die "Error when opening $fai_file: $!";
    while (<$FAI>) {
        @fai_fields = split( "\t", $_ );
        $chr_sizes{ $fai_fields[0] } = $fai_fields[1];
    }
    close $FAI
      or die "Error when closing $fai_file: $!";
    my $this_chr;

    foreach my $in_clus ( @{$input} ) {
        ( $this_chr, $this_start, $this_stop ) = @{$in_clus};

        $this_padded_start = $this_start - $pad;
        if ( $this_padded_start < 1 ) {
            $this_padded_start = 1;
        }

        $this_padded_stop = $this_stop + $pad;
        if ( $this_padded_stop > $chr_sizes{$this_chr} ) {
            $this_padded_stop = $chr_sizes{$this_chr};
        }

        # special first case
        if ( $last_chr eq 'null' ) {
            $last_padded_start = $this_padded_start;
            $last_padded_stop  = $this_padded_stop;
            $last_start        = $this_start;
            $last_stop         = $this_stop;
        }
        elsif ( $this_chr ne $last_chr ) {
            my @entry = ( $last_chr, $last_start, $last_stop );
            push( @output, \@entry );
            $last_padded_start = $this_padded_start;
            $last_padded_stop  = $this_padded_stop;
            $last_start        = $this_start;
            $last_stop         = $this_stop;
        }
        else {
            if ( $this_padded_start > $last_padded_stop ) {
                ## no overlap between these padded clusters.
                ## Report the last one, trimming off its dangling pads
                my @entry = ( $this_chr, $last_start, $last_stop );
                push( @output, \@entry );
                $last_padded_start = $this_padded_start;
                $last_padded_stop  = $this_padded_stop;
                $last_start        = $this_start;
                $last_stop         = $this_stop;
            }
            else {

   # here, same chr, this_padded_start is <= last_padded_stop, so we are merging
                if ( $this_padded_start < $last_padded_start ) {
                    $last_padded_start = $this_padded_start;
                }
                if ( $this_start < $last_start ) {
                    $last_start = $this_start;
                }
                if ( $this_padded_stop > $last_padded_stop ) {
                    $last_padded_stop = $this_padded_stop;
                }
                if ( $this_stop > $last_stop ) {
                    $last_stop = $this_stop;
                }
            }
        }
        $last_chr = $this_chr;

    }    # foreach locus
    ## last one
    my @entry = ( $last_chr, $last_start, $last_stop );
    push( @output, \@entry );
    return @output;
}

=method extend_cluster

Extend the cluster into an appropriate region

=cut

sub extend_cluster {
    my ( $self, @args ) = @_;
    my $cluster = shift @args;
    my ( $chr, $start, $stop ) = @{$cluster};
    my $locus_size = $stop - $start + 1;
    my $get_size = 300;
    # get the coordinates .. centered on the middle of the cluster
    # first get the middle of the cluster
    my $middle = int(($locus_size / 2) + $start);
    my $get_start = int ($middle - ($get_size / 2));
    my $get_stop = int ($middle + ($get_size / 2));
    # ensure that the get_start is not less than 1
    if($get_start < 1) {
        $get_start = 1;
    }
    return ( $chr, $get_start, $get_stop );
}

1;
