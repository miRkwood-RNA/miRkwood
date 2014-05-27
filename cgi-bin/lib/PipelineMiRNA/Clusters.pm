package PipelineMiRNA::Clusters;

# ABSTRACT: Handling the cluster making from a BAM file.

use strict;
use warnings;

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
        $fai_hash{ $fields[0] } = $fields[1];
        $genome_sum += $fields[1];
    }
    close $FAI;

    my @islands = ();

    my $last_start;
    my $last_ok;

    my $island;

    foreach my $chr (@chrs) {

        $last_start = -1;
        $last_ok    = -1;

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

                    $island = "$chr" . ':' . "$last_start" . '-' . "$last_ok";
                    push( @islands, $island );
                    $last_start = $fields[1];
                    $last_ok    = $fields[1];
                }
            }
        }
        close $DEPTH;
        ## close the last one, if there is one
        if ( $last_ok != -1 ) {
            $island = "$chr" . ':' . "$last_start" . '-' . "$last_ok";
            push( @islands, $island );
        }

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
    if (! -e $fai_file ) {
        warn(
"\nFatal in sub-routine get_folding_regions : expected fai file $fai_file does not exist\n"
        );
        exit;
    }
    my @fai_fields = ();
    open( my $FAI,  '<', "$fai_file" )
        or die "Error when opening $fai_file: $!";
    while (<$FAI>) {
        @fai_fields = split( "\t", $_ );
        $chr_sizes{ $fai_fields[0] } = $fai_fields[1];
    }
    close $FAI
        or die "Error when closing $fai_file: $!";
    my $this_chr;
    my $entry;
    foreach my $in_clus (@{$input}  ) {
        if ( $in_clus =~ /^(\S+):(\d+)-(\d+)$/ ) {
            $this_chr          = $1;
            $this_start        = $2;
            $this_padded_start = $this_start - $pad;
            if ( $this_padded_start < 1 ) {
                $this_padded_start = 1;
            }
            $this_stop        = $3;
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
                $entry = "$last_chr" . ':' . "$last_start" . '-' . "$last_stop";
                push( @output, $entry );
                $last_padded_start = $this_padded_start;
                $last_padded_stop  = $this_padded_stop;
                $last_start        = $this_start;
                $last_stop         = $this_stop;
            }
            else {
                if ( $this_padded_start > $last_padded_stop ) {
                    ## no overlap between these padded clusters.
                    ## Report the last one, trimming off its dangling pads
                    $entry =
                      "$this_chr" . ':' . "$last_start" . '-' . "$last_stop";
                    push( @output, $entry );
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
        }
        else {
            warn(
"\nFATAL: in sub-routine \'merge_clusters\' : failed to parse initial locus $in_clus\n"
            );
            exit;
        }
    }
    ## last one
    $entry = "$last_chr" . ':' . "$last_start" . '-' . "$last_stop";
    push( @output, $entry );
    return @output;
}

1;
