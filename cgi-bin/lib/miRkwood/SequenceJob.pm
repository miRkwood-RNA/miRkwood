package miRkwood::SequenceJob;

# ABSTRACT: Job processing a given sequence

use strict;
use warnings;

use miRkwood::Parsers;
use miRkwood::Programs;
use miRkwood::Utils;

use File::Spec;
use File::Basename;
use Log::Message::Simple qw[msg error debug];

=method new

Constructor

my $sequence_job = miRkwood::SequenceJob->new($sequence_dir, $name, $sequence, '+');

=cut

sub new {
    my ( $class, @args ) = @_;
    my ($directory, $name, $sequence, $strand) = @args;
    my $self = {
        directory => $directory,
        name => $name,
        sequence => $sequence,
        strand => $strand
    };
    bless $self, $class;
    return $self;
}

=method get_strand

Accessor to the strand

=cut

sub get_strand {
    my ($self, @args)  = @_;
    return $self->{'strand'};
}

sub is_opposite_strand {
    my $self = shift;
    return ($self->get_strand() eq '-');
}

sub get_sequence {
    my ($self, @args)  = @_;
    return $self->{'sequence'};
}

=method get_sequence_length

Return the length of the sequence

=cut

sub get_sequence_length{
    my ( $self, @args ) = @_;
    return length $self->get_sequence();
}

=method get_raw_candidates_for_sequence

Usage: my $candidates = $sequence_job->get_raw_candidates_for_sequence();
Returns: A refernce to an array of candidates (as hashes)

=cut

sub get_raw_candidates_for_sequence {
    my ( $self, @args ) = @_;
    my $rnalfold_output = $self->run_rnalfold_on_sequence();

    my ($rnastemloop_out_stemloop, $rnastemloop_out_optimal) =
      $self->run_RNAstemloop_on_rnalfold_output( $rnalfold_output );

    my $rnaeval_out_optimal =
      $self->run_RNAeval_on_RNAstemloop_optimal_output( $rnastemloop_out_optimal );

    my $rnaeval_out_stemloop =
      $self->run_RNAeval_on_RNAstemloop_stemloop_output( $rnastemloop_out_stemloop );

    my $seq_length = length $self->{'sequence'};
    my $candidates = $self->process_RNAstemloop_on_filenames(
                                                $rnastemloop_out_stemloop,
                                                $rnaeval_out_optimal,
                                                $rnaeval_out_stemloop );
    return $candidates;
}

=method run_rnalfold_on_sequence

 Usage : run_rnalfold_on_sequence($sequence);
 Return: $rnalfold_output

=cut

sub run_rnalfold_on_sequence {
    my ( $self, @args ) = @_;
    debug( 'Running RNALfold', miRkwood->DEBUG() );
    my $rnalfold_output = File::Spec->catfile( $self->{'directory'}, 'RNALfold.out' );
    my $temp_file = File::Spec->catfile( $self->{'directory'}, 'tempFile.txt' );
    miRkwood::Programs::run_rnalfold( $self->{'name'}, $self->{'sequence'}, $temp_file,
        $rnalfold_output )
      or die("Problem when running RNALfold: $!");
    return $rnalfold_output;
}

=method run_RNAstemloop_on_rnalfold_output

 Usage : run_RNAstemloop_on_rnalfold_output( $rnalfold_output, $sequence_dir );
 Return: ($rnastemloop_out_stemloop, $rnastemloop_out_optimal);

=cut

sub run_RNAstemloop_on_rnalfold_output {
    my ( $self, @args ) = @_;
    my $rnalfold_output = shift @args;

    my $rnastemloop_out_optimal = File::Spec->catfile( $self->{'directory'}, 'rnastemloop_optimal.out' );
    my $rnastemloop_out_stemloop =
      File::Spec->catfile( $self->{'directory'}, 'rnastemloop_stemloop.out' );
    debug( "Running RNAstemloop on $rnalfold_output", miRkwood->DEBUG() );
    miRkwood::Programs::run_rnastemloop( $rnalfold_output,
        $rnastemloop_out_stemloop, $rnastemloop_out_optimal )
      or die("Problem when running RNAstemloop");
    return ($rnastemloop_out_stemloop, $rnastemloop_out_optimal);
}

=method run_RNAeval_on_RNAstemloop_optimal_output

=cut

sub run_RNAeval_on_RNAstemloop_optimal_output {
    my ($self, @args)  = @_;
    my $rnastemloop_out_optimal = shift @args;
    return $self->run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_optimal, 'optimal' );
}

=method run_RNAeval_on_RNAstemloop_stemloop_output

=cut

sub run_RNAeval_on_RNAstemloop_stemloop_output {
    my ($self, @args)  = @_;
    my $rnastemloop_out_stemloop = shift @args;
    return $self->run_RNAeval_on_RNAstemloop_output( $rnastemloop_out_stemloop, 'stemloop' );
}

=method run_RNAeval_on_RNAstemloop_output


=cut

sub run_RNAeval_on_RNAstemloop_output {
    my ( $self, @args ) = @_;
    my ( $rnastemloop_out, $suffix ) = @args    ;
    my $current_sequence_dir = dirname($rnastemloop_out);
    debug( "Processing RNAstemloop output for $suffix $rnastemloop_out", miRkwood->DEBUG() );
    my $rnaeval_out =
      File::Spec->catfile( $current_sequence_dir, "rnaeval_$suffix.out" );

    debug( "Running RNAeval in $rnaeval_out", miRkwood->DEBUG() );
    miRkwood::Programs::run_rnaeval( $rnastemloop_out, $rnaeval_out )
      or die("Problem when running RNAeval");

    return $rnaeval_out;
}

=method process_RNAstemloop_on_filenames

Pass-through method for process_RNAstemloop

=cut

sub process_RNAstemloop_on_filenames {
    my ($self, @args)  = @_;
    my ($rnastemloop_out_stemloop)      = shift @args;
    my ($rnaeval_out_optimal)  = shift @args;
    my ($rnaeval_out_stemloop) = shift @args;

    open( my $STEM_FH, '<', $rnastemloop_out_stemloop ) or die "Error opening $rnastemloop_out_stemloop: $!";
    open( my $EVAL_OPT_FH, '<', $rnaeval_out_optimal ) or die $!;
    open( my $EVAL_STEM_FH, '<', $rnaeval_out_stemloop ) or die $!;
    my $msg = "Processing RNAstemloop ( $rnastemloop_out_stemloop, $rnaeval_out_optimal, $rnaeval_out_stemloop )";
    debug( $msg, miRkwood->DEBUG() );
    my $candidates =
      $self->process_RNAstemloop($STEM_FH,
        $EVAL_OPT_FH, $EVAL_STEM_FH );
    close($STEM_FH);
    close($EVAL_OPT_FH);
    close($EVAL_STEM_FH);
    return $candidates;
}

=method process_RNAstemloop

Process the outputs of RNAstemloop + RNAeval (both optimal and stemloop)
Returns a list of candidates.

=cut

sub process_RNAstemloop {
    my ($self, @args)  = @_;
    my ($STEM_FH)      = shift @args;
    my ($EVAL_OPT_FH)  = shift @args;
    my ($EVAL_STEM_FH) = shift @args;

    my ($line_eval_opt, $line_eval_stem);
    my ( $nameSeq, $dna, $structure_stemloop );
    my @candidates_array = ();

    while ( my $stem_line = <$STEM_FH> ) {

        if ( miRkwood::Utils::is_fasta_header( $stem_line )) {
            $nameSeq = substr ($stem_line, 1, -1);
        }
        elsif ( miRkwood::Utils::is_fasta_line_relaxed($stem_line ) ) {
            $dna = substr $stem_line, 0, -1;
            $line_eval_opt = substr( <$EVAL_OPT_FH>, 0, -1 );    # the sequence as well
            if ( miRkwood::Utils::is_fasta_header( $line_eval_opt ) ) {
                $line_eval_opt = substr( <$EVAL_OPT_FH>, 0, -1 );
            }
            $line_eval_stem = substr( <$EVAL_STEM_FH>, 0, -1 );    # the sequence as well
             if ( miRkwood::Utils::is_fasta_header( $line_eval_stem )) {
                $line_eval_stem = substr( <$EVAL_STEM_FH>, 0, -1 );
            }
            if ( $dna ne $line_eval_opt || $dna ne $line_eval_stem ) {
                warn ('The sequences differ in RNAeval and RNAstemloop output');
            }
        }
        elsif ( ( $stem_line =~ /(.*)/ ) ) {
            $structure_stemloop = $1;
            $line_eval_opt = <$EVAL_OPT_FH>;    # the structure as well, and the energy
            $line_eval_stem = <$EVAL_STEM_FH>;

            my ( $structure_optimal, $energy_optimal ) =
                miRkwood::Parsers::parse_Vienna_line($line_eval_opt);
            my ( $structure_stemloop, $energy_stemloop ) =
                miRkwood::Parsers::parse_Vienna_line($line_eval_stem);
            if ($structure_optimal)
            {                        # We have a structure

                if ( $nameSeq =~ /.*__(\d*)-(\d*)$/ ) {
                    my ($mfei, $amfe) =
                      miRkwood::Utils::compute_mfei_and_amfe( $dna, $energy_optimal );
                    my ( $start, $end );
                    if ( $self->is_opposite_strand() ) {
                        my $res =
                          miRkwood::Utils::get_position_from_opposite_strand
                          ( $1, $2, $self->get_sequence_length() );
                        ( $start, $end ) = @{$res};
                    }
                    else {
                        ( $start, $end ) = ( $1, $2 );
                    }
                    my $res = {
                        "name"      => $nameSeq,
                        "start"     => $start,
                        "end"       => $end,
                        "mfei"      => $mfei,
                        "amfe"      => $amfe,
                        "dna"       => $dna,
                        "structure_optimal" => $structure_optimal,
                        "structure_stemloop" => $structure_stemloop,
                        "strand"    => $self->get_strand(),
                        "energy_optimal" => $energy_optimal,
                        "energy_stemloop" => $energy_stemloop
                    };
                    push @candidates_array, $res;

                }
            }
            else {
                warn( "No structure found in $line_eval_opt" );
            }    # if $line2
        }
        else {
            warn( "Unrecognized line " );
        }    #if $line1
    }    #while $line=<IN>
    return \@candidates_array;
}

1;