package miRkwood::MultiAlignments_Isa;

# ABSTRACT: Module to create multiple alignments from 2 by 2 alignments

use strict;
use warnings;
use miRkwood;
use miRkwood::Utils;
use Data::Dumper;

sub fillTabTemp2D{
    my (@args) = @_;
    my $alignments = shift @args;
    my $id_candidate = shift @args;
    my $candidate_seq = shift @args;
    my @keys = sort {
        ( miRkwood::Utils::get_element_of_split( $a, '-', 0 )
              <=> miRkwood::Utils::get_element_of_split( $b, '-', 0 ) )
          || ( miRkwood::Utils::get_element_of_split( $a, '-', 1 )
            <=> miRkwood::Utils::get_element_of_split( $b, '-', 1 ) )
    } keys %{$alignments};
    my $additional_data = {};

    foreach my $key (@keys) {
        my @tab = $alignments->{$key};
        my $start_aln = miRkwood::Utils::get_element_of_split( $key, '-', 0 );
        my $end_aln = miRkwood::Utils::get_element_of_split( $key, '-', 1 );
        my $length = $end_aln - $start_aln + 1;
        my $global_query = substr( $candidate_seq, $start_aln - 1, $length);
        $global_query = uc( $global_query );
        my $final_query = [];
        my $final_aln = [];

        $additional_data->{ 'start_aln' } = $start_aln;
        $additional_data->{ 'end_aln' } = $end_aln;

        my $query = miRkwood::Utils::get_element_of_split( $tab[0][0]{'alignment'}, '\n', 0 );
        my $target = miRkwood::Utils::get_element_of_split( $tab[0][0]{'alignment'}, '\n', 2 );
        if ( $tab[0][0]{'end_target'} < $end_aln ){
            for (my $i = $tab[0][0]{'end_target'}; $i < $end_aln; $i++){
                $query = $query . '#';
                $target = $target . '#';
            }
        }

        $additional_data->{ 'def_query' }[0] = $tab[0][0]{'def_query'};

        if ( @{$tab[0]} == 1 ){
            # only one alignment
            my @final_query = split( //, $query );
            $final_query = \@final_query;
            my @target = split( //, $target );
            $final_aln->[0] = \@target;
        }
        else {
            # multiple alignments
            my $positionTab = [];
            my $targets = [];
            $positionTab = init_positionTab( $query, $target );
            $target =~ s/-//g;
            push @{ $targets->[0] }, split( //, $target );

            #~ print STDERR Dumper( $positionTab );
            #~ write_tab ($positionTab);

            for ( my $i = 1; $i < @{$tab[0]}; $i++){
                my $query = miRkwood::Utils::get_element_of_split( $tab[0][$i]{'alignment'}, '\n', 0 );
                my $target = miRkwood::Utils::get_element_of_split( $tab[0][$i]{'alignment'}, '\n', 2 );
                $additional_data->{ 'def_query' }[$i] = $tab[0][$i]{'def_query'};
                if ( $tab[0][$i]{'begin_target'} > $start_aln ){
                    for (my $j = $start_aln; $j < $tab[0][$i]{'begin_target'}; $j++){
                        $query = '#' . $query;
                        $target = '#' . $target;
                    }
                }
                if ( $tab[0][$i]{'end_target'} < $end_aln ){
                    for (my $j = $tab[0][$i]{'end_target'}; $j < $end_aln; $j++){
                        $query = $query . '#';
                        $target = $target . '#';
                    }
                }

                $positionTab = modify_positionTab( $positionTab, $query, $target );
                #~ write_tab ($positionTab);

                $target =~ s/-//g;
                push @{ $targets->[$i] }, split( //, $target );
            }
            ($final_query, $final_aln) = construct_multiple_aln ( $positionTab, $global_query, $targets );
            #~ print STDERR Dumper( $final_aln );
        }
        $additional_data->{ 'final_query' } = $final_query;
        my $stars_line = construct_stars_line( $final_aln, $final_query );
        write_in_file( $id_candidate, $final_aln, $stars_line, $additional_data );
    }
    return;
}

sub init_positionTab {
    my (@args) = @_;
    my $query = shift @args;
    my $target = shift @args;
    my $position_tab = [];
    my @query = split ( //, $query );
    my @target = split ( //, $target );

    # First position on the position tab is given by the number
    # of gaps in query before the query really starts ( + 1 )
    my $i = 0;
    my $nb_left_gaps = 0;
    while ( $query[$i] eq '-' ){
        $nb_left_gaps++;
        $i++;
    }
    for (my $j = 0; $j < $nb_left_gaps; $j++){
        shift @target;
        shift @query;
    }

    $position_tab->[0][0] = $nb_left_gaps + 1;
    my $last_non_null_value = $position_tab->[0][0];

    for (my $i = 1; $i < scalar(@query); $i++ ){
        if ( $target[$i] eq '-' ){
            push @{$position_tab->[0]}, 0;
        }
        else {
            if ( $query[$i] eq '-' ){
                $last_non_null_value++;
            }
            else {
                $last_non_null_value++;
                push @{$position_tab->[0]}, $last_non_null_value;
            }
        }
    }
    return $position_tab;
}

sub construct_multiple_aln {
    my (@args) = @_;
    my $position_tab = shift @args;
    my $query = shift @args;
    my $targets = shift @args;
    #~ print STDERR "construct_multiple_aln is called with :\n";
    #~ print STDERR "position_tab : \n";
    #~ print STDERR Dumper($position_tab);
    #~ print STDERR "query : $query\n";
    #~ print STDERR "targets : \n";
    #~ print STDERR Dumper($targets);
    my $final_aln = [];
    $query =~ s/-//g;
    my @query = split( //, $query );
    my $final_query;
    my $nb_gaps;
    my $diff_max = 0;
    my $nb_aln = scalar( @{$position_tab} );

    ### initialise the final query and targets.
    ### deal with the fact that some targets may start earlier than the query.
    my $left_gaps = 0;
    # count how much gaps should be added on the left of the alignments
    for (my $aln = 0; $aln < scalar( @{$position_tab} ); $aln++){
        if ( $position_tab->[$aln][0] - 1 > $left_gaps) {
            $left_gaps = $position_tab->[$aln][0] - 1;
        }
    }
    # add gaps for the query
    for (my $i = 0; $i < $left_gaps; $i++){
        push @{ $final_query }, '-';
    }
    push @{ $final_query }, $query[0];
    # add the appropriate number of gaps for the targets
    for (my $aln = 0; $aln < scalar( @{$position_tab} ); $aln++){
        for (my $i = 0; $i < $left_gaps - $position_tab->[$aln][0] + 1; $i++){
            push @{ $final_aln->[$aln] }, '-';
        }
        for (my $i = 0; $i < $position_tab->[$aln][0]; $i++){
            push @{ $final_aln->[$aln] }, shift @{ $targets->[$aln] };
        }
    }

    # main part of the alignment
    for (my $pos = 1; $pos < scalar( @query ); $pos++){
        $diff_max = 0;
        for (my $aln = 0; $aln < $nb_aln; $aln++){
            my $diff = $position_tab->[$aln][$pos] - $position_tab->[$aln][$pos-1];
            if ( $position_tab->[$aln][$pos-1] == 0 ){
                $diff = 1;
            }
            $nb_gaps->[$aln] = $diff;
            if ( $diff > $diff_max ){
                $diff_max = $diff;
            }
        }

        if ( $diff_max == 1 ){
            # no indel on this position
            for (my $aln = 0; $aln < $nb_aln; $aln++){
                if ( $nb_gaps->[$aln] < 0 ){
                    # gap in target
                    push @{ $final_aln->[$aln] }, '-';
                }
                else {
                    # match / mismatch
                    push @{ $final_aln->[$aln] }, shift @{ $targets->[$aln] };
                }
            }
        }
        else {  # one or several targets have an insertion
            for (my $j = 0; $j < $diff_max - 1; $j++){
                push @{ $final_query }, '-';
            }
            for (my $aln = 0; $aln < $nb_aln; $aln++){
                # add the needed nucleotides or gaps
                for (my $j = 0; $j < $diff_max - 1; $j++){
                    if ( $nb_gaps->[$aln] == $diff_max ){
                        push @{ $final_aln->[$aln] }, shift @{ $targets->[$aln] };
                    }
                    else {
                        push @{ $final_aln->[$aln] }, '-';
                    }
                }
                # add the current letter
                if ( $nb_gaps->[$aln] < 0 ){
                    push @{ $final_aln->[$aln] }, '-';
                }
                else {
                    push @{ $final_aln->[$aln] }, shift @{ $targets->[$aln] };
                }
            }
        }
        push @{ $final_query }, $query[$pos];
    }
    ### deal with possible gaps at the right end of the query
    # look for the longest target
    my $length_max = 0;
    for (my $aln = 0; $aln < $nb_aln; $aln++){
        if ( $targets->[$aln] ){
            push @{ $final_aln->[$aln] }, @{$targets->[$aln]};
        }
        if ( scalar( @{ $final_aln->[$aln] } ) > $length_max ){
            $length_max = scalar( @{ $final_aln->[$aln] } );
        }
    }
    # add gaps at the right end of the query
    while ( scalar( @{ $final_query } ) < $length_max ){
        push @{ $final_query }, '-';
    }
    # add gaps at the right end of targets
    for (my $aln = 0; $aln < $nb_aln; $aln++){
        while ( scalar( @{ $final_aln->[$aln] } ) < $length_max ){
            push @{ $final_aln->[$aln] }, '-';
        }
    }
    return ($final_query, $final_aln);
}


sub modify_positionTab {
    my (@args) = @_;
    my $position_tab = shift @args;
    my $query = shift @args;
    my $target = shift @args;
    my $query_without_gaps = $query;
    $query_without_gaps =~ s/-//g;
    my @query = split ( //, $query );
    my @target = split ( //, $target );
    my $curr_aln = scalar( @{$position_tab} );

    my $position_target = 1;
    my $index_seq = 0;
    for (my $i = 0; $i < length( $query_without_gaps ); $i++ ){
        if ( $target[ $index_seq ] eq '-' ){
            $position_tab->[$curr_aln][$i] = 0;
            #~ print STDERR "$query[$i] is paired with a gap\n";
        }
        else{
            if ( $query[ $index_seq ] eq '-' ){
                $index_seq++;
                #~ print STDERR "There is an indel around position $i\n";
            }
            $position_tab->[$curr_aln][$i] = ($index_seq + 1);
            #~ print STDERR "$query[$i] is paired with the ".($index_seq + 1). " position\n";
        }
        $index_seq++;
    }

    return $position_tab;
}


sub write_tab {
    my @args = @_;
    my $position_tab = shift @args;
    for (my $i = 0; $i < scalar( @{$position_tab} ); $i++){
        print STDERR "aln $i | ". join( ' ', @{$position_tab->[$i]})."\n";
    }
    return;
}

sub construct_stars_line {
    my @args = @_;
    my $aln_tab = shift @args;
    my $query = shift @args;
    my $length = scalar( @{$query} );
    my $stars_line = [];
    for (my $pos = 0; $pos < $length; $pos++ ){
        my $char = '*';
        for (my $aln = 0; $aln < scalar( @{$aln_tab} ); $aln++ ){
            if ( $aln_tab->[$aln][$pos] ne $query->[$pos] ){
                $char = ' ';
            }
        }
        push @{$stars_line}, $char;
    }
    return $stars_line;
}

sub write_in_file {
    my @args = @_;
    my $id_candidate = shift @args;
    my $aln_tab = shift @args;
    my $stars_line = shift @args;
    my $additional_data = shift @args;
    my $cfg = miRkwood->CONFIG();
    my $aln_dir = miRkwood::Paths::get_dir_alignments_path_from_job_dir( $cfg->param('job.directory') );
    my $output = '';
    my $output_tab = {};

    $output_tab->{ 'query' } = '    query       '.$additional_data->{'start_aln' }.'  ';
    # all alignments must start at the same position
    for (my $aln = 0; $aln < scalar( @{$aln_tab} ); $aln++) {
        $output_tab->{ 'miRBase ' . ($aln + 1) } = '    miRBase ' . ($aln + 1);
        while ( length( $output_tab->{ 'miRBase ' . ($aln + 1) } ) < length( $output_tab->{ 'query' } ) ){
            $output_tab->{ 'miRBase ' . ($aln + 1) } .= ' ';
        }
    }
    $output_tab->{ 'stars' } = '';
    while ( length( $output_tab->{ 'stars' } ) < length( $output_tab->{ 'query' } ) ){
        $output_tab->{ 'stars' } .= ' ';
    }
    $output_tab->{ 'stars' } .= join( '', @{$stars_line} );

    # add the alignments
    $output_tab->{ 'query' } .= join( '', @{$additional_data->{ 'final_query' }} ) . '  ' . $additional_data->{'end_aln'};
    for (my $aln = 0; $aln < scalar( @{$aln_tab} ); $aln++) {
        $output_tab->{ 'miRBase ' . ($aln + 1) } .= join( '', @{$aln_tab->[$aln]} );
        $output_tab->{ 'miRBase ' . ($aln + 1) } =~ s/#/-/g;
    }

    # write on output
    $output .= "Prediction : $additional_data->{'start_aln' }-$additional_data->{'end_aln' }\n\n";
    $output .= $output_tab->{ 'query' }."\n";
    for (my $aln = 0; $aln < scalar( @{$aln_tab} ); $aln++) {
        $output .= $output_tab->{ 'miRBase ' . ($aln + 1) }."\n";
    }
    $output .= $output_tab->{ 'stars' }."\n\n";
    for (my $aln = 0; $aln < scalar( @{$aln_tab} ); $aln++) {
        $output .= 'miRBase ' . ($aln + 1) . ": $additional_data->{'def_query'}[$aln]\n";
    }
    $output .= "\n";

    #~ print STDERR $output;
    
    open(my $FILE, '>>', "$aln_dir/${id_candidate}_aln.txt") or die "ERROR when creating $aln_dir/${id_candidate}_aln.txt: $!";
    print $FILE $output;
    close($FILE);

    return;

}


1;
