package PipelineMiRNA::WebFunctions;

use strict;
use warnings;
use Data::Dumper;

sub resultstruct2csv{
    my ($self, @args) = @_;
    my $size = shift @args;
    my $results = shift @args;
    
    my @headers= ('Name',
                  'Positions',
                  'MFEI',
                  'MFE',
                  'AMFE',
                  'Pvalue',
                  'sequence',
                  'secondary structure',
                  );
    my $res =  join( ',', @headers ) . "\n";
    for (my $i=0; $i<$size; $i++)
    {
        my @array = ($results->names($i),
                     $results->positions($i),
                     $results->mfeis($i),
                     $results->mfes($i),
                     $results->amfes($i),
                     $results->pvalues($i),
                     $results->DNASequence($i),
                     $results->Vienna($i),
                     $results->selfContain($i),
                     );
        $res .= join( ',', @array ) . "\n";
    }
    return $res;
}

1;
