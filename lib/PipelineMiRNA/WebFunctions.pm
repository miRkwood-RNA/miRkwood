package PipelineMiRNA::WebFunctions;

use strict;
use warnings;
use Data::Dumper;
use File::Spec;
use PipelineMiRNA::Paths;
use PipelineMiRNA::Parsers;

my @headers = ('name', 'position', 'mfei', 'mfe', 'amfe', 'p_value', 'self_contain', 'alignment', 'image', 'Vienna', 'DNASequence');

=method jobId_to_jobPath

Get the job path from a job identifier

=cut

sub jobId_to_jobPath {
    my ( $self, @args ) = @_;
    my $id_job      = shift @args;
    my $dirJob_name = 'job' . $id_job;
    my $jobPath = File::Spec->catdir( 'results', $dirJob_name);
    return $jobPath;
}

=method is_valid_jobID

Test whether a jobID is valid - ie if there are results for it.

=cut

sub is_valid_jobID {
    my ( $self, @args ) = @_;
    my $id_job          = shift @args;
    my $jobPath = $self->jobId_to_jobPath($id_job);
    my $full_path = PipelineMiRNA::Paths->get_absolute_path($jobPath);
    return (-e $full_path);
}

=method get_structure_for_jobID

Get the results structure of a given job identifier

=cut

sub get_structure_for_jobID {
    my ( $self, @args ) = @_;
    my $jobId   = shift @args;
    my $job = $self->jobId_to_jobPath($jobId);
    my $job_dir = PipelineMiRNA::Paths->get_absolute_path($job);
    my %myResults = ();

    opendir DIR, $job_dir;    #ouverture répertoire job
    my @dirs;
    @dirs = readdir DIR;
    closedir DIR;
    foreach my $dir (@dirs)    # parcours du contenu
    {
        my $full_dir = File::Spec->catdir( $job_dir, $dir );
        if (    $dir ne "."
             && $dir ne ".."
             && -d $full_dir )    #si fichier est un répertoire
        {
            opendir DIR, $full_dir;    # ouverture du sous répertoire
            my @files;
            @files = readdir DIR;
            closedir DIR;
            foreach my $subDir (@files) {
                my $subDir_full = File::Spec->catdir( $job_dir, $dir, $subDir );
                if (    ( $subDir ne "." )
                     && ( $subDir ne ".." )
                     && -d $subDir_full ) # si le fichier est de type repertoire
                {
                    $myResults{$subDir} = ();

                    $myResults{$subDir}{'name'} = $dir;    #récupération nom séquence
                    my @position = split( /__/, $subDir );
                    $myResults{$subDir}{'position'} = $position[1]; # récupération position
                      ;    # récupération noms + positions
                    my $pvalue =
                      File::Spec->catfile( $subDir_full, 'pvalue.txt' );
                    if ( -e $pvalue )    # si fichier existe
                    {
                        $myResults{$subDir}{'p_value'} = PipelineMiRNA::Parsers::parse_pvalue($pvalue);
                    }

                    #Récupération valeur MFEI
                    my $mfei_out =
                      File::Spec->catfile( $subDir_full, 'outMFEI.txt' );
                    if ( -e $mfei_out )                 # si fichier existe
                    {
                        my @mfeis = PipelineMiRNA::Parsers::parse_mfei($mfei_out);
                        $myResults{$subDir}{'mfei'} = $mfeis[0];
                        $myResults{$subDir}{'mfe'} = $mfeis[1];
                        $myResults{$subDir}{'amfe'} = $mfeis[2];
                    }

                    #Récupération valeur self contain
                    my $selfcontain_out =
                      File::Spec->catfile( $subDir_full, 'selfContain.txt' );
                    if ( -e $selfcontain_out )
                    {
                        $myResults{$subDir}{'self_contain'} = PipelineMiRNA::Parsers::parse_selfcontain($selfcontain_out);
                    }

                    #Récupération séquence et format Vienna
                    my $vienna_out = File::Spec->catfile( $subDir_full,
                                                       'outViennaTraited.txt' );
                    if ( -e $vienna_out )                  # si fichier existe
                    {
                        my @vienna_res = PipelineMiRNA::Parsers::parse_vienna($vienna_out);
                        $myResults{$subDir}{'DNASequence'} = $vienna_res[0];
                        $myResults{$subDir}{'Vienna'} = $vienna_res[1];
                    }

                    #Récupération alignement avec mirBase
                    my $file_alignement =
                      File::Spec->catfile( $subDir_full, 'alignement.txt' );
                    if ( -e $file_alignement )                # si fichier existe
                    {
                        $myResults{$subDir}{'alignment'} = PipelineMiRNA::Parsers::parse_alignment($file_alignement);
                    }
                    my $image_path = PipelineMiRNA::Paths->get_server_path($job,  $dir, $subDir, 'image.png');
                    $myResults{$subDir}{'image'} = $image_path;
                }
            }
        }
    }

    return %myResults;
}

=method resultstruct2csv

Convert the results structure to CSV

=cut

sub resultstruct2csv {
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my %results = %{$results};
    my @headers = ('name', 'position', 'mfei', 'mfe', 'amfe', 'p_value', 'self_contain', 'Vienna', 'DNASequence');

    my $result = join( ',', @headers ) . "\n";
    while ( my ($key, $value) = each %results )
    {
        for my $header (@headers){
            $result .= "${$value}{$header},";
        }
        $result .= "\n";
    }
    return $result;
}

=method resultstruct2table

Convert the results structure to HTML table

=cut

sub resultstruct2table {
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my %results = %{$results};

    my $HTML_results = <<'END_TXT';
            <div class="titreDiv"> Identification of miRNA/miRNA hairpins results:</div>
            <div id="table" ></div>
END_TXT

    my $row = 0;
    my $column = -1;
    $HTML_results .= "<table>\n<tbody>";
    $HTML_results .= "<tr>";
    for my $header (@headers){
        $column += 1;
        my $th_content = "id='cell-$row-$column' width='100' onclick='showCellInfo($row, $column)' onmouseover='colorOver($row, $column)' onmouseout='colorOut($row, $column)'";
        $HTML_results .= "<th $th_content>$header</th>\n";
    }
    $HTML_results .= "</tr>\n";
    while ( my ($key, $value) = each %results )
    {
        $row += 1;
        $column = -1;
      $HTML_results .= '<tr>';
      for my $header (@headers){
          $column += 1;
          my $td_content = "id='cell-$row-$column' onmouseover='colorOver($row, $column)' onmouseout='colorOut($row, $column)' onclick='showCellInfo($row, $column)'";
          $HTML_results .= "<th $td_content>${$value}{$header}</th>\n";
      }
      $HTML_results .= "\n</tr>\n";
    }
    $HTML_results .= "</tbody></table>";
    return $HTML_results;
}

=method resultstruct2pseudoXML

Convert the results structure to to pseudo XML format

=cut

sub resultstruct2pseudoXML {
    my ( $self, @args ) = @_;
    my $results = shift @args;
    my %results = %{$results};
    my $result = "<results id='all'>\n";
    while ( my ($key, $value) = each %results )
    {
        $result .= "<Sequence";
        for my $header (@headers){
            $result .= " $header='${$value}{$header}'";
        }
        $result .= "></Sequence>\n";
    }
    $result .= "</results>";
    return $result;
}

1;
