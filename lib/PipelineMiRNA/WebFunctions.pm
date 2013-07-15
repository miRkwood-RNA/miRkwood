package PipelineMiRNA::WebFunctions;

use strict;
use warnings;
use Data::Dumper;
use File::Spec;
use PipelineMiRNA::Paths;

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
                        open( PVALUE, $pvalue ) || die "$!";
                        while ( my $line = <PVALUE> ) {
                            if ( $line =~ /(.*)\t(.*)\t(.*)/ ) {
                                $myResults{$subDir}{'p_value'} = $3;  # récupération pvalues
                            }
                        }
                        close PVALUE;
                    }

                    #Récupération valeur MFEI
                    my $mfei_out =
                      File::Spec->catfile( $subDir_full, 'outMFEI.txt' );
                    if ( -e $mfei_out )                 # si fichier existe
                    {
                        open( MFEI, $mfei_out ) || die "$!";
                        while ( my $line = <MFEI> ) {
                            if ( $line =~ /(.*)\t(.*)\t(.*)\t(.*)/ ) {
                                $myResults{$subDir}{'mfei'} = $2;    # récupération mfei
                                $myResults{$subDir}{'mfe'} = $3;    # récupération mfei
                                $myResults{$subDir}{'amfe'} = $4;    # récupération mfei
                            }
                        }
                        close MFEI;
                    }

                    #Récupération valeur self contain
                    my $selfcontain_out =
                      File::Spec->catfile( $subDir_full, 'selfContain.txt' );
                    if ( -e $selfcontain_out )          # si fichier existe
                    {
                        open( SC, $selfcontain_out ) || die "$!";
                        while ( my $line = <SC> ) {
                            if ( $line =~ /(.*) (.*)/ ) {
                                $myResults{$subDir}{'self_contain'} = $2; # récupération mfei
                            }
                        }
                        close SC;
                    }

                    #Récupération séquence et format Vienna
                    my $vienna_out = File::Spec->catfile( $subDir_full,
                                                       'outViennaTraited.txt' );
                    if ( -e $vienna_out )                  # si fichier existe
                    {
                        open( my $VIENNA, '<', $vienna_out ) || die "$!";
                        while ( my $line = <$VIENNA> ) {
                            if ( $line =~ /(.*)\t(.*)\t(.*)/ ) {
                                $myResults{$subDir}{'DNASequence'} = $2;
                                  ;    # récupération sequence
                                $myResults{$subDir}{'Vienna'} = $3;    # récupération Vienna

                            }
                        }
                        close $VIENNA;
                    }

                    #Récupération alignement avec mirBase
                    my $file_alignement =
                      File::Spec->catfile( $subDir_full, 'alignement.txt' );
                    if ( -e $file_alignement )                # si fichier existe
                    {
                        open( my $ALIGN, $file_alignement ) || die "$!";
                        my $align = "none";
                        while ( my $line = <$ALIGN> ) {
                            if ( $line =~ /^C4/ ) {
                                $align = $file_alignement;
                                last;
                            }
                        }
                        $myResults{$subDir}{'alignment'} = $align;
                        close $ALIGN or die ("Error closing $file_alignement: $!");
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
