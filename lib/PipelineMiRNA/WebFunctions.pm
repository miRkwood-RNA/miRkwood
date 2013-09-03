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

Usage:
my %results = PipelineMiRNA::WebFunctions->get_structure_for_jobID($jobId);

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
                    my %candidate = retrieve_candidate_information($job, $dir, $subDir);
                    $myResults{$subDir} = \%candidate;
                }
            }
        }
    }
    return %myResults;
}

=method retrieve_candidate_information

Get the result for a given candidate

Arguments:
- $job - the job identifier
- $dir - the sequence name
- $subDir - the candidate name
=cut

sub retrieve_candidate_information {
    my @args = @_;
    my $job = shift @args;
    my $dir = shift @args;
    my $subDir = shift @args;

    my $job_dir = PipelineMiRNA::Paths->get_absolute_path($job);
    my $subDir_full = File::Spec->catdir( $job_dir, $dir, $subDir );
    my %result = ();

    $result{'name'} = $dir;    #récupération nom séquence
    my @position = split( /__/, $subDir );
    $result{'position'} = $position[1]; # récupération position
      ;    # récupération noms + positions
    my $pvalue =
      File::Spec->catfile( $subDir_full, 'pvalue.txt' );
    if ( -e $pvalue )    # si fichier existe
    {
        $result{'p_value'} = PipelineMiRNA::Parsers::parse_pvalue($pvalue);
    }

    #Récupération valeur MFEI
    my $mfei_out =
      File::Spec->catfile( $subDir_full, 'outMFEI.txt' );
    if ( -e $mfei_out )                 # si fichier existe
    {
        my @mfeis = PipelineMiRNA::Parsers::parse_mfei($mfei_out);
        $result{'mfei'} = $mfeis[0];
        $result{'mfe'} = $mfeis[1];
        $result{'amfe'} = $mfeis[2];
    }

    #Récupération valeur self contain
    my $selfcontain_out =
      File::Spec->catfile( $subDir_full, 'selfContain.txt' );
    if ( -e $selfcontain_out )
    {Dumper
        $result{'self_contain'} = PipelineMiRNA::Parsers::parse_selfcontain($selfcontain_out);
    }

    #Récupération séquence et format Vienna
    my $vienna_out = File::Spec->catfile( $subDir_full,
                                       'outViennaTraited.txt' );
    if ( -e $vienna_out )                  # si fichier existe
    {
        my @vienna_res = PipelineMiRNA::Parsers::parse_vienna($vienna_out);
        $result{'DNASequence'} = $vienna_res[0];
        $result{'Vienna'} = $vienna_res[1];
    }

    #Récupération alignement avec mirBase
    my $file_alignement = PipelineMiRNA::Paths->get_absolute_path($job,  $dir, $subDir, 'alignement.txt');
    if ( -e $file_alignement )                # si fichier existe
    {
        #my $align_res = PipelineMiRNA::Parsers::parse_alignment($file_alignement);
        $result{'alignment'} = $file_alignement;
    }
    my $image_path = File::Spec->catfile($job,  $dir, $subDir, 'image.png');
    $result{'image'} = $image_path;

    return %result;
}

=method resultstruct2csv

Convert the results structure to CSV

=cut

sub resultstruct2csv {
    my ( $self, @args ) = @_;
    my $results = shift @args;
	my @tab = shift @args;
	my %results = %{$results};
    my @headers = ('name', 'position', 'mfei', 'mfe', 'amfe', 'p_value', 'self_contain', 'Vienna', 'DNASequence');
	my $result = join( ',', @headers ) . "\n";
    while ( my ($key, $value) = each %results )
    {
    	if (  $key ~~ \@tab ) 
    	{
 	        for my $header (@headers)
 	        {
	            $result .= "${$value}{$header},";
	        }
	        $result .= "\n";
    	}
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

    my @headers1 = ('name', 'position', 'mfei', 'mfe', 'amfe', 'p_value', 'self_contain', 'alignment');
    my @headers2 = ('Vienna', 'DNASequence');

    while ( my ($key, $value) = each %results )
    {
        $result .= "<Sequence";
        for my $header (@headers1){
            $result .= " $header='${$value}{$header}'";
        }
        my $img = PipelineMiRNA::Paths->get_server_path(${$value}{'image'});
        $result .= " image='$img'";
        for my $header (@headers2){
            $result .= " $header='${$value}{$header}'";
        }
        $result .= "></Sequence>\n";
    }
    $result .= "</results>";
    return $result;
}


=method make_Vienna_viz

Make a nicer Vienna display by cutting too long lines.

Usage:
my $string = make_Vienna_viz($Vienna, $DNASequence)

=cut

sub make_Vienna_viz {
    my ($self, @args) = @_;
    my $Vienna = shift @args;
    my $DNASequence = shift @args;

    my $viennaString   = q{};
    my $sequenceString = q{};
    my $string         = q{};
    for ( 1 .. length($Vienna) ) {

        $viennaString   .= substr $Vienna,      $_ - 1, 1;
        $sequenceString .= substr $DNASequence, $_ - 1, 1;
        if ( $_ % 50 == 0 ) {

            $string .= $viennaString . "\n" . $sequenceString . "\n\n";
            $viennaString   = q{};
            $sequenceString = q{};
        }
        if ( ( $viennaString ne q{} ) && ( $_ == length($Vienna) ) ) {
            $string .= $viennaString . "\n" . $sequenceString . "\n\n";
        }
    }
    return $string
}

1;
