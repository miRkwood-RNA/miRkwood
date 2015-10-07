package miRkwood::ResultsExporterMaker;

# ABSTRACT: Builder for Results exporter

use strict;
use warnings;

# A PseudoXMLExporter module was created (by Jean-Fred ?)
# It was never used so this module, the function candidate_as_pseudoXML
# and the corresponding tests are deleted in r13371 and r13464.
# If you want to use them, be aware that this code is out of date.

sub make_fasta_results_exporter {
    my ($self, @args) = @_;
    my $mirna_type = shift @args;
    require miRkwood::ResultsExporter::FastaExporter;
    return miRkwood::ResultsExporter::FastaExporter->new( $mirna_type );
}

sub make_dotbracket_results_exporter {
    my ($self, @args) = @_;
    my $mirna_type = shift @args;
    require miRkwood::ResultsExporter::DotBracketExporter;
    return miRkwood::ResultsExporter::DotBracketExporter->new( $mirna_type );
}

sub make_gff_results_exporter {
    my ($self, @args) = @_;
    my $mirna_type = shift @args;
    require miRkwood::ResultsExporter::GFFExporter;
    return miRkwood::ResultsExporter::GFFExporter->new( $mirna_type );
}

sub make_csv_results_exporter {
    my ($self, @args) = @_;
    my $pipeline_type = shift @args;
    my $mirna_type    = shift @args;
    require miRkwood::ResultsExporter::CSVExporter;
    require miRkwood::ResultsExporter::CSVExporterGenomic;
    require miRkwood::ResultsExporter::CSVExporterSmallRNAseqKnown;
    require miRkwood::ResultsExporter::CSVExporterSmallRNAseqNew;
    if ( $pipeline_type eq 'abinitio' ){
        return miRkwood::ResultsExporter::CSVExporterGenomic->new( $mirna_type );
    }
    else{
        if ( $mirna_type eq 'novel_miRNA' ){
            return miRkwood::ResultsExporter::CSVExporterSmallRNAseqNew->new( $mirna_type );
        }
        else{
            return miRkwood::ResultsExporter::CSVExporterSmallRNAseqKnown->new( $mirna_type );
        }
    }

}

sub make_opendocument_results_exporter {
    my ($self, @args) = @_;
    my $mirna_type = shift @args;
    require miRkwood::ResultsExporter::OpenDocumentExporter;
    return miRkwood::ResultsExporter::OpenDocumentExporter->new( $mirna_type );
}

sub make_html_results_exporter {
    my ($self, @args) = @_;
    my $pipeline_type = shift @args;
    my $mirna_type    = shift @args;
    require miRkwood::ResultsExporter::HTMLExporter;
    require miRkwood::ResultsExporter::HTMLExporterGenomic;
    require miRkwood::ResultsExporter::HTMLExporterSmallRNAseqKnown;
    require miRkwood::ResultsExporter::HTMLExporterSmallRNAseqNew;
    if ( $pipeline_type eq 'abinitio' ){
        return miRkwood::ResultsExporter::HTMLExporterGenomic->new( $mirna_type );
    }
    else{
        if ( $mirna_type eq 'novel_miRNA' ){
            return miRkwood::ResultsExporter::HTMLExporterSmallRNAseqNew->new( $mirna_type );
        }
        else{
            return miRkwood::ResultsExporter::HTMLExporterSmallRNAseqKnown->new( $mirna_type );
        }
    }

}

sub make_reads_clouds_results_exporter {
    my ($self, @args) = @_;
    my $mirna_type = shift @args;
    require miRkwood::ResultsExporter::ReadsExporter;
    return miRkwood::ResultsExporter::ReadsExporter->new( $mirna_type );
}



1;
