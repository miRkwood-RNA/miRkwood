package miRkwood::ResultsExporterMaker;

# ABSTRACT: Builder for Results exporter

use strict;
use warnings;

sub make_fasta_results_exporter {
    require miRkwood::ResultsExporter::FastaExporter;
    return miRkwood::ResultsExporter::FastaExporter->new();
}

sub make_dotbracket_results_exporter {
    require miRkwood::ResultsExporter::DotBracketExporter;
    return miRkwood::ResultsExporter::DotBracketExporter->new();
}

sub make_gff_results_exporter {
    require miRkwood::ResultsExporter::GFFExporter;
    return miRkwood::ResultsExporter::GFFExporter->new();
}

sub make_csv_results_exporter {
    my ($self, @args) = @_;
    my $html_exporter = shift @args;
    require miRkwood::ResultsExporter::CSVExporter;
    require miRkwood::ResultsExporter::CSVExporterGenomic;
    require miRkwood::ResultsExporter::CSVExporterSmallRNAseqKnown;
    require miRkwood::ResultsExporter::CSVExporterSmallRNAseqNew;
    if ( $html_exporter eq 'smallRNAseqKnown' ){
       return miRkwood::ResultsExporter::CSVExporterSmallRNAseqKnown->new();
    }
    elsif( $html_exporter eq 'smallRNAseqNew' ){
        return miRkwood::ResultsExporter::CSVExporterSmallRNAseqNew->new();
    }
    else{
        return miRkwood::ResultsExporter::CSVExporterGenomic->new();
    }
}

sub make_opendocument_results_exporter {
    require miRkwood::ResultsExporter::OpenDocumentExporter;
    return miRkwood::ResultsExporter::OpenDocumentExporter->new();
}

sub make_pseudoxml_results_exporter {
    require miRkwood::ResultsExporter::PseudoXMLExporter;
    return miRkwood::ResultsExporter::PseudoXMLExporter->new();
}

sub make_html_results_exporter {
    my ($self, @args) = @_;
    my $html_exporter = shift @args;
    require miRkwood::ResultsExporter::HTMLExporter;
    require miRkwood::ResultsExporter::HTMLExporterGenomic;
    require miRkwood::ResultsExporter::HTMLExporterSmallRNAseqKnown;
    require miRkwood::ResultsExporter::HTMLExporterSmallRNAseqNew;
    if ( $html_exporter eq 'smallRNAseqKnown' ){
       return miRkwood::ResultsExporter::HTMLExporterSmallRNAseqKnown->new();
    }
    elsif( $html_exporter eq 'smallRNAseqNew' ){
        return miRkwood::ResultsExporter::HTMLExporterSmallRNAseqNew->new();
    }
    else{
        return miRkwood::ResultsExporter::HTMLExporterGenomic->new();
    }
}

sub make_reads_clouds_results_exporter {
    my ($self, @args) = @_;
    my $mirna_type = shift @args;
    require miRkwood::ResultsExporter::ReadsExporter;
    return miRkwood::ResultsExporter::ReadsExporter->new( $mirna_type );
}



1;
