SYNOPSIS
    miRkwood is an application that allows for the fast and easy identification of microRNAs. It is specifically designed for plant microRNAs.


INSTALL
    See file miRkwood_installation.txt.


REQUIRED FILES
    If you have annotation files for CDS, r/t/snoRNA and want to filter out
    CDS and/or r/t/snoRNA, make sure to store the corresponding files into
    the {miRkwood_path}/cgi-bin/data/annotations/ directory.
    These files should be named after your genome file. 
    E.g if your genome file name is "species.fa" or "species.fasta" you should
    name your annotation files "species_CDS.gff" and "species_tRNA_rRNA_snoRNA.gff".
    
    If you have a miRBase file for your species, miRkwood can detect known miRNAs.
    Make sure your miRBase file is stored in the {miRkwood_path}/cgi-bin/data/miRBase/
    directory and is named after your genome file.
    E.g if your genome file name is "species.fa" or "species.fasta" you should
    name your annotation file "species_miRBase.gff3".
    

USAGE
    miRkwood comes in two distinct pipelines, according to the input data type.
    
    -mirkwood.pl: scans a genomic sequence and finds all potential microRNA precursors.
        Input: a FASTA file.
    
    -mirkwood-bed.pl: analyses small RNA deep sequencing data and find all potential microRNAs.
        Input : a BED file.


OPTIONS
    -mirkwood.pl: ./mirkwood.pl [options] [FASTA files]
        --output
                Output directory. If non existing it will be created. The
                directory must be empty.

        --both-strands
                Scan both strands.

        --species-mask
                Mask coding regions against the given organism.

        --shuffles
                Compute thermodynamic stability (shuffled sequences).

        --filter-mfei
                Select only sequences with MFEI < -0.6.

        --filter-rrna
                Filter out ribosomal RNAs (using RNAmmer).

        --filter-trna
                Filter out tRNAs (using tRNAscan-SE).

        --align Flag conserved mature miRNAs (alignment with miRBase + miRdup).

        --no-varna
                Disable the structure generation using Varna.

        -help   Print a brief help message and exits.

        -man    Prints the manual page and exits.


    -mirkwood-bed.pl:./mirkwood-bed.pl [options] [BED file]
      Mandatory options:
        --genome
            Path to the genome (fasta format).

        --output
            Output directory. If non existing it will be created. The directory
            must be empty.

      Additional options:
        --shuffles
            Compute thermodynamic stability (shuffled sequences).

        --align
            Flag conserved mature miRNAs (alignment with miRBase + miRdup).

        --no-filter-mfei
            Don't filter out sequences with MFEI >= -0.6. Default : only keep
            sequences with MFEI < -0.6.

        --no-filter-CDS
            Don't filter out CDS. Default: if an annotation GFF file is
            available CDS are filtered out.

        --no-filter-r-t-snoRNA
            Don't filter out rRNA, tRNA, snoRNA. Default: if an annotation GFF
            file is available rRNA, tRNA, snoRNA are filtered out.

        --no-filter-multimapped
            Don't filter out multimapped reads. Default: reads that map at more
            than 5 positions are filtered out.

        --no-varna
            Disable the structure generation using Varna.

        -help
            Print a brief help message and exits.

        -man
            Prints the manual page and exits.
    

OUTPUT
    read_clouds: contains all text files for the candidates read clouds.
    
    results: contains all results files, in several formats (csv, fa, gff, html and txt)
    
    YML: contains all candidates data in YAML format.

    Depending on the options you chose for your job you may find some of
    the following files:
    
    your_bed_CDS.tar.gz: a BED containing all reads from your 
        input BED file corresponding to CDS.
    
    your_bed_tRNA_rRNA_snoRNA.tar.gz: a BED containing all reads from your 
        input BED file corresponding to t/r/snoRNA.

    your_bed_multimapped.tar.gz: a BED containing all reads from your 
        input BED file mapping at more than 5 positions.

    your_bed_miRNAs.tar.gz: a BED containing all reads from your 
        input BED file corresponding to miRNAs present in miRBase.

    your_bed_orphan_clusters.tar.gz: a BED containing all reads from your 
        input BED file that fall into a peak but that don't correspond to
        a valid miRNA candidate.

    your_bed_filtered.bed: a BED containing all reads from your 
        input BED file that have not been filtered out in one of the
        previous categories.





