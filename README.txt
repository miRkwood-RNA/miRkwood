SYNOPSIS

    miRkwood is an application that allows for the fast and easy identification of microRNAs. It is specifically designed for plant microRNAs.


INSTALL

    See file INSTALL.md.


USAGE

    miRkwood comes in two distinct pipelines, according to the input data type.

    -mirkwood.pl (abinitio pipeline): scans a genomic sequence and finds all potential microRNA precursors.
        Input: a FASTA file.

    -mirkwood-bed.pl (smallRNAseq pipeline): analyses small RNA deep sequencing data and find all potential microRNAs.
        Input : a BED file.


OPTIONS

    -mirkwood.pl: perl -I/{miRkwood_path}/cgi-bin/lib/ mirkwood.pl [options]
          Mandatory options:
            --input
                Path to the fasta file.

            --output
                Output directory. If non existing it will be created. The directory
                must be empty.

          Additional options:
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

            --align
                Flag conserved mature miRNAs (alignment with miRBase + miRdup).

            --varna
                Allow the structure generation using Varna.

            --help
                Print a brief help message and exits.

            --man
                Prints the manual page and exits.


    -mirkwood-bed.pl: perl -I/{miRkwood_path}/cgi-bin/lib/ mirkwood-bed.pl [options]
          Mandatory options:
            --input
                Path to the BED file (created with our script mirkwood-bam2bed.pl).

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

            --mirbase
                If you have a gff file containing known miRNAs for this assembly,
                use this option to give the path to this file.

            --gff
                List of annotation files (gff or gff3 format). Reads matching with
                an element of these files will be filtered out. For instance you can
                filter out CDS by providing a suitable GFF file.

            --no-filter-bad-hairpins
                By default the candidates with a quality score of 0 and no
                conservation are discarded from results and are stored in a BED
                file. Use this option to keep all results.

            --min-read-positions-nb
                Minimum number of positions for each read to be kept. Default : 0.

            --max-read-positions-nb
                Maximum number of positions for each read to be kept. Default : 5
                (reads that map at more than 5 positions are filtered out).

            --varna
                Allow the structure generation using Varna.

            --help
                Print a brief help message and exits.

            --man
                Prints the manual page and exits.


OUTPUT

    For both pipelines:

        alignments : folder containing all alignments files
            (only if option --align is on).

        images: folder containing images created by VARNA
            (only if option --varna is on).

        results: folder containing all results files, in several 
            formats (csv, fa, gff, html and txt).

        sequences: folder containing sequences for each candidate 
            in fasta and dotbracket format, alternatives sequences 
            if they exist and optimal structure if it is different 
            from the stemloop structure.

        YML: folder containing all candidates data in YAML format.

        basic_candidates.yml: contains a summary of all candidates
            with basic informations (this file is needed to create
            the results files).

        log.log: log file (hey, what did you expect?)

        run_options.cfg: config file with the chosen options.

    ab initio pipeline only:

        masks: folder containing results of BlastX, rnammer and tRNAscan-SE.

        input_sequences.fas: your sequences.

    smallRNAseq pipeline only:

        read_clouds: folder containing all text files for the candidates 
            read clouds.

        bed_sizes.txt: tabulated file with the number of reads in each BED file.

        summary.txt: contains a summary of your options and of results.


        Depending on the options you chose for your job you may find 
        some of the following files:

            your_bed_your_GFF.tar.gz: a compressed BED containing all reads matching
                to features from your GFF file, for each GFF file that you
                provided.

            your_bed_multimapped.tar.gz: a compressed BED containing all reads from your 
                input BED file mapping at less than --min-read-positions-nb positions
                or more than --max-read-positions-nb positions.

            your_bed_miRNAs.tar.gz: a compressed BED containing all reads from your 
                input BED file corresponding to miRNAs present in miRBase.

            your_bed_orphan_clusters.tar.gz: a compressed BED containing all reads from your 
                input BED file that fall into a peak but that don't correspond to
                a valid miRNA candidate.

            your_bed_orphan_hairpins.tar.gz: a compressed BED containing all candidates
                with a quality score of 0 and no conservation. By default 
                these candidates are excluded from final results, but you can
                change this behaviour with flag option --no-filter-bad-hairpins.

            your_bed_filtered.bed: a BED containing all reads from your 
                input BED file that have not been filtered out in one of the
                previous categories.


