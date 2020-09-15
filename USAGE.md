# General usage


## Abinitio pipeline:

Here is an example of command line for abinitio pipeline:

`./mirkwood.pl \`
`--input sequenceSomething.fas \`
`--output results_abinitio \`
`--filter-rrna \`
`--filter-trna \`
`--varna`

Type the following command for detailed help on the command options:
`./mirkwood.pl --help`


## SmallRNAseq pipeline:

miRkwood smallRNAseq pipeline takes as input a BED file with the following syntax:

`1	18092	18112	AAACGTGTAGAGAGAGACTCA	1	-`
`1	18094	18118	GATTCTTTTGTTTGCCACT	2	+`
`1	18096	18119	TCGATAGGATCAAGTACATCT	1	+`
`1	18100	18124	AAGAAGAAAAAGAAGAAGAAGAAG	9	+`

In this file, each line is a unique read.
The fields are, from left to right: name of the chromosome, starting position,
ending position, read identifier, number of occurrences of the read in the data, strand.
Positions follow the BED numbering convention: the first base of the chromosome is
considered position 0 (0-based position) and the feature does not include the stop position.

You can convert a BAM file into the needed BED format with our custom script:
`./mirkwood-bam2bed.pl \`
`--in input.bam \`
`--bed output.bed \`
`--min 18 \`
`--max 25`

Type the following command for detailed help on the command options:
`./mirkwood-bam2bed.pl --help`


Here is an example of command line for smallRNAseq pipeline:

`./mirkwood-bed.pl \`
`--input sample.bed \`
`--output results_smallRNAseq \`
`--genome my_genome.fasta \`
`--mirbase my_mirbase_file.gff3 \`
`--gff my_annotations_file1.gff \`
`--gff my_annotations_file2.gff \`
`--min-read-positions-nb 0 \`
`--max-read-positions-nb 5 \`
`--align`

Type the following command for detailed help on the command options:
`./mirkwood-bed.pl --help`


# Usage with Docker

miRkwood scripts are stored in `/home/mirkwood/cgi-bin/bin`.

There are several ways to use a Docker container.


## Use the same container


### Create a Docker container

`sudo docker run -dit -v /your/mapping/repertory/:/MAPPING/ --name mirkwood iguigon/mirkwood:latest`

With the `-v /your/mapping/repertory/:/MAPPING/` parameter,
Docker will mount your local folder `/your/mapping/repertory/` into the container under `/MAPPING/`.
Make sure you have stored all your needed input files in this folder.


### Launchs jobs

Example with SmallRNAseq pipeline:


### Outside the container

`sudo docker exec -it mirkwood mirkwood-bed.pl \`
`--input /MAPPING/sample.bed \`
`--output /MAPPING/results_smallRNAseq \`
`--genome /MAPPING/my_genome.fasta \`
`--mirbase /MAPPING/my_mirbase_file.gff3 \`
`--gff /MAPPING/my_annotations_file1.gff \`
`--gff /MAPPING/my_annotations_file2.gff \`
`--min-read-positions-nb 0 \`
`--max-read-positions-nb 5 \`
`--align`


## Use a new container for each job

`sudo docker run \`
`-v /your/mapping/repertory/:/MAPPING/ \`
`iguigon/mirkwood:latest mirkwood-bed.pl \`
`--input /MAPPING/sample.bed \`
`--output /MAPPING/results_smallRNAseq \`
`--genome /MAPPING/my_genome.fasta \`
`--mirbase /MAPPING/my_mirbase_file.gff3 \`
`--gff /MAPPING/my_annotations_file1.gff \`
`--gff /MAPPING/my_annotations_file2.gff \`
`--min-read-positions-nb 0 \`
`--max-read-positions-nb 5 \`
`--align`
