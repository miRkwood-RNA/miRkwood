miRkwood
========

miRkwood is a computational pipeline for the identification of miRNAs and their hairpin precursors.

It is constituted of:
- a back-end in Perl, running many third-party programs, stored in `cgi-bin`
- a front-end in Perl/JavaScript, stored in `html`


miRkwood comes with a fair number of dependencies. It was tested on Ubuntu 12.04.
The easiest way to deploy miRkwood is to create a virtual machine (VM) and to use
the configuration management software Ansible.


Installation in a VM
--------------------

1. Create the VM

- Install VirtualBox <https://www.virtualbox.org/wiki/Downloads>

- Install Vagrant in its most recent version : <http://www.vagrantup.com/downloads.html>
  (tested on Vagrant 1.4.3 and Vagrant 1.6.5)

Vagrant is a tool to create and configure virtual development environments.
It can be considered a wrapper around virtualization software such as VirtualBox
and configuration management software such as Chef, Salt and Puppet − or Ansible in our case.


2. Install Ansible

Ansible is an IT automation tool. It can configure systems, deploy software, and orchestrate 
more advanced IT tasks such as continuous deployments or zero downtime rolling updates.

- Install Ansible in its most recent version (at least 1.6) <http://docs.ansible.com/intro_installation.html>
(tested with Ansible 1.6 and 1.7)


3. Install miRkwood dependencies

- Clone the miRkwood repository on the gitHub repository
    `git clone https://github.com/miRkwood-RNA/miRkwood.git`

- Run Vagrant at the miRkwood repository root
    `vagrant up`

- Vagrant will
    - download an Ubuntu ISO
    - mount the local development directory in the virtual machine
    - provision it using Ansible
    - set up network forwarding

This step can take up to 10 minutes.

Congratulations! A miRkwood instance is now running at <http://192.168.33.20/mirkwood>


Manual installation
-------------------

If you prefer to install miRkwood manually or if you don't want to use
a VM, please check that the following dependencies are installed 
according to your needs, and check that the paths in /cgi-bin/lib/programs.cfg
are corrects.
"${local_programs}" refers to /{miRkwood_path}/cgi-bin/programs/.


1. Mandatory dependencies

- Perl modules
Config::Simple
YAML::XS
File::Which
MIME::Lite (for Web version only)
Inline::CPP (for smallRNAseq pipeline only)
LWP::UserAgent (for smallRNAseq pipeline only)
Bio::DB::Fasta (for smallRNAseq pipeline only)


Inline::CPP may sometimes be difficult to install. In that case, please 
install first the following modules:
Inline::C
Inline
Pegex
ExtUtils::MakeMaker
File-ShareDir-Install
Capture::Tiny
ExtUtils::CppGuess


- RNAstemloop
RNAstemloop is an in-house program that takes an output file from
RNALfold and extract candidates stemloop from the secondary structures.
RNAstemloop is stored in ${local_programs}.
You can create a symbolic link to the program that suits your architecture,
or simply change the corresponding line in the programs.cfg file.

- Vienna package
The ViennaRNA Package consists of a C code library and several stand-alone 
programs for the prediction and comparison of RNA secondary structures.
See http://www.tbi.univie.ac.at/RNA/ for the installation.
Tested with version 2.1.6-1.

Make sure that the following programs are correctly installed:
b2ct
RNAfold
RNAlfold
RNAeval

Be careful that RNAfold, RNAlfold and RNAeval are usually installed
in /usr/bin, and b2ct is usually installed in /usr/share/ViennaRNA/bin.

- bedtools
Install package 'bedtools' with your usual package manager,
for instance with
`sudo apt-get install bedtools`
or
`sudo yum install bedtools`


2. Optional dependencies for both pipelines

- miRdup
miRdup is a tool for the validation of pre-miRNAs predictions.
Download the archive at
http://www.cs.mcgill.ca/~blanchem/mirdup/
(tested with version 1.2)
However, if you don't want or don't manage to install it, you can skip
it. It will only affect the quality score of some candidates.

- rnashuffles
rnashuffles is an in-house program that computes the thermodynamic
stability.
Ensuire Python pip is installed.
Copy the sources where you want it to be (the sources are given in
/{miRkwood_path}/provisioning/roles/mirkwood-software/files/) and
then build it with pip.
`pip install /path/rnashuffles`

- VARNA (Visualization Applet for RNA secondary structure)
Make sure that the Java Runtime Environment is installed.
Install package 'default-jre' with your usual package manager,
for instance with
`sudo apt-get install default-jre`
or
`sudo yum install default-jre`

You can download VARNA jar on http://varna.lri.fr/bin/VARNAv3-91.jar and
then create a symbolic link or change the corresponding line in the 
programs.cfg file.

- piccolo
piccolo is an in-house program that performs alignments.
It is stored in the ${local_programs} directory.


3. Optional dependencies for abinitio pipeline

- tRNAscan-SE
This program searches for tRNA genes in genomic sequences.
Download the archive at http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz.
Extract the archive and edit the Makefile to replace '$(HOME)' by the 
path where you want to install it, then compile it.

- rnammer
This program predicts 5s/8s, 16s/18s and 23s/28s ribosomal
RNA in genomic sequences.
# Install hmmer
wget --directory-prefix=/tmp/  http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz
tar xf /tmp/hmmer-2.3.2.tar.gz --directory /opt
cd /opt/hmmer-2.3.2/
./configure --enable-threads
make --directory=/opt/hmmer-2.3.2/
ln -s /opt/hmmer-2.3.2/src/hmmsearch /usr/bin/hmmsearch23

# Install Perl dependency
sudo apt-get install libxml-simple-perl

# Copy RNAmmer archive in /opt
cp $ROOT_PATH/provisioning/roles/mirkwood-software/files/rnammer-1.2.src.tar.Z /tmp/rnammer-1.2.src.tar.Z

# Create RNAmmer directory
mkdir /opt/RNAmmer

# Extract RNAmmer
tar xf /tmp/rnammer-1.2.src.tar.Z --directory /opt/RNAmmer

# Add necessary module import to RNAmmer perl executable
sed -re 's/(use Getopt::Long;)/use File::Basename;\n\1/' -i /opt/RNAmmer/rnammer

# Update self-path in RNAmmer perl executable
sed -re 's/"\/usr\/cbs\/bio\/src\/rnammer-1.2"/dirname(__FILE__)/' -i /opt/RNAmmer/rnammer

# Update paths to HMMER in RNAmmer perl executable
sed -re 's/\$HMMSEARCH_BINARY\s?=.*/$HMMSEARCH_BINARY="\/usr\/bin\/hmmsearch23"/' -i /opt/RNAmmer/rnammer

# Make relevant user/group
chown -R www-data:www-data "/opt/RNAmmer"

- blastX
BLAST finds regions of local similarity between sequences.
We use it to mask CDS regions in the sequences given by the user.
Install package 'ncbi-blast+' with your usual package manager,
for instance with
`sudo apt-get install ncbi-blast+`
or
`sudo yum install ncbi-blast+`


Usage with Docker
-----------------

If you choose to use Docker you won't have to deal with installing 
miRkwood dependencies. All you will have to install is Docker
(https://docs.docker.com/engine/installation/).

A Docker image for miRkwood is stored at: https://hub.docker.com/r/iguigon/mirkwood/
You can download miRkwood image with:
sudo docker pull iguigon/mirkwood


Usage:

1. Abinitio pipeline:

Here is an example of command line for abinitio pipeline:

sudo docker run -v /your/mapping/repertory/:/home/mapping/ \
    iguigon/mirkwood \
    ./mirkwood.pl \
    --input /home/mapping/input/sequenceSomething.fas \
    --output /home/mapping/output/results_abinitio \
    --filter-rrna --filter-trna --varna

Type the following command for detailed help on the command options:
sudo docker run iguigon/mirkwood perl ./mirkwood.pl --help


2. SmallRNAseq pipeline:

miRkwood smallRNAseq pipeline takes as input a BED file with the following syntax:

1    18092    18112    SRR051927.5475072    1    -
1    18094    18118    SRR051927.2544175    2    +
1    18096    18119    SRR051927.3033336    1    +
1    18100    18124    SRR051927.172198     9    +

In this file, each line is a unique read.
The fields are, from left to right: name of the chromosome, starting position,
ending position, read identifier, number of occurrences of the read in the data, strand.
Positions follow the BED numbering convention: the first base of the chromosome is
considered position 0 (0-based position) and the feature does not include the stop position.

You can convert a BAM file into the needed BED format with our custom script:
sudo docker run -v /your/mapping/repertory/:/home/mapping/ \
    iguigon/mirkwood \
    ./mirkwood-bam2bed.pl \
    --in /home/mapping/input.bam \
    --bed /home/mapping/output.bed \
    --min 18 \
    --max 25

Type the following command for detailed help on the command options:
sudo docker run iguigon/mirkwood ./mirkwood-bam2bed.pl --help


Here is an example of command line for smallRNAseq pipeline:

sudo docker run -v /your/mapping/repertory/:/home/mapping/ \
    iguigon/mirkwood \
    ./mirkwood-bed.pl \
    --input /home/mapping/input/sample.bed \
    --output /home/mapping/output/results_smallRNAseq \
    --genome /home/mapping/input/my_genome.fasta \
    --mirbase /home/mapping/input/my_mirbase_file.gff3 \
    --gff /home/mapping/input/my_annotations_file.gff \
    --min-read-positions-nb 0 --max-read-positions-nb 5 --align

Type the following command for detailed help on the command options:
sudo docker run iguigon/mirkwood ./mirkwood-bed.pl --help


With the -v /your/mapping/repertory/:/home/mapping/ parameter,
Docker will mount your local folder /your/mapping/repertory/ into the container under /home/mapping/.
Make sure you have stored all your needed input files in this folder.
Results files will be stored here as well.


