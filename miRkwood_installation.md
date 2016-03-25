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

- Clone the miRkwood repository on the Inria Sequoia forge
    `svn checkout svn+ssh://scm.gforge.inria.fr/svnroot/sequoia/pipelineMiRNA/web/ mirkwood`

(Note that miRkwood uses SVN externals to fetch some provisionning dependencies)

- Run Vagrant at the miRkwood repository root
    `vagrant up`

- Vagrant will
    - download an Ubuntu ISO
    - mount the local development directory in the virtual machine
    - provision it using Ansible
    - set up network forwarding

This step can take up to 10 minutes.

Congratulations! A miRkwood instance is now running at <http://192.168.33.20/mirkwood>


Installation in a VM without Ansible
------------------------------------

If you are unlucky enough that you cannot install Ansible (this problem may
occur on Mac systems), you will have to install yourself all dependencies.

1. Create the VM

We still recommand to run miRkwood on an Ubuntu 12.04.

# Need details here to install Ubuntu on the VM.


2. Install dependencies

All these steps take place on the VM.

Run script install-dependencies.sh in super-user mode.


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
MIME::Lite (pour la version Web)
Inline::CPP (for smallRNAseq pipeline only)
Bio::DB::Fasta (for smallRNAseq pipeline only)

(pas utiles pour la version ligne de commande, tester pour la version web)
Image::Size
File::Size
File::Slurp
Archive::Zip


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
This program takes an output file from RNALfold and extract 
candidates stemloop from the secondary structures.
You can create a symbolic link to the program that suits your architecture,
or simply change the corresponding line in the programs.cfg file.

- Vienna package
The ViennaRNA Package consists of a C code library and several stand-alone 
programs for the prediction and comparison of RNA secondary structures.
http://www.tbi.univie.ac.at/RNA/
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
rnashuffles is a "homemade" program that computes the thermodynamic
stability.
Ensuire Python pip is installed.
Copy the sources where you want it to be (the sources are given in
/{miRkwood_path}/provisioning/roles/mirkwood-software/files/) and
then build it with pip.
pip install /path/rnashuffles

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
piccolo is a "homemade" program that performs alignments.
It is stored in the ${local_programs} directory.

3. Optional dependencies for abinitio pipeline
- tRNAscan-SE
This program searches for tRNA genes in genomic sequences.
Download the archive at http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz.
Extract the archive and replace '$(HOME)' by the path where you want to 
install it, then compile it.

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

- Pandoc
We use Pandoc to convert results report from ORG format to PDF format.
Check
http://pandoc.org/installing.html
for the installation guide.

You will also need LaTex, please install package 'texlive' with
your usual package manager.
