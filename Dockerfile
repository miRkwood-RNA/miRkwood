FROM ubuntu:14.04
MAINTAINER miRkwood team (mirkwood@univ-lille1.fr)


##### Check prerequisites ############################################
RUN apt-get -y update
RUN apt-get -y install make
RUN apt-get -y install gcc
RUN apt-get -y install g++
RUN apt-get -y install unzip
RUN apt-get -y install cpanminus
RUN apt-get -y install wget
RUN apt-get -y install default-jre --fix-missing


##### Install Perl modules ###########################################
RUN apt-get -y install libconfig-simple-perl
RUN apt-get -y install libyaml-libyaml-perl
RUN apt-get -y install libfile-which-perl
RUN apt-get -y install libmime-lite-perl
RUN apt-get -y install libimage-size-perl
RUN apt-get -y install libfile-type-perl
RUN apt-get -y install libfile-slurp-perl
RUN apt-get -y install libarchive-zip-perl

RUN cpanm LWP::UserAgent
RUN cpanm Bio::DB::Fasta
RUN cpanm Inline::CPP
RUN cpanm Log::Message::Simple


##### Install dependencies not provided with miRkwood ################

### Install bedtools
RUN apt-get -y install bedtools

### Install NCBI Blast+
RUN apt-get -y install ncbi-blast+

### Download VARNA jar
RUN wget --directory-prefix=/opt/ http://varna.lri.fr/bin/VARNAv3-91.jar

### Install tRNAscan-SE
# Download and extract tRNAscan-SE archive
RUN wget --directory-prefix=/tmp/ http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz
RUN tar xf /tmp/tRNAscan-SE.tar.gz --directory /tmp

# Update tRNAscan-SE paths
RUN sed -re 's/\$\(HOME\)\//\/opt\/tRNAscan-SE\//g' -i /tmp/tRNAscan-SE-1.3.1/Makefile

# Build and install tRNAscan-SE
RUN make --directory=/tmp/tRNAscan-SE-1.3.1/
RUN make install --directory=/tmp/tRNAscan-SE-1.3.1/

# Make relevant user/group and file mode
RUN chown -R www-data:www-data "/opt/tRNAscan-SE"
RUN chmod -R +x /opt/tRNAscan-SE/bin/tRNAscanSE/


##### Copy miRkwood code
RUN rm -Rf /home/mirkwood/
COPY ./cgi-bin/ /home/mirkwood/cgi-bin/
COPY ./provisioning/ /home/mirkwood/provisioning/


##### Install dependencies provided with miRkwood ####################

### Install Vienna package
RUN cp /home/mirkwood/provisioning/roles/mirkwood-software/files/vienna-rna_2.1.6-1_amd64.deb /tmp/vienna-rna_2.1.6-1_amd64.deb
RUN dpkg -i /tmp/vienna-rna_2.1.6-1_amd64.deb
RUN ln -s /usr/share/ViennaRNA/bin/b2ct /usr/bin/b2ct

### Install RNAstemloop
RUN rm -f /home/mirkwood/cgi-bin/programs/RNAstemloop
RUN ln -s /home/mirkwood/cgi-bin/programs/RNAstemloop-x86_64 /home/mirkwood/cgi-bin/programs/RNAstemloop

### Install miRdup
RUN wget --directory-prefix=/tmp/ http://www.cs.mcgill.ca/~blanchem/mirdup/miRdup_1.4.zip
RUN unzip -qq /tmp/miRdup_1.4.zip -d /opt/miRdup

### Install rnammer
# Install hmmer
RUN wget --directory-prefix=/tmp/  http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz
RUN mkdir /opt/hmmer-2.3.2/
RUN tar xf /tmp/hmmer-2.3.2.tar.gz --directory /opt
RUN cd /opt/hmmer-2.3.2/ && ./configure && make && make install
RUN ln -s /opt/hmmer-2.3.2/src/hmmsearch /usr/bin/hmmsearch23

# Install Perl dependency
RUN apt-get -y install libxml-simple-perl

# Copy RNAmmer archive in /opt
RUN cp /home/mirkwood/provisioning/roles/mirkwood-software/files/rnammer-1.2.src.tar.Z /tmp/rnammer-1.2.src.tar.Z

# Create RNAmmer directory
RUN mkdir /opt/RNAmmer

# Extract RNAmmer
RUN tar xf /tmp/rnammer-1.2.src.tar.Z --directory /opt/RNAmmer

# Add necessary module import to RNAmmer perl executable
RUN sed -re 's/(use Getopt::Long;)/use File::Basename;\n\1/' -i /opt/RNAmmer/rnammer

# Update self-path in RNAmmer perl executable
RUN sed -re 's/"\/usr\/cbs\/bio\/src\/rnammer-1.2"/dirname(__FILE__)/' -i /opt/RNAmmer/rnammer

# Update paths to HMMER in RNAmmer perl executable
RUN sed -re 's/\$HMMSEARCH_BINARY\s?=.*/$HMMSEARCH_BINARY="\/usr\/bin\/hmmsearch23";/' -i /opt/RNAmmer/rnammer

# Make relevant user/group
RUN chown -R www-data:www-data "/opt/RNAmmer"

### Install RNAshuffles
# Ensure Python pip is installed
RUN sudo apt-get -y install python-pip

# Copy RNAshuffles files
RUN cp -r /home/mirkwood/provisioning/roles/mirkwood-software/files/rnashuffles /opt/

# Build RNAshuffles
RUN pip install /opt/rnashuffles


##### Create symbolic links for programs #############################
RUN ln -s /opt/VARNAv3-91.jar /home/mirkwood/cgi-bin/programs/VARNA.jar
RUN ln -s /opt/miRdup /home/mirkwood/cgi-bin/programs/miRdup-1.4
RUN ln -s /opt/tRNAscan-SE /home/mirkwood/cgi-bin/programs/tRNAscan-SE
RUN ln -s /opt/RNAmmer /home/mirkwood/cgi-bin/programs/rnammer


##### Deploy miRkwood data ###########################################
RUN sh /home/mirkwood/cgi-bin/install-data.sh /home/mirkwood/cgi-bin/data

