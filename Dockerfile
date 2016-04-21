FROM ubuntu:14.04
MAINTAINER miRkwood team (mirkwood@univ-lille1.fr)


RUN apt-get -y update && apt-get -y install cpanminus
RUN apt-get -y install make
RUN apt-get -y install gcc
RUN apt-get -y install g++


# Install Perl modules
RUN cpanm Config::Simple
RUN cpanm YAML::XS
RUN cpanm File::Which
RUN cpanm MIME::Lite
RUN cpanm LWP::UserAgent
RUN cpanm Bio::DB::Fasta

# RUN apt-get -y install libextutils-cppguess-perl

# RUN cpanm Pegex
# RUN cpanm ExtUtils::MakeMaker
# RUN cpanm File::ShareDir::Install
# RUN cpanm Capture::Tiny
# RUN cpanm ExtUtils::CppGuess

# RUN cpanm Math::Simple
# RUN cpanm Boo::Far

# RUN apt-get -y install libinline-perl
# RUN cpanm --force Inline::C 

# RUN cpanm Inline

# RUN cpanm --force Inline::CPP

# Copy miRkwood code
COPY . /home/mirkwood/

# Install Vienna package
RUN cp /home/mirkwood/provisioning/roles/mirkwood-software/files/vienna-rna_2.1.6-1_amd64.deb /tmp/vienna-rna_2.1.6-1_amd64.deb
RUN dpkg -i /tmp/vienna-rna_2.1.6-1_amd64.deb
RUN ln -s /usr/share/ViennaRNA/bin/b2ct /usr/bin/b2ct

# install RNAstemloop
# RUN rm -f /home/mirkwood/cgi-bin/programs/RNAstemloop
# RUN mv /home/mirkwood/cgi-bin/programs/RNAstemloop-x86_64 /home/mirkwood/cgi-bin/programs/RNAstemloop

RUN sed -i -re 's/RNAstemloop/RNAstemloop-x86_64/g'  /home/mirkwood/cgi-bin/lib/programs.cfg
# RUN ln -s /home/mirkwood/cgi-bin/programs/RNAstemloop-x86_64 /home/mirkwood/cgi-bin/programs/RNAstemloop




# CMD ["--help"]
# ENTRYPOINT ["perl", "-I/home/miRkwood/cgi-bin/lib", "/home/miRkwood/cgi-bin/bin/mirkwood.pl"]


