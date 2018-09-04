#!/bin/sh

# Script to install all dependencies for miRkwood on Ubuntu
# if the install with Ansible didn't work.
# Takes as argument the path to your miRkwood directory.
# You need to run this script in super-user mode.



########## Get options ###################################################################
WEB_ONLY=
CLI_ONLY=
GIVEN_PATH=
ROOT_PATH=$(dirname $(readlink -f $0));

while getopts 'c' OPTION
do
    case $OPTION in
    c)  CLI_ONLY=1
        ;;
    ?)  printf "Usage : %s [-c]\n" "$0"
        printf "   -c : install only requirements for CLI version\n"
        exit 2
        ;;
    esac
done


cd $ROOT_PATH


########## Look for the architecture #####################################################
ARCH=`uname -m`


########## Update & upgrade ##############################################################
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Update & upgrade"

sudo apt-get -y update
sudo apt-get -y upgrade


########## Prerequisites #################################################################
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Check prerequisites"

##### Make sure 'make' is installed
MAKE=/usr/bin/make
if [ ! -e $MAKE ]
then
    echo "..... Install make ........................................."
    sudo apt-get -y install make
else
    echo "..... make is already installed............................."
fi


##### Make sure compiler 'g++' is installed
GCOMP=/usr/bin/g++
if [ ! -e $GCOMP ]
then
    echo "..... Install g++ .........................................."
    sudo apt-get -y install g++
else
    echo "..... g++ is already installed ............................."
fi


##### Make sure 'unzip' is installed
UNZIP=/usr/bin/unzip
if [ ! -e $UNZIP ]
then
    echo "..... Install unzip ........................................"
    sudo apt-get -y install unzip
else
    echo "..... unzip is already installed ..........................."
fi


##### Make sure 'wget' is installed
WGET=/usr/bin/wget
if [ ! -e $WGET ]
then
    echo "..... Install wget ........................................."
    sudo apt-get -y install wget
else
    echo "..... wget is already installed ............................"
fi



########## Install dependencies for both web version and CLI version of miRkwood #########
echo
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "Install dependencies for both web version and CLI version of miRkwood"

##### Install Vienna package
echo "..... Install Vienna package ..............................."
cp provisioning/roles/mirkwood-software/files/vienna-rna_2.1.6-1_amd64.deb /tmp/vienna-rna_2.1.6-1_amd64.deb
sudo dpkg -i /tmp/vienna-rna_2.1.6-1_amd64.deb
ln -s /usr/share/ViennaRNA/bin/b2ct /usr/bin/b2ct


##### Install build-essential
echo "..... Install build-essential .............................."
sudo apt-get -y install build-essential


##### Install bedtools
echo "..... Install bedtools ....................................."
sudo apt-get -y install bedtools


##### Install NCBI Blast+
echo "..... Install ncbi-blast+ .................................."
sudo apt-get -y install ncbi-blast+


##### Install miRdup
echo "..... Install miRdup ......................................."
wget --directory-prefix=/tmp/ http://www.cs.mcgill.ca/~blanchem/mirdup/miRdup_1.4.zip
unzip -qq /tmp/miRdup_1.4.zip -d /opt/miRdup


##### Install VARNA
echo "..... Install VARNA ........................................"
# Ensure the Java Runtime Environment is installed
sudo apt-get -y install default-jre

# Download VARNA jar
sudo wget --directory-prefix=/opt/ http://varna.lri.fr/bin/VARNAv3-91.jar


##### Install tRNAscan-SE (only needed for ab initio pipeline)
echo "..... Install tRNAscan-SE .................................."
# Download and extract tRNAscan-SE archive
wget --directory-prefix=/tmp/ http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz
tar xf /tmp/tRNAscan-SE.tar.gz --directory /tmp

# Update tRNAscan-SE paths
sed -re 's/\$\(HOME\)\//\/opt\/tRNAscan-SE\//g' -i /tmp/tRNAscan-SE-1.3.1/Makefile

# Build and install tRNAscan-SE
make --directory=/tmp/tRNAscan-SE-1.3.1/
make install --directory=/tmp/tRNAscan-SE-1.3.1/

# Make relevant user/group and file mode
chown -R www-data:www-data "/opt/tRNAscan-SE"
sudo chmod -R +x /opt/tRNAscan-SE/bin/tRNAscanSE/


##### Install RNAmmer (only needed for ab initio pipeline)
echo "..... Install RNAmmer ......................................"
# Install hmmer
wget --directory-prefix=/tmp/  http://eddylab.org/software/hmmer/2.3.2/hmmer-2.3.2.tar.gz
sudo mkdir /opt/hmmer-2.3.2/
sudo tar xf /tmp/hmmer-2.3.2.tar.gz --directory /opt
cd $ROOT_PATH
make --directory=/opt/hmmer-2.3.2/
sudo ln -s /opt/hmmer-2.3.2/src/hmmsearch /usr/bin/hmmsearch23

# Install Perl dependency
sudo apt-get -y install libxml-simple-perl

# Copy RNAmmer archive in /opt
cp $ROOT_PATH/provisioning/roles/mirkwood-software/files/rnammer-1.2.src.tar.Z /tmp/rnammer-1.2.src.tar.Z

# Create RNAmmer directory
sudo mkdir /opt/RNAmmer

# Extract RNAmmer
sudo tar xf /tmp/rnammer-1.2.src.tar.Z --directory /opt/RNAmmer

# Add necessary module import to RNAmmer perl executable
sudo sed -re 's/(use Getopt::Long;)/use File::Basename;\n\1/' -i /opt/RNAmmer/rnammer

# Update self-path in RNAmmer perl executable
sudo sed -re 's/"\/usr\/cbs\/bio\/src\/rnammer-1.2"/dirname(__FILE__)/' -i /opt/RNAmmer/rnammer

# Update paths to HMMER in RNAmmer perl executable
sudo sed -re 's/\$HMMSEARCH_BINARY\s?=.*/$HMMSEARCH_BINARY="\/usr\/bin\/hmmsearch23";/' -i /opt/RNAmmer/rnammer

# Make relevant user/group
sudo chown -R www-data:www-data "/opt/RNAmmer"


##### Install RNAshuffles
echo "..... Install RNAshuffles .................................."
# Ensure Python pip is installed
sudo apt-get -y install python-pip

# Copy RNAshuffles files
cp -r $ROOT_PATH/provisioning/roles/mirkwood-software/files/rnashuffles /opt/

# Build RNAshuffles
pip install /opt/rnashuffles


##### Install Perl dependencies
echo "..... Install Perl dependencies via apt ...................."
sudo apt-get -y install libconfig-simple-perl
sudo apt-get -y install libyaml-libyaml-perl
sudo apt-get -y install libfile-which-perl
sudo apt-get -y install libmime-lite-perl
sudo apt-get -y install libimage-size-perl
sudo apt-get -y install libfile-type-perl
sudo apt-get -y install libfile-slurp-perl
sudo apt-get -y install libarchive-zip-perl

# Ensure cpanm is available
echo "..... Install cpanm ........................................"
sudo apt-get -y install cpanminus

# Install Perl dependencies from CPAN
echo "..... Install Perl dependencies from CPAN .................."
sudo cpanm Inline::CPP
sudo cpanm LWP::UserAgent
sudo cpanm Bio::DB::Fasta
sudo cpanm Log::Message::Simple


##### Instal local RNAstemloop
echo "..... Install local RNAstemloop ............................"
RNAstemloop=$ROOT_PATH"/cgi-bin/programs/RNAstemloop"
ln -s $RNAstemloop"-"$ARCH $RNAstemloop


##### Create symbolic links for programs
echo "..... Create symbolic links for programs ..................."
ln -s "/opt/VARNAv3-91.jar" $ROOT_PATH"/cgi-bin/programs/VARNA.jar"
ln -s "/opt/miRdup" $ROOT_PATH"/cgi-bin/programs/miRdup-1.4"
ln -s "/opt/tRNAscan-SE" $ROOT_PATH"/cgi-bin/programs/tRNAscan-SE"
ln -s "/opt/RNAmmer" $ROOT_PATH"/cgi-bin/programs/rnammer"


##### Deploy miRkwood data
echo "..... Deploy miRkwood data ................................."
sh $ROOT_PATH/cgi-bin/install-data.sh $ROOT_PATH/cgi-bin/data



########## Requirements to run miRkwood web-service locally ##############################

if [ ! "$CLI_ONLY" ]
then
    echo
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "Requirements to run miRkwood web-service locally"


    ##### Create directories
    echo "..... Create directories ..................................."

    # Create cgi-bin directory
    sudo mkdir -p /bio1/www/cgi-bin/

    # Link it to vagrant cgi-bin directory
    sudo ln -s /vagrant/cgi-bin /bio1/www/cgi-bin/mirkwood

    # Create html directory
    sudo mkdir -p /bio1/www/html/

    # Link it to vagrant html directory
    sudo ln -s /vagrant/html /bio1/www/html/mirkwood

    # Create results directory
    mkdir /home/vagrant/results
    sudo chmod -R 777 ~/results

    # Link results directory
    sudo ln -s /home/vagrant/results /bio1/www/html/mirkwood/results/


    ##### Install Apache (only needed for the Web service)
    echo "..... Install Apache ......................................."
    sudo apt-get -y install apache2


    ##### Install PHP (only needed for the Web service)
    echo "..... Install PHP .........................................."
    sudo apt-get -y install libapache2-mod-php5


    ##### BioInfo virtual host
    echo "..... BioInfo virtual host ................................."
    sudo cp $ROOT_PATH/provisioning/roles/bioinfo/files/bioinfo.conf /etc/apache2/sites-available/bioinfo.conf
    sudo service apache2 restart


    ##### Enable BioInfo virtualhost
    echo "..... Enable BioInfo virtualhost ..........................."
    sudo ln -s /etc/apache2/sites-available/bioinfo.conf /etc/apache2/sites-enabled/bioinfo
    sudo service apache2 restart


    ##### Disable default Apache virtualhost
    echo "..... Disable default Apache virtualhost ..................."
    sudo rm -Rf /etc/apache2/sites-enabled/000-default
    sudo service apache2 restart


    ##### Make style directory
    echo "..... Make style directory ................................."
    sudo mkdir /bio1/www/html/Style


    ##### Get BioInfo CSS
    echo "..... Get BioInfo CSS ......................................"
    sudo wget --directory-prefix=/bio1/www/html/Style/css/ http://bioinfo.lifl.fr/Style/css/bioinfo.css


    ##### Make style sub directory
    echo "..... Make style sub directory ............................."
    sudo mkdir /bio1/www/html/Style/css
    sudo mkdir /bio1/www/html/Style/css/theme


    ##### Get BioInfo CSS
    echo "..... Get BioInfo CSS ......................................"
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/box_homepage.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/divers.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/footer.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/header.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/menu_central.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/page_theme.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css http://bioinfo.lifl.fr/Style/css/table.css
    sudo wget --directory-prefix=/bio1/www/html/Style/css/theme http://bioinfo.lifl.fr/Style/css/theme/rna.css


    ##### Get BioInfo main page
    echo "..... Get BioInfo main page ................................"
    sudo wget --directory-prefix=/bio1/www/html/ http://bioinfo.lifl.fr/index.php

fi

echo
echo "miRkwood has been successfully installed."

