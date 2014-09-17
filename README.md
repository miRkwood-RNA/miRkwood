miRkwood
========

miRkwood is a computational pipeline for the identification of miRNAs and their hairpin precursors.

It is constituted of:
- a back-end in Perl, running many third-party programs, stored in `cgi-bin`
- a front-end in Perl/JavaScript, stored in `html`


Installation for usage
----------------------

The easiest way to deploy miRkwood is to use the configuration management software Ansible.


Installation for developers
---------------------------

The easiest way to deploy a miRkwood development environment is to use Vagrant,
in conjunction with Ansible.

Vagrant is a tool to create and configure virtual development environments.
It can be considered a wrapper around virtualization software such as VirtualBox
and configuration management software such as Chef, Salt and Puppet − or Ansible in our case.

Steps are:

- Install VirtualBox <https://www.virtualbox.org/wiki/Downloads>

- Install Vagrant in its most recent version : <http://www.vagrantup.com/downloads.html>
  (tested on Vagrant 1.4.3 and Vagrant 1.6.5)

- Install Ansible in its most recent version (at least 1.6) <http://docs.ansible.com/intro_installation.html>
  (tested with Ansible 1.6 and 1.7)

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

- Congratulations! A miRkwood instance is now running at <http://192.168.33.20/mirkwood>
