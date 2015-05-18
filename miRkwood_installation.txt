miRkwood
========

miRkwood is a computational pipeline for the identification of miRNAs and their hairpin precursors.

It is constituted of:
- a back-end in Perl, running many third-party programs, stored in `cgi-bin`
- a front-end in Perl/JavaScript, stored in `html`


miRkwood comes with a fair number of dependencies. It was tested on Ubuntu 12.04.
The easiest way to deploy miRkwood is to create a virtual machine (VM) and to use
the configuration management software Ansible.


1. Recommended installation

1.1. Create the VM

- Install VirtualBox <https://www.virtualbox.org/wiki/Downloads>

- Install Vagrant in its most recent version : <http://www.vagrantup.com/downloads.html>
  (tested on Vagrant 1.4.3 and Vagrant 1.6.5)


1.2. Install Ansible

Ansible is an IT automation tool. It can configure systems, deploy software, and orchestrate 
more advanced IT tasks such as continuous deployments or zero downtime rolling updates.

- Install Ansible in its most recent version (at least 1.6) <http://docs.ansible.com/intro_installation.html>
(tested with Ansible 1.6 and 1.7)


1.3. Install miRkwood dependencies

- Extract miRkwood archive and cd to the newly created directory.

- Run Vagrant at the miRkwood repository root
    `vagrant up`

- Vagrant will
    - download an Ubuntu ISO
    - mount the local development directory in the virtual machine
    - provision it using Ansible
    - set up network forwarding

This step can take up to 10 minutes.

Congratulations! A miRkwood instance is now running at <http://192.168.33.20/mirkwood>


2. Alternative installation

If you are unlucky enough that you cannot install Ansible (this problem may
occur on Mac systems), you will have to install yourself all dependencies.

2.1 Create the VM

We still recommand to run miRkwood on an Ubuntu 12.04.

# Need details here to install Ubuntu on the VM.


2.2 Install dependencies

All these steps take place on the VM.

Run script install-dependencies.sh in super-user mode.

