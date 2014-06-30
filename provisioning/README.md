MiRkwood/Ansible configuration
==============================

To deploy miRkwood locally
--------------------------

  ansible-playbook -i locals mirkwood-local.yml --ask-sudo-pass


To play only the after deployment parts
---------------------------------------

  ansible-playbook -i locals mirkwood-local.yml --ask-sudo-pass --tags "post-deploy"
