Deployment of miRkwood
======================

miRkwood is to be run on the BioInfo server. This procedure is tailored to that usecase.

It was decided that miRkwood should be deployed there using Subversion.
This is well adapted as we are using an interpreted language which does not need compilation.


On a regular server
-------------------

Deployment is made through Subversion: the two directories, cgi-bin and html,
are checked-out from the Inria forge to the right place in the filesystem:

- /bio1/www/cgi-bin/mirkwood --> cgi-bin directory
- /bio1/www/html/mirkwood --> html directory

The Shell script `deploy-mirkwood-locally` can be used for that.
It takes as optional argument the number of the Subversion revision to deploy.

This is done on the test and qualification servers.

To deploy the latest versuib on a distant server, just run:

    ssh <server> < deploy-mirkwood-locally

On BioInfoDev, you may use the wrapper script `deploy-mirkwood-to-bioinfodev`,
which also takes as optional argument the number of the Subversion revision to deploy.

    ./deploy-mirkwood-to-bioinfodev

    ./deploy-mirkwood-to-bioinfodev 11460



On BioInfo
----------

BioInfo does not accept connections to the outside world, so we cannot check out the Forge repository.

The solution is to mirror the contents of the qualification server using scp.

The Shell script `deploy-mirkwood-to-bioinfo` can be used for that.

You must use it from bioinfodev, or use the wrapper `mirror-deploy-mirkwood-from-local-to-bioinfo`.
