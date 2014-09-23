Browser tests
=============

Dependencies
------------

To run these tests you will need Ruby and RubyGems.

To install the dependencies:
 
    gem install bundler
    bundle install


Running the tests
-----------------

Just run cucumber:

    cucumber

Or through `bundle`:

    bundle exec cucumber

To run tests with a specific browser, use the `BROWSER` variable:

    bundle exec cucumber BROWSER=firefox
    bundle exec cucumber BROWSER=chrome
    bundle exec cucumber BROWSER=phantomjs

To run tests in Headless mode, use the `HEADLESS` variable:

    bundle exec cucumber HEADLESS=1

To change the miRkwood instance, use the `MIRKWOOD_URL` variable:

    bundle exec cucumber MIRKWOOD_URL=http://127.0.0.1/cgi-bin/mirkwood/web_scripts/

or you may use the helper script `run-browser-tests.sh`.


Helper script examples
----------------------

The helper script `run-browser-tests.sh` takes the target environment to test as argument.

Tu run the tests locally:

    ./run-browser-tests.sh LOCAL

To run the tests against the Vagrant virtual machine:

    ./run-browser-tests.sh VAGRANT

To run the tests against BioinfoTest:

    ./run-browser-tests.sh BIOINFOTEST

To run the tests against Bioinfo:

    ./run-browser-tests.sh BIOINFO

You may still use the `BROWSER` variable:

    BROWSER=chrome ./run-browser-tests.sh BIOINFOTEST


Continuous integration
----------------------

To run the tests CI-style, ie to run in Headless and output results as HMTL and Junit,
set the `CI_MODE` variable before using `run-browser-tests.sh`

    ./run-browser-tests.sh BIOINFO

