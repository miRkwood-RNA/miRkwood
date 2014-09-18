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

    cucumber BROWSER=firefox
    cucumber BROWSER=chrome
    cucumber BROWSER=phantomjs

To change the miRkwood instance, use the `MIRKWOOD_URL` variable:

    cucumber MIRKWOOD_URL=http://127.0.0.1/cgi-bin/mirkwood/web_scripts/


Tu run the tests locally:

    bundle exec cucumber BROWSER=chrome HEADLESS=True MIRKWOOD_URL=http://127.0.0.1/cgi-bin/mirkwood/web_scripts/


To run tests CI-style:

    bundle exec cucumber --format pretty --format junit --out . --format html --out browsertest-report.html

