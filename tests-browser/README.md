Browser tests
=============

Dependencies
------------

To run these tests you will need Ruby and RubyGems.

(Research whether it is better to use RVM or APT)


To install the dependencies:
 
    gem install bundler
    bundle install


Running the tests
-----------------

Just run cucumber:

    cucumber

Or through `bundle`:

    bundle exec cucumber

To run tests with a specific browser:

    cucumber BROWSER=firefox
    cucumber BROWSER=chrome
    cucumber BROWSER=phantomjs

To run tests CI-style:

    bundle exec cucumber --format pretty --format junit --out . --format html --out browsertest-report.html

