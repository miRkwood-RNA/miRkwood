Feature: Navigation on the miRkwood home page

    Scenario: Going to web server page from menu
        Given I am on miRkwood home page
        When I select web server in the menu
        Then I should land on miRkwood interface page

    Scenario: Going to web server page from link
        Given I am on miRkwood home page
        When I select web server in the text
        Then I should land on miRkwood interface page

    Scenario: Going to help page from menu
        Given I am on miRkwood home page
        When I select help in the menu
        Then I should land on miRkwood help page

    Scenario: Going to ID page from menu
        Given I am on miRkwood home page
        When I select retrieve result in the menu
        Then I should land on miRkwood ID page

