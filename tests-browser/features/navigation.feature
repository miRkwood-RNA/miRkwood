Feature: Navigation on the miRkwood home page

    Background:
        Given I am on miRkwood ab initio home page

    Scenario: Going to web server page from menu
        When I select web server in the menu
        Then I should land on miRkwood interface page

    Scenario: Going to web server page from link
        When I select web server in the text
        Then I should land on miRkwood interface page

    Scenario: Going to help page from menu
        When I select help in the menu
        Then I should land on miRkwood help page

    Scenario: Going to ID page from menu
        When I select retrieve result in the menu
        Then I should land on miRkwood ID page

