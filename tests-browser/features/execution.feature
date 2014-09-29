@execution
Feature: miRkwood execution

    Scenario: Results on example
        Given I am on miRkwood interface page
        And I use the Example feature
        When I launch the pipeline
        Then I should land on miRkwood waiting page
        And I should land on miRkwood results page
