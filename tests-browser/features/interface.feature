Feature: miRkwood interface page

    Scenario: Example sequence filling
        Given I am on miRkwood interface page
        When I use the Example feature
        Then a sequence gets filled

    Scenario: Clearing area
        Given I am on miRkwood interface page
        And I use the Example feature
        When I use the Clear feature
        Then the sequence area is clear

    Scenario: Warning if no sequence
        Given I am on miRkwood interface page
        Then a no sequence warning is provided when I launch the pipeline

    @execution
    Scenario: Results on example
        Given I am on miRkwood interface page
        And I use the Example feature
        When I launch the pipeline
        Then I should land on miRkwood waiting page
        And I should land on miRkwood results page
