Feature: MiRNA

    Scenario: Example sequence filling
        Given I am on MiRNA interface page
        When I use the Example feature
        Then a sequence gets filled

    Scenario: Clearing area
        Given I am on MiRNA interface page
        And I use the Example feature
        When I use the Clear feature
        Then the sequence area is clear

    Scenario: Warning if no sequence
        Given I am on MiRNA interface page
        Then a no sequence warning is provided when I launch the pipeline

    Scenario: Results on example
        Given I am on MiRNA interface page
        And I use the Example feature
        When I launch the pipeline
        Then I get the waiting page
        And I get the results page
