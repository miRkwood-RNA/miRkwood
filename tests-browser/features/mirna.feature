Feature: MiRNA

    Scenario: Example sequence filling
        Given I am on MiRNA interface page
        When I use the Example feature
        Then a sequence gets filled

    Scenario: Warning if no sequence
        Given I am on MiRNA interface page
        When I launch the pipeline
        Then a no sequence warning is provided

    Scenario: Results on example
        Given I am on MiRNA interface page
        And I use the Example feature
        When I launch the pipeline
        Then I get the waiting page
        And I get the results page
