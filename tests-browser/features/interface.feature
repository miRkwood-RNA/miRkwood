Feature: miRkwood interface page
    Background:
        Given I am on miRkwood interface page

    @screenshot
    Scenario: Example sequence filling
        When I use the Example feature
        Then a sequence gets filled

    Scenario: Clearing area
        Given I use the Example feature
        When I use the Clear feature
        Then the sequence area is clear

    Scenario: Warning if no sequence
        Then a no sequence warning is provided when I launch the pipeline

