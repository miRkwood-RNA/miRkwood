@dummy-results-setup
Feature: Navigation on the miRkwood home page

    Scenario Outline:
        Given I am on miRkwood results page with ID <id>
        Then the job ID is <id>
        And there are <count> candidates found
    Examples:
       | id       | count |
       | TestData | 2     |

    Scenario:
        Given I am on miRkwood results page with ID TestData
        When I select all the candidates in the form
        And I select the Fasta export

