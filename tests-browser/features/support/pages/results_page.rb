class ResultsPage
  include PageObject

  div('results', :id => "table")

  expected_element(:results, 60)

end
