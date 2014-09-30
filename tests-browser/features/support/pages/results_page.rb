class ResultsPage
  include PageObject

  @params = { :run_id => '' }
  include URL
  page_url "#{URL.url()}resultsWithID.pl?run_id=<%=params[:run_id]%>"

  div('results', :id => "table")
  expected_element(:results, 60)

  p('job_id', :id => 'job_id')
  p('precursors_count', :id => 'precursors_count')

  a('select_all', :id => 'select-all')
  a('deselect_all', :id => 'deselect-all')

  radio_button('fasta', :id => 'export-fas')

  button('export', :id => 'export-button')

end
