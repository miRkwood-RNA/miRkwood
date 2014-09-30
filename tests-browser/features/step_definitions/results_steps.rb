Given(/^I am on miRkwood results page with ID (.*)$/) do |id|
  visit_page(ResultsPage, :using_params => {:run_id => id})
end

Then(/^the job ID is (.*)$/) do |id|
  text = "Job ID: #{id}"
  expect(on(ResultsPage).job_id).to eq(text)
end

Then(/^there are (\d+) candidates found$/) do |count|
  text = "#{count} miRNA precursor(s) found"
  expect(on(ResultsPage).precursors_count).to eq(text)
end

When(/^I select all the candidates in the form$/) do
  on(ResultsPage).select_all
end

When(/^I select the Fasta export$/) do
  on(ResultsPage).select_fasta
end

When(/^I export the results$/) do
  on(ResultsPage).export
end

