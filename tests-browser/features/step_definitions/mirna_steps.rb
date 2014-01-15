Given(/^I am on MiRNA interface page$/) do
  @browser.goto("http://bioinfotest.lifl.fr/")
  @browser.a(:text => "MiRNA").click
end

When(/^I use the Example feature$/) do
  @browser.button(:id => "seq_button").click
end

When(/^I launch the pipeline$/) do
  @browser.button(:id => "upload").click
end

Then(/^a sequence gets filled$/) do
  @browser.textarea(:id => "seqArea").value.should include("contig15750") 
end

Then(/^a no sequence warning is provided$/) do
  @browser.alert.should exist
  @browser.alert.text.should include("You must provide sequences")
  @browser.alert.ok
end

Then(/^I get the waiting page$/) do
  @browser.div(:class => "waitMessage").should exist
end

Then(/^I get the results page$/) do
  @browser.div(:id => "table").wait_until_present
  @browser.div(:id => "table").should exist
end

