Given(/^I am on MiRNA interface page$/) do
  visit InterfacePage
end

When(/^I use the Example feature$/) do
  on(InterfacePage).example_button
end

When(/^I launch the pipeline$/) do
  on(InterfacePage).run_button
end

Then(/^a sequence gets filled$/) do
  on(InterfacePage).sequence_area.should include("sample")
end

When(/^I use the Clear feature$/) do
  on(InterfacePage).area_clear
end

Then(/^the sequence area is clear$/) do
  on(InterfacePage).sequence_area.should eq('')
end

Then(/^a no sequence warning is provided when I launch the pipeline$/) do
  on(InterfacePage) do |page|
    message = page.alert do
      page.run_button
    end
  message.should == "You must provide sequences"
  end
end

Then(/^I get the waiting page$/) do
  on(WaitingPage).has_expected_element?
end

Then(/^I get the results page$/) do
  on(ResultsPage).has_expected_element?
end

Given(/^I am on miRkwood home page$/) do
  visit HomePage
end

When(/^I select web server in the menu$/) do
  on(HomePage).menu_web_server
end

Then(/^I should land on miRkwood interface page$/) do
  on(InterfacePage).has_expected_element?
end

Then(/^I should land on miRkwood help page$/) do
  on(HelpPage).has_expected_element?
end

When(/^I select help in the menu$/) do
  on(HomePage).menu_help
end

When(/^I select web server in the text$/) do
  on(HomePage).link_web_server
end

When(/^I select retrieve result in the menu$/) do
  on(HomePage).menu_retrieve_result
end

Then(/^I should land on miRkwood ID page$/) do
  on(IDPage).has_expected_element?
end

