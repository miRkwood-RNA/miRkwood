Given(/^I am on miRkwood interface page$/) do
  visit InterfacePage
end

Given(/^I am on miRkwood home page$/) do
  visit HomePage
end


Then(/^I should land on miRkwood waiting page$/) do
  on(WaitingPage).has_expected_element?
end

Then(/I should land on miRkwood results page$/) do
  on(ResultsPage).has_expected_element?
end

Then(/^I should land on miRkwood interface page$/) do
  on(InterfacePage).has_expected_element?
end

Then(/^I should land on miRkwood ID page$/) do
  on(IDPage).has_expected_element?
end

Then(/^I should land on miRkwood help page$/) do
  on(HelpPage).has_expected_element?
end

