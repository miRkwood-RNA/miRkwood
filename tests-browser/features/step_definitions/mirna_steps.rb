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
  on(WaitingPage).loaded?
end

Then(/^I get the results page$/) do
  on(ResultsPage).wait_until do
     on(ResultsPage).results?
  end
end

