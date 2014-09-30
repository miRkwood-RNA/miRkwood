When(/^I use the Example feature$/) do
  on(InterfacePage).example_button
end

When(/^I launch the pipeline$/) do
  on(InterfacePage).run_button
end

Then(/^a sequence gets filled$/) do
  expect(on(InterfacePage).sequence_area).to include("sample")
end

When(/^I use the Clear feature$/) do
  on(InterfacePage).area_clear
end

Then(/^the sequence area is clear$/) do
  expect(on(InterfacePage).sequence_area).to eq('')
end

Then(/^a no sequence warning is provided when I launch the pipeline$/) do
  on(InterfacePage) do |page|
    message = page.alert do
      page.run_button
    end
  expect(message).to eq("You must provide sequences")
  end
end

