When(/^I select help in the menu$/) do
  on(HomePage).menu_help
end

When(/^I select retrieve result in the menu$/) do
  on(HomePage).menu_retrieve_result
end

When(/^I select web server in the menu$/) do
  on(HomePage).menu_web_server
end


When(/^I select web server in the text$/) do
  on(HomePage).link_web_server
end

