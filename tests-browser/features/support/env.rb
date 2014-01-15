require "watir-webdriver"
require "rspec/expectations"

browser = nil
Browser = Watir::Browser

if ENV['FIREFOX']
  browser = Browser.new :ff
elsif ENV['CHROME']
  Selenium::WebDriver::Chrome.path = '/usr/bin/chromium-browser'
  browser = Browser.new :chrome
elsif ENV['PHANTOM']
  browser = Watir::Browser.new :phantomjs
else
  browser = Browser.new
end

Before do
  @browser = browser
end

at_exit do
  browser.close
end
