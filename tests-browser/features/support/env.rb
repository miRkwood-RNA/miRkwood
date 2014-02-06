require "watir-webdriver"
require "rspec/expectations"
require 'page-object'
require 'page-object/page_factory'

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

World(PageObject::PageFactory)

Before do
  @browser = browser
end

at_exit do
  browser.close
end
