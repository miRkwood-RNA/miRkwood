require "watir-webdriver"
require "rspec/expectations"
require 'page-object'
require 'page-object/page_factory'
require 'fileutils'

browser = nil
Browser = Watir::Browser

if ENV['HEADLESS']
  require 'headless'
  headless = Headless.new
  headless.start
  at_exit do
    headless.destroy
  end
end

if ENV['BROWSER'] == 'firefox'
  browser = Browser.new :ff
elsif ENV['BROWSER'] == 'chrome'
  Selenium::WebDriver::Chrome.path = '/usr/bin/chromium-browser'
  browser = Browser.new :chrome
elsif ENV['BROWSER'] == 'phantomjs'
  browser = Watir::Browser.new :phantomjs
else
  browser = Browser.new
end

World(PageObject::PageFactory)

Before do
  @browser = browser
end

After('@screenshot') do |scenario|
  embed_screenshot()
end

After do |scenario|
  if(scenario.failed?)
    embed_screenshot()
  end
end

def embed_screenshot
  encoded_img = @browser.driver.screenshot_as(:base64)
  embed("data:image/png;base64,#{encoded_img}", 'image/png')
end

at_exit do
  browser.close
end
