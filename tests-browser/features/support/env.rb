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

After do |scenario|
  FileUtils.mkpath('html-report/screenshots') if not File.directory?('html-report/screenshots')
  filename = "screenshots/Screenshot_#{scenario.name.gsub(' ','_').gsub(/[^0-9A-Za-z_]/, '')}.png"
  browser.screenshot.save "./html-report/#{filename}"
  embed "./#{filename}", 'image/png'
end

at_exit do
  browser.close
end
