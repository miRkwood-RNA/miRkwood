class MirkwoodPage
  include PageObject

  a('menu_home', :text => 'home')
  a('menu_web_server', :text => 'web server')
  a('menu_help', :text => 'help')
  a('menu_retrieve_result', :text => 'retrieve result with an ID')

end

class HomePage < MirkwoodPage
  include URL
  page_url URL.home_url()

  a('link_help', :text => 'help')
  a('link_web_server', :text => 'use web server')

  image('mirkwood_diagram', :alt => 'miRkwood diagram')
  expected_element(:mirkwood_diagram)
end

class HelpPage < MirkwoodPage
  include URL
  page_url URL.home_url() + 'help.php'

  div('toc', :class=> 'table-of-contents')
  a('link_web_server', :text => 'miRkwood website')

  expected_element(:toc)
  expected_title "miRkwood - MicroRNA identification - Help"

end

class IDPage < MirkwoodPage
  include URL
  page_url URL.home_url() + 'id.php'
  text_field('id_input', :id=> 'run_id')

  expected_element(:id_input)
end

