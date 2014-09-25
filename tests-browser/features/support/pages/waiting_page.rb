class WaitingPage
  include PageObject

  div('waiting', :class => "waitMessage")

  expected_element(:waiting)
end

