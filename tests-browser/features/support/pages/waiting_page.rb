class WaitingPage
  include PageObject

  div('waiting', :class => "waitMessage")

  def loaded?
    waiting?
   end
end

