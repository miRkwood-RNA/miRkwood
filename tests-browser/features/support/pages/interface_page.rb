class InterfacePage
  include PageObject

  include URL
  page_url URL.url() + 'interface.pl'

  button('example_button', :id => "seq_button")
  button('run_button', :id => "upload")
  textarea('sequence_area', :id => "seqArea")

  def loaded?
    run_button?
  end
end


