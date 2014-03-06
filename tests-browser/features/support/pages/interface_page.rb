class InterfacePage
  include PageObject

  include URL
  page_url URL.url() + 'interface.pl'

  a('example_button', :id => "seq_button")
  a('area_clear', :id => "area_clear")
  button('run_button', :id => "upload")
  textarea('sequence_area', :id => "seqArea")

  def loaded?
    run_button?
  end
end


