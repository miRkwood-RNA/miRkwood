class InterfacePage
  include PageObject

  page_url 'http://bioinfotest.lifl.fr/cgi-bin/MiRNA/web_scripts/interface.pl'

  button('example_button', :id => "seq_button")
  button('run_button', :id => "upload")
  textarea('sequence_area', :id => "seqArea")

  def loaded?
    run_button?
  end
end


