module URL
  def self.url()
    if ENV["MIRKWOOD_URL"]
      bioinfo_url = ENV["MIRKWOOD_URL"]
    else
      bioinfo_url = 'http://bioinfotest.lifl.fr/cgi-bin/mirkwood/web_scripts/'
    end
    bioinfo_url
  end
end
