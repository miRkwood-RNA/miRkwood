module URL
  def self.url()
    if ENV["MIRKWOOD_URL"]
      bioinfo_url = ENV["MIRKWOOD_URL"]
    else
      bioinfo_url = 'http://bioinfotest.lifl.fr/cgi-bin/mirkwood/web_scripts/'
    end
    bioinfo_url
  end

  def self.home_url()
    if ENV["MIRKWOOD_HOME_URL"]
      mirkwood_home_url = ENV["MIRKWOOD_HOME_URL"]
    else
      mirkwood_home_url = 'http://bioinfotest.lifl.fr/mirkwood/'
    end
    mirkwood_home_url
  end
end
