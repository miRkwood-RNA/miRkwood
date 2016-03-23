module URL
  def self.url()
    if ENV["MIRKWOOD_URL"]
      mirkwood_url = ENV["MIRKWOOD_URL"]
    else
      mirkwood_url = 'http://bioinfo.lifl.fr/cgi-bin/mirkwood/web_scripts/'
    end
    mirkwood_url
  end

  def self.home_url()
    if ENV["MIRKWOOD_HOME_URL"]
      mirkwood_home_url = ENV["MIRKWOOD_HOME_URL"]
    else
      mirkwood_home_url = 'http://bioinfo.lifl.fr/mirkwood/'
    end
    mirkwood_home_url
  end

  def self.ab_initio_home_url()
    if ENV["MIRKWOOD_HOME_URL"]
      mirkwood_home_url = ENV["MIRKWOOD_HOME_URL"]
    else
      mirkwood_home_url = 'http://bioinfo.lifl.fr/mirkwood/abinitio/'
    end
    mirkwood_home_url
  end

end
