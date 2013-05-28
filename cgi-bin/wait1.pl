#!/usr/bin/perl -w
use Class::Struct;
use CGI; 




print <<DATA
Content-type: text/html


<html>
<head>

    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
    <meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
    <link title="test" type="text/css" rel="stylesheet" href="../css/script.css" />
    <title>carnac: job running...</title>
  </head>

  <body>
    <div class="logo"></div>
    <div class="theme-border"></div>
    <!-- debut du fichier qui represente le menu colore vertical de gauche -->

<div class="bloc_gauche">
  <div class="nav-top"></div>
  <div class="menu">
    
    <!-- SEQUOIA -->
    <ul class="menu_partie">
      <li class ="case_titre Research_titre"> 
	<a href="/index.php">bioinfo.lifl.fr</a>
      </li>
      <li class ="case Research"> 
	<a href="http://www.lifl.fr/bonsai">Bonsai</a>
      </li>
    </ul>
    
    <!-- Mreps -->
    <ul class="menu_partie"> 
      <li class ="case_titre Mreps_titre">
	<a href="/mreps/index.php">mreps</a>
      </li>
    </ul>
    
    <!-- Yass -->
    <ul class="menu_partie">
      <li class ="case_titre Yass_titre">
	<a href="/yass/index.php">YASS</a>
      </li>
    </ul> 

    <!-- Magnolia -->
    <ul class="menu_partie">
      <li class="case_titre Magnolia_titre">
	<a href="/magnolia/index.php">Magnolia</a>
      </li>
    </ul>
    

    <ul class="menu_partie">
      <li class="case_titre proteins_titre">
	<a href="/proteins/index.php">Proteins</a>
      </li>
    <!-- Path -->
      <li class ="case proteins">
        <a href="/path/index.php">Path</a>
      </li>
    <!-- Protea -->
      <li class="case proteins">
	<a href="/protea/index.php">Protea</a>
      </li>    
    <!-- Reblosum -->
      <li class="case proteins">
	<a href="/reblosum/index.php">ReBLOSUM</a>
      </li>
    </ul>
        
    <!-- RNA -->
    <ul class="menu_partie"> 
      <li class="case_titre Rna_titre">
	<a href="/RNA/index.php">RNA</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/carnac/index.php">Carnac</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/RNAfamily/rnafamily.php">RNAfamily</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/gardenia/index.php">Gardenia</a>
      </li>
      <li class="case Rna">
	<a href="/RNA/regliss/index.php">Regliss</a>
      </li>
      <li class="case Rna">
	<a href="/CGseq/index.php">CG-seq</a>
      </li>

    </ul>
    
    <!-- TFM -->
    <ul class="menu_partie"> 
      <li class="case_titre Tfm_titre">
	<a href="/TFM/">TFM</a>
      </li>
      <li class ="case Tfm">
	<a href="/TFM/TFME/">TFM-Explorer</a>
      </li>
      <li class ="case Tfm">
	<a href="/TFM/TFMscan/">TFM-Scan</a>
      </li>
      <li class="case Tfm">
	<a href="/TFM/TFMpvalue/">TFM-Pvalue</a>
      </li>      
      <li class="case Tfm">
	<a href="/TFM/TFM-CUDA/">TFM-CUDA</a>
      </li>      
    </ul>

    <!-- NRPS -->
    <ul class="menu_partie no-margin">
      <li class ="case_titre Nrps_titre">
	<a href="/norine/">NRPS</a>
      </li>
      <li class ="case Nrps"> 
	<a href="/norine/form.jsp">Norine</a>
      </li>

    </ul>
    
  </div> <!-- left menu -->
  <div class="nav-bottom"></div>  
</div> <!-- bloc gauche-->


<script type="text/javascript">
var pkBaseURL = (("https:" == document.location.protocol) ? "https://bioinfo.lifl.fr/piwik/" : "http://bioinfo.lifl.fr/piwik/");
document.write(unescape("%3Cscript src='" + pkBaseURL + "piwik.js' type='text/javascript'%3E%3C/script%3E"));
</script>
<script type="text/javascript">
try {
var piwikTracker = Piwik.getTracker(pkBaseURL + "piwik.php", 1);
piwikTracker.trackPageView();
piwikTracker.enableLinkTracking();
} catch( err ) {}
</script>
<noscript><p><img src="http://bioinfo.lifl.fr/piwik/piwik.php?idsite=1" style="border:0" alt=""/></p></noscript>


<!-- End Piwik Tag -->

<!-- fin du fichier menu.txt -->




    <div class="bloc_droit">
      <div class="frametitle">
	<h1>
carnac
</h1>
	<div class="menu_central">
  <a href="/RNA/carnac/index.php">home</a>
  <a href="/cgi-bin/RNA/carnac/carnac.py">web server</a>
  <a href="/RNA/carnac/help.php">help</a>
  <a href="/RNA/carnac/examples.php">examples</a>
  <a href="/RNA/carnac/id.php">retrieve result with an ID</a>
</div>

<!-- Piwik -->
<script type="text/javascript">
var pkBaseURL = (("https:" == document.location.protocol) ? "https://bioinfo.lifl.fr/piwik/" : "http://bioinfo.lifl.fr/piwik/");
document.write(unescape("%3Cscript src='" + pkBaseURL + "piwik.js' type='text/javascript'%3E%3C/script%3E"));
</script><script type="text/javascript">
try {
var piwikTracker = Piwik.getTracker(pkBaseURL + "piwik.php", 3);
piwikTracker.trackPageView();
piwikTracker.enableLinkTracking();
} catch( err ) {}
</script><noscript><p><img src="http://bioinfo.lifl.fr/piwik/piwik.php?idsite=3" style="border:0" alt=""/></p></noscript>
<!-- End Piwik Tag -->

      </div>

      <div class="main">
	

<div class="dialog">
  <br/>
  <br/>

  Your request has been successfully submitted to CARNAC. 


  <br />

  Your ID is <B>May_28_2012_19_52_03_17008</B>. <br />

  <br /><br />

  This page is updated every five seconds.

  <br /><br /><br /><br />

  <b><font size="+1">Please wait.</font></b>

  <br /><br /><br />
</div>


<script type="text/javascript"> 
<!--
window.setTimeout('window.location.replace("/cgi-bin/RNA/carnac/carnac.py?command=result&amp;run_id=May_28_2012_19_52_03_17008")', 1000*5);
--> 
</script>
	
      </div>
      <div class="footer">
	For questions about <b>carnac</b> or for bug reports, please contact 
	
<script type="text/javascript">
  
  <!--
      function escramble()
      {
      var a,b,c,d,e,f;
      a=  '<a href="';
      b=  'mai'; 
      c= 'car'; 
      b+= 'lto:'; 
      d=  'nac';
      i= '@univ-';
      j= 'lille1.fr';
      e=  '">';
      
      g=  '<';
      h=  '/a>';
      document.write(a+b+c+d+i+j+'?Subject=[carnac:Suggestion/Comment/Critics]'+e+'carnac'+g+h);
      }
      escramble();
      //-->
  
</script>

      </div>
      <div class="copyright">
	2009 - <a href="http://www.lifl.fr/SEQUOIA">SEQUOIA </a>
      </div>
      
    </div>

  </body>
</html>






DATA
