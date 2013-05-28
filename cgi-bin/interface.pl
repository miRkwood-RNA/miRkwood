#!/usr/bin/perl -w




print <<DATA
Content-type: text/html


<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
<meta name="keywords" content="RNA, ARN, mfold, fold, structure, prediction, secondary structure" />
<link title="test" type="text/css" rel="stylesheet" href="/arn/css/script.css" />

        <script src="/arn/js/miARN.js" type="text/javascript" LANGUAGE="JavaScript"></script>
       
        <title>miREST :: identification of miRNA/miRNA hairpins in plants</title>


</head>
<body>
<div class="theme-border"></div>
<div class="logo"></div>
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



<!-- End Piwik Tag -->

<!-- fin du fichier menu.txt -->


<div class="bloc_droit">
<div class="frametitle">
 <h1>miREST :: identification of miRNA/miRNA hairpins in plants</h1>
   <div class="menu_central">
     <a href="./interface.pl">home</a> 
	<a href="./help.pl">help</a>

     <a href="./id.pl">retrieve result with an ID</a>
  	
     </div>
</div>

<div class="main">

        <form  name="form" onsubmit="verifySequence()" onsubmit="wainting()" method="post" action="./results.pl" enctype="multipart/form-data">
            <fieldset id="fieldset">    
				 <div class="forms">
				<tr>
					<td class="label"> 
						Enter a name for your job (optional)</i>: 
						<input type="text" name="job" size="20">
					</td>
				  </tr>
				</div>
				<div class="forms">
					<label class="textDiv" for="seqArea"><legend><h2>Paste your sequence(s)[<a href="./help.pl">?</a>]:</h2></legend> </label><br />
					<textarea name="seqArea"  rows="10" cols="150" ></textarea>            
					<label class="textDiv" style="color:#333" for="seq">Or, upload file :</label><br />
					<input type="file" name="seqFile" id="file" />
				
				</div>
				  <div class="forms">
					 <label class="checkbox" for="db"><h2>Mask coding regions [<a href="./help.pl">?</a>] <i><small>(may be slow) </small></i> :
					<input  id ="CDS" type="checkbox" name="check" value="checked" onclick="showHideBlock()"> </h2>               
					<div id="menuDb"> 
					<label class="choixDiv" for="db">Choose organism database :</label>
						<select class="db" name="db">
							<option class="db" selected>ATpepTAIR10</option>
							<option class="db">plante</option>
						</select>
					</div>
               </div>
               <div class="forms">
				<h2>Select additional features:</h2>
				<P>
                <input class="checkbox" type="checkbox" checked="checked" name="randfold" value="randfoldChecked">Compute thermodynamic stability [<a href="./help.pl">?</a>]</input>
                    <br />
                <input class="checkbox" type="checkbox" checked="checked" name="mfei" value="mfeiChecked">Compute MFE/MFEI/AMFE (minimal folding energy)[<a href="./help.pl">?</a>]</input>
                  <br />
                <input class="checkbox" type="checkbox" checked="checked" name="align" value="alignChecked">Align against mature microRNAs miRBase [<a href="./help.pl">?</a>]</input>
                   <br />
                <input class="checkbox" type="checkbox" checked="checked" name="selfContain" value="SCChecked">Compute self_containment index [<a href="./help.pl">?</a>]</input>
                </div>
               </P>
                
                </label>
				<div class="forms">
				<tr>
					<td class="label"> 
						Enter your <b>E-mail</b> address <i>(optional)</i>: 
						<input type="text" name="mail" size="20">
					</td>
				  </tr>
				</div>
              <input style="margin-left:632px" type="submit" name="upload" id="upload" value="Submit JOB">
                
            
            
            </fieldset>
        
        </form>
        
        
        
        
        
        
</div><!-- main -->

 <div class="copyright">
          2004 - <a href="http://www.lifl.fr/bonsai">Bonsai bioinformatics</a>
         </div>
</div><!-- bloc droit-->

</body>
</html>

DATA
###End###
