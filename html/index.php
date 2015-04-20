<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='../Style/bioinfo.css' />
        <link type='text/css' rel='stylesheet' href='./style/script.css' />
        <link type='text/css' rel='stylesheet' href='./style/rna.css' />
        <script type='text/javascript' src='/js/miARN.js'></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
<body>
    <div class="theme-border"></div>
    <a href="/">
        <div class="logo"></div>
    </a>
    <? include("./static/bioinfo_menu.txt") ?>
    <div class="bloc_droit">
        <? include("./static/header_menu.txt") ?>
        <div class="main">

<br />
	  
            <p>miRkwood is a web application that allows for the fast and easy identification of microRNAs. It is specifically designed for plant microRNAs. It offers an user-friendly interface to navigate in the data, as well as many export options to allow the user to conduct further analyses on a local computer.</p>

            <p>miRkwood comes in two versions: </p>


            <ul>
                <li><b>mirkwood ab initio</b>, for the analysis of raw genomic sequence,</li>

                <li><b>mirkwood small RNA-seq</b>, for the analysis of deep sequencing data.</li>
            </ul>

            <br /><br />

	    <center>
            <div>
                <a href="./abinitio/index.php"><span class="urlbox">Use miRkwood <i>ab initio</i></span></a>
                <a href="./smallRNAseq/index.php"><span class="urlbox">Use miRkwood small RNA-seq</span></a>
            </div>
            </center>
            <br /><br /><br />

	    <h2> Where does the name come from ? </h2> 

	    <ul>
	      <li> <a href="name.php">visit Tolkien's legendarium</a></li> 
           </ul>
	    
	    
            <h2> Authors</h2>
            <p>Sylvain Legrand (Evo-Eco-Paleo, UMR CNRS 8198 University of Lille)<br/>
            Isabelle Guigon, Jean-Fr&eacute;d&eacute;ric Berthelot, Mohcen Benmounah, H&eacute;l&egrave;ne Touzet, <a href='http://www.lifl.fr/bonsai'>Bonsai</a> (CRIStAL, UMR CNRS 9189 University of Lille,  and Inria Lille Nord Europe)</p>


          


            <h2>Contact</h2>
            <p><script type="text/javascript">
                <!--
                function escramble()
                {
                var a,b,c,d,e,f;
                a=  '<a href="';
                b=  'mai';
                c= 'mirkw';
                b+= 'lto:';
                d=  'ood@u';
                i= 'niv-l';
                j= 'ille1.fr';
                e=  '">';
                g=  '<';
                h=  '/a>';
                document.write(a+b+c+d+i+j+'?Subject=[miRkwood]'+e+c+d+i+j+g+h);
                }
                escramble();

                //-->
            </script></p>

	    <br />
        <p>This project is funded by <a href="https://www.france-genomique.org/">France Genomique.</a></p>
        </div><!-- main-->
    </div><!-- bloc droit-->
    <? include("./static/footer.txt") ?>
</body>

</html>
