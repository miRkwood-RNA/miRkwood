<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <link type='text/css' rel='stylesheet' href='./style/mirkwood.css' />
        <link type='text/css' rel='stylesheet' href='../Style/css/bioinfo.css' />
        <script type="text/javascript" src="../libs/jquery-1.11.3.min.js"></script>
        <script type="text/javascript" src="./js/header.js"></script>
        <title>miRkwood - MicroRNA identification</title>
    </head>
    <body>
        <div class="frametitle">
            <h1 id="title" onclick="loadLink('/mirkwood/index.php');">miRkwood</h1>                 
        </div>
        
        <div class="tabs" id="menu_central" style="display:inline-block"><?php include("./static/header_menu.txt")?></div>
        <div id="arborescence"></div>

        <div id="main">

            <br />
<h2>Fast and easy identification of plant microRNAs</h2>

 <p>miRkwood is a software package for the discovery of microRNAs and their hairpin precursors in plant genomes. It combines multiple evidences to support the prediction: thermodynamical stability, conservation, miRNA:miRNA* duplex,... mirkwood has a user-friendly interface to navigate in the data, as
well as many export options to conduct further
analyses. </p>

<p>It comes in two versions.</p>

<br />
       
<p>
<div class="mirkwood_button">
 
     <a href="./abinitio/index.php"><b>Run mirkwood ab initio</b>, for the analysis of  raw
  genomic sequences</a>

</div>
</p>
<p>
<div class="mirkwood_button">
  
     <a href="./smallRNAseq/index.php"><b>Run mirkwood small RNA-seq</b>, for the analysis of deep
  sequencing data</a>

</div>
</p>

            
            <br />

            <h3> Where does the name come from ? </h3> 

            <p> Visit <a href="name.php">Tolkien's legendarium</a>. </p>
            
        
            <h3> Authors</h3>
            <p>Sylvain Legrand (Evo-Eco-Paleo, UMR CNRS 8198 University of Lille)<br/>
            Isabelle Guigon, Jean-Fr&eacute;d&eacute;ric Berthelot, Mohcen Benmounah, H&eacute;l&egrave;ne Touzet, <a href='http://www.lifl.fr/bonsai'>Bonsai</a> (CRIStAL, UMR CNRS 9189 University of Lille,  and Inria Lille Nord Europe)
            </p>

            <h3>Contact</h3>
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

           
<h3> Funding</h3>
            <p>This project was supported  by <a href="https://www.france-genomique.org/">France Genomique</a>.</p>

        </div><!-- main-->
        
    <?php require("/bio1/www/html/lib.inc")?>
    <?php footer("miRkwood","miRkwood", "mirkwood@univ-lille1.fr","2013"); ?>

    </body>

</html>
