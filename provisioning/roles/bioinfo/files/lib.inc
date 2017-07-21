<?php

/* pour afficher les couleurs php en html */
function fromRGB($R, $G, $B) {
  $R=dechex($R);
  if (strlen($R)<2)
    $R='0'.$R;
  $G=dechex($G);
  if (strlen($G)<2)
    $G='0'.$G;
  $B=dechex($B);
  if (strlen($B)<2)
    $B='0'.$B;
  return "".$R."".$G."".$B;
}


function footer($soft, $responsable, $mail, $year){ 
  $len=(int)(strlen($mail)/4);
  $a = substr($mail, 0      ,$len);
  $b = substr($mail, $len   ,$len);
  $c = substr($mail,($len*2),$len);
  $d = substr($mail,($len*3),strlen($mail)-($len*3));
  echo "<div id='footer'>\n";
  echo "<table id=\"footerTable\"><tr><td id='footer_left'>";
  echo " For questions about <b>$soft</b> or for bug reports, please contact ";
  echo "
      	<script type=\"text/javascript\">
  	<!--
      
	    function escramble()
	    {
	    var a,b,c,d,e,f;
	    a=  '<a href=\"';
	    b=  'mai'; 
	    c= '$a'; 
	    b+= 'lto:'; 
	    d=  '$b';
	    i= '$c';
	    j= '$d';
	    e=  '\">';
	   
	    g=  '<';
	    h=  '/a>';
	    document.write(a+b+c+d+i+j+'?Subject=[$soft]'+e+'$responsable'+g+h);
	    }
	    escramble();

            //-->
	    </script>

     ";
    echo "</td><td id='footer_right'>";
    echo "$year - <a href=\"http://cristal.univ-lille.fr/bonsai\">Bonsai bioinformatics</a>";
    echo "</td></tr></table>";
    echo"\n </div> <!-- footer -->\n"    ;
}

?>



