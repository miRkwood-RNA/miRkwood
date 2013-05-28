function readFile(url)  // fonction ajax qui permet de lire un fichier txt 
{
     if(window.XMLHttpRequest) // FIREFOX
          xhr_object = new XMLHttpRequest(); 
     else if(window.ActiveXObject) // IE
          xhr_object = new ActiveXObject("Microsoft.XMLHTTP"); 
     else 
          return(false); 
     xhr_object.open("GET", url, false); 
     xhr_object.send(null); 
     if(xhr_object.readyState == 4){ var f = xhr_object.responseText; document.getElementById('alignement').innerText = f;return f;}
     else{ return(false);}
}

function displayFile(url)
{console.log("one")
		var file = readFile(url);
}
