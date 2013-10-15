/**
 * Instantiation  de classe results.js
 */
function main(id)
{
	myResults = new results(id);
	rowsNumber = myResults.getSequencesNamesList().length; //récupération de la longueur du tableau à afficher(ligne)
	columnsNumber = myResults.getFactorsNamesList().length;//nom de collone
	createGrid('table',rowsNumber+1,columnsNumber+1);// fonction permettant de créer le tableau par rapport au résultats
	//#######plugin imgPreview########
	jQuery.noConflict();
	(function($){  
		
	$('tbody#cases a').imgPreview({
		containerID: 'imgPreviewWithStyles',
		imgCSS: {
		// Limit preview size:
		height: 200
		},
		// When container is shown:
		onShow: function(link){
		// Animate link:
		$(link).stop().animate({opacity:0.4});
		// Reset image:
		$('img', this).stop().css({opacity:0});
		},
		// When image has loaded:
		onLoad: function(){
		// Animate image
		$(this).animate({opacity:1}, 300);
		},
		// When container hides: 
		onHide: function(link){
		// Animate link:
		$(link).stop().animate({opacity:1});
		}
		});

	})(jQuery);

}

/**
 * Fonction permettant d'afficher la valeur d'une cellule
 */
function showCellInfo(i,j)
{
	var id_job = document.getElementById('id_job').innerHTML;
	if ((i!=0)&&(j!=0))
	{
		var factor = myResults.getFactorsNamesList()[j-1];
		if ((factor != "image") &&(factor != "alignment") &&(factor != "mfei")&&(factor != "mfe")&&(factor != "amfe")&&(factor != "position")&&(factor != "p_value"))
		{
			window.open("./resultsCell.pl?typePage=simpleCell&nameSeq="+myResults.getSequenceNameByIndex(i-1)+"&factor="+ myResults.getFactorNameByIndex(j-1)+"&value="+myResults.getValueByIndices(i-1,j-1)+"&position="+myResults.getValueByIndices(i-1,0));
			
			//window.open("./resultsCell.pl?typePage=alignement&url="+myResults.getValueByIndices(i-1,j-1));
		}
		// 
		//{ 
	//		window.open("./resultsCell.pl?typePage=simpleCell&nameSeq="+myResults.getSequenceNameByIndex(i-1)+"&factor="+ myResults.getFactorNameByIndex(j-1)+"&value="+myResults.getValueByIndices(i-1,j-1)+"&position="+myResults.getValueByIndices(i-1,0));
			//document.getElementById("singleCell").innerHTML="<div id = 'showInfo'> <h2 class='titre'><u>Sequence Informations </u></h2><br/> <li><b>Name sequence :</b> " + myResults.getSequenceNameByIndex(i-1) + "</li><br/> <li> <b>" + myResults.getFactorNameByIndex(j-1) + " </b>: "+ myResults.getValueByIndices(i-1,j-1) + " </li> <br/> </div>"
		
		//}
	} 
	if ((i!=0)&&(j==0)) 
	{
		var nameTemp = myResults.getSequenceNameByIndex(i-1);
		var positions = myResults.getValuesByFactorName("position");
		var factorsTemp = myResults.getSequenceByNameFactors(nameTemp,positions[i-1]);
		window.open("./resultsRow.pl?jobID="+id_job+"&nameSeq="+nameTemp+"&position="+factorsTemp.position);
		//document.getElementById("singleCell").innerHTML="<div id = 'showInfo'> <h2 class='titre'><u>Sequence Informations </u></h2><br/> <li><b>Name sequence :</b> " + nameTemp + "</li><br/><li> <b>MFEI : </b>" + factorsTemp.mfei + "</li> <br/> <li><b>P_value :</b> " + factorsTemp.p_value + "</li><br/> <li><b>Initial position : </b>" + factorsTemp.position +"</li> <br/> </div>"
	} 
	/**
	if ((i==0)&&(j!=0))  
	{	
		var factor = myResults.getFactorsNamesList()[j-1];
		if ((factor != "image") &&(factor != "quality") &&(factor != "alignment")&&(factor != "amfe")&&(factor != "mfe")&&(factor != "position")) // quand on clique pas sur image et alignement 
		{
			var values= "";
			var factorsTemp = myResults.getFactorNameByIndex(j-1);
			var valuesFactor = myResults.getValuesByFactorName(factorsTemp);
			//for (var i=0;i<valuesFactor.length;i++)
			//{
			//	values = values + "<li>"+valuesFactor[i]+"</li>" ; 
			//}
			var names = myResults.getSequencesNamesList();
			window.open("./resultsCol.pl?factor="+factorsTemp+"&values="+ valuesFactor+"&length="+ valuesFactor.length+"&names="+ names);
		//document.getElementById("singleCell").innerHTML="<div id = 'showInfo'> <h2 class='titre'><u>List of " + factorsTemp +"s</u></h2><br/>"+ values + " <br/> </div>"
		}
	}  */
}

/**
 * Gérer couleur
 */
function colorOver(a,b)
{
	for (var i=0;i<rowsNumber+1;i++)
	{
	
		if (b == 8) //gerer que la couleur de la colonne de l'image  
		{
			document.getElementById('cell-'+i+'-'+b).setAttribute('bgcolor','#EDEDED'); 
		}
	}
	for (var j=0;j<columnsNumber+1;j++)
	{
	
		document.getElementById('cell-'+a+'-'+j).setAttribute('bgcolor','#EDEDED');
	}
}

/**
 * Gérer couleur en dehors de la sélection
 */
function colorOut(a,b)
{
	for (var i=0;i<rowsNumber+1;i++)
	{
		document.getElementById('cell-'+i+'-'+b).setAttribute('bgcolor','white');
	}
	for (var j=0;j<columnsNumber+1;j++)
	{
		document.getElementById('cell-'+a+'-'+j).setAttribute('bgcolor','white');
	}
}

/**
 * création du tableau avec les résultats
 */
function createGrid(id,rowsNumber,columnsNumber)
{
	
	//var div  = document.getElementById('select');
	
	//div.id = "select";
	//div.innerHTML = "<p style='font-size:15px' >Export selected entries \( <a onclick='selectAll("+rowsNumber+")' href='#' class='myButton'>Select all<\/a> /          <a  onclick='deSelectAll("+rowsNumber+")'  href='#' class='myButton'>Deselect all</a> \)</p ><p>rrrrfdd</p><p>rrrrfdd</p>";
	
	var tar=document.getElementById(id); // div "table"
	//tar.appendChild(div);
	var table=document.createElement('table');
	
	
	table.border='1';
	var tbdy=document.createElement('tbody');
	tbdy.id = 'cases';
	table.appendChild(tbdy); //ajouter au tableau table
	for (var i=0;i<rowsNumber;i++) // parcours nombre de ligne
	{
		var tr=document.createElement('tr'); 
		tbdy.appendChild(tr); //ajouter ligne
		if (i!=0) 	
		{
			var checkbox = document.createElement('input');
			checkbox.type = "checkbox";
			checkbox.name = "name";
			checkbox.value = "value";
			checkbox.id = "checkbox"+i;
			checkbox.setAttribute("class" , "allCheckbox");
			var td=document.createElement('td');
			tr.appendChild(td);
			td.appendChild(checkbox);
		}else {

			var checkbox = document.createElement('input');
			checkbox.type = "checkbox";
			checkbox.name = "name";
			checkbox.value = "value";
			checkbox.id = "checkboxSelect";
			checkbox.setAttribute('onclick',"checkboxSelect()");
			var td=document.createElement('td');
			td.appendChild(checkbox);
			
			tr.appendChild(td);
		
		
		}
		for (var j=0;j<columnsNumber;j++) // parcours nombre de colonnes
		{
			if (i==0) 
			{
				var td=document.createElement('th');
				td.setAttribute("class" , "factors");
			}else if (j==0) {
				var td=document.createElement('td');
				td.setAttribute("class" , "names");
			}
			else
			{
				var td=document.createElement('td'); 
			}
			//var td=document.createElement('td');  //création colonne
			td.setAttribute('id',"cell-"+i+"-"+j);
			td.setAttribute('onclick',"showCellInfo("+i+","+j+")");
			td.setAttribute('onmouseover',"colorOver("+i+","+j+")");
			td.setAttribute('onmouseout',"colorOut("+i+","+j+")");
			td.width='100';
			if ((i==0)&&(j!=0))
			{
				
				
				var value = myResults.getFactorsNamesList()[j-1]; // ajouter critère
				//if ((myResults.getFactorsNamesList()[j-1]) == 'quality' ) {value = value + '<h5>fff</h5>'} ; 
				if  (value.toString() == 'image')  {
					td.innerHTML = '<h3 style="text-transform:capitalize;">2D structure</h3>' ;
				}
				else if (value.toString() == 'quality') 
				{
					td.innerHTML = '<h3 onclick ="sortingTable(\'all2\')"   style="text-decoration:underline;text-transform:capitalize;">'+value.toString()+'</h3>';
				}else 
				{
					td.innerHTML = '<h3 style="text-transform:capitalize;">'+value.toString()+'</h3>' ; // ajouter la valeur à la cellule
				}			
				
			}
			if ((i!=0)&&(j==0)) // ajouter nom séquence 
			{
				 var value = myResults.getSequencesNamesList()[i-1];
				 td.innerHTML = "<a class='name'>"+value + "</a>"; // ajouter la valeur à la cellule
			}
			if ((i!=0)&&(j!=0))  // cas cellule
			{
				var factor = myResults.getFactorsNamesList()[j-1];
				

				if ( factor =='quality') 
				{	
					var value = myResults.getValueByIndices(i-1,j-1); // appel fonction qui définit la valeur à partir des indices 				
					var string = repeat("<img src='/arn/images/Star.png' alt='arobas' style='width:15px; height:15px;' /> 	 ", parseInt(value))
				
					td.innerHTML = string;
				}
				else if ( factor =='image') //cas lien image
				{	
					var value = myResults.getValueByIndices(i-1,j-1);
			
					//td.innerHTML = "<a target='_blank' href='"+ value + "'>"+myResults.getSequencesNamesList()[i-1]+"</a>  "; // ajouter la valeur à la cellule
					td.innerHTML = "<a  href='" +value + "'><a target='_blank' href='./resultsCell.pl?typePage=image&amp;url=" +value + "'><img src='/arn/images/loupe.png' alt='arobas' style='width:15px; height:15px;' /></a></a>  "; // ajouter la valeur à la cellule
					
					//var a=document.createElement('a'); 
					//var value = a.setAttribute("href","../images/".myResults.getValueByIndices(i-1,j-1));
					//var value = document.createTextNode(document.location.href= '../images/contig15750__ss.ps'); // appel fonction qui définit la valeur à partir des indices 
				}
				else if( factor =='alignment') //cas lien alignement 
				{
					var value = myResults.getValueByIndices(i-1,j-1);
					if (value == 1)
					{
						
						//td.innerHTML = "<a target='_blank' href='./resultsCell.pl?typePage=alignement&amp;url="+value+ "'>View alignment</a>  "; // ajouter la valeur à la cellule
						td.innerHTML = "<img src='/arn/images/check.png' alt='arobas' style='width:15px; height:15px;' />"; // ajouter la valeur à la cellule
					}
					else
					{
						td.innerHTML = "<img src='/arn/images/cross.svg' alt='arobas' style='width:15px; height:15px;' />"; // ajouter la valeur à la cellule
					
					}
				}
				else
				{
					var value = document.createTextNode(myResults.getValueByIndices(i-1,j-1)); // appel fonction qui définit la valeur à partir des indices 
					td.appendChild(value); // ajouter la valeur à la cellule
				}
			}
			if ((i==0)&&(j==0))
			{
				var value = 'name';
				td.innerHTML = '<h3 onclick ="sortingTable(\'all\')"   style="text-decoration:underline;text-transform:capitalize;">'+value.toString()+'</h3>';
			}
			
			tr.appendChild(td); // ajoute colonne à la ligne
		}
	}
	tar.appendChild(table);

}




/**
 * 
 */
function checkboxSelect()
{
	if (document.getElementById('checkboxSelect').checked == true)
	{
		selectAll();
	}else
	{
		deSelectAll();
	}
}

/**
 * 
 */
function selectAll()
{
	rowsNumber = myResults.getSequencesNamesList().length; //récupération de la longueur du tableau à afficher(ligne)
	for (var i=1;i<rowsNumber+1;i++)
	{
		document.getElementById('checkbox'+i).checked = true;	
	}
}

/**
 * 
 */
function deSelectAll()
{	
	rowsNumber = myResults.getSequencesNamesList().length; //récupération de la longueur du tableau à afficher(ligne)
	for (var i=1;i<rowsNumber+1;i++)
	{
		document.getElementById('checkbox'+i).checked = false;	
	}
}

/**
 * Function to repeat String
 */
function repeat(str, times) {
    return new Array(times + 1).join(str);
}


/**
 * Retrieve the names of the sequences selected
 * @return {Array} An array of the selected name
 */
function getChecked()
{
	var tab=new Array();
	for (var i=1;i<rowsNumber;i++)
	{
		if (document.getElementById('checkbox'+i).checked == true )
		{
			var nameTemp = myResults.getSequenceNameByIndex(i-1);
			var positions = myResults.getValuesByFactorName("position");
			var factorsTemp = myResults.getSequenceByNameFactors(nameTemp,positions[i-1]);
			tab.push(nameTemp+"__"+factorsTemp.position);
		}
	}
	return tab;
}


function exportTo(id )
{
	var radios = document.getElementsByName('export');
	for (var i = 0, length = radios.length; i < length; i++) {
    		if (radios[i].checked) {
			// do whatever you want with the checked radio
			var checked = radios[i].value;
		
			// only one radio can be logically checked, don't check the rest
			break;
		}
	}
	if (checked == 'csv') {exportCSV(id);};	
	if (checked == 'odf') {exportODT(id);};
	if (checked == 'gff') {exportGFF(id);};	
}


/**
 * Export the selected sequences as CSV
 * @param {String} id Identifier of the job
 */
function exportCSV(id)
{	
	
	var tab=new Array();
	tab = getChecked();
	if(tab.length ==0)
	{
	   alert("You must select at least one element to export!! ");
	   return;
	}
	window.location = "/cgi-bin/resultsAsCSV.pl?" + "data=" + tab.join(',') + "&run_id=" + id;
}

/**
 * Export the selected sequences as ODT
 * @param {String} id Identifier of the job
 */
function exportODT(id)
{
	window.location = "/cgi-bin/resultsAsODT.pl?" + "&run_id=" + id;
}

/**
 * Export the selected sequences as GFF
 * @param {String} id Identifier of the job
 */
function exportGFF(id)
{
	window.location = "/cgi-bin/resultsAsGFF.pl?" + "&run_id=" + id;
}


function sortingTable(id)
{
	var table = document.getElementById('table');
	table.innerHTML = "";
	main(id);
}


function changeValue()
{
	if (document.getElementById('sort').value == 'Sort by position')
		
	{
		document.getElementById('sort').value = 'Sort by quality';
		sortingTable('all');
	} else 
	{
		document.getElementById('sort').value = 'Sort by position';
		sortingTable('all2');
	}
	 	
}

