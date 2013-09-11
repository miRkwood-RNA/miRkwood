/**
 * Instantiation  de classe results.js
 */
function main()
{
	myResults = new results();
	rowsNumber = myResults.getSequencesNamesList().length; //récupération de la longueur du tableau à afficher(ligne)
	columnsNumber = myResults.getFactorsNamesList().length;//nom de collone
	createGrid('table',rowsNumber+1,columnsNumber+1);// fonction permettant de créer le tableau par rapport au résultats
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
		if ((factor != "image") &&(factor != "alignment") &&(factor != "mfei")&&(factor != "mfe")&&(factor != "amfe")&&(factor != "position")&&(factor != "p_value")&&(factor != "self_contain"))
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
	if ((i==0)&&(j!=0))  
	{	
		var factor = myResults.getFactorsNamesList()[j-1];
		if ((factor != "image") &&(factor != "alignment")&&(factor != "amfe")&&(factor != "mfe")&&(factor != "position")&&(factor != "alignment")&&(factor != "alignment")) // quand on clique pas sur image et alignement 
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
	} 
}

/**
 * Gérer couleur
 */
function colorOver(a,b)
{
	for (var i=0;i<rowsNumber+1;i++)
	{
		document.getElementById('cell-'+i+'-'+b).setAttribute('bgcolor','#EDEDED'); 
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
	
	var div  = document.createElement('div');
	div.id = "select";
	div.innerHTML = "<a onclick='selectAll("+rowsNumber+")' href='#' class='myButton'>Select all</a>           <a  onclick='deSelectAll("+rowsNumber+")'  href='#' class='myButton'>Deselect all</a>";

	var tar=document.getElementById(id); // div "table"
	tar.appendChild(div);
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
			
			var td=document.createElement('td');
			tr.appendChild(td);
			td.appendChild(checkbox);
		}else {

	
			var td=document.createElement('td');


			tr.appendChild(td);
		
		
		}
		for (var j=0;j<columnsNumber;j++) // parcours nombre de colonnes
		{
			if (i==0) 
			{
				var td=document.createElement('th');
				td.setAttribute("class" , "factors");
			}else if (j==0) {
				var td=document.createElement('th');
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
				var value = document.createTextNode(myResults.getFactorsNamesList()[j-1]); // ajouter critère
				td.appendChild(value); // ajouter la valeur à la cellule
			}
			if ((i!=0)&&(j==0)) // ajouter nom séquence 
			{
				 var value = myResults.getSequencesNamesList()[i-1];
				 td.innerHTML = "<a class='name'>"+value + "</a>"; // ajouter la valeur à la cellule
			}
			if ((i!=0)&&(j!=0))  // cas cellule
			{
				var factor = myResults.getFactorsNamesList()[j-1];
				if ( factor =='image') //cas lien image
				{	
					var value = myResults.getValueByIndices(i-1,j-1);
					//td.innerHTML = "<a target='_blank' href='"+ value + "'>"+myResults.getSequencesNamesList()[i-1]+"</a>  "; // ajouter la valeur à la cellule
					td.innerHTML = "<a target='_blank' href='./resultsCell.pl?typePage=image&amp;url=" +value + "'>View structure</a>  "; // ajouter la valeur à la cellule
					
					//var a=document.createElement('a'); 
					//var value = a.setAttribute("href","../images/".myResults.getValueByIndices(i-1,j-1));
					//var value = document.createTextNode(document.location.href= '../images/contig15750__ss.ps'); // appel fonction qui définit la valeur à partir des indices 
				}
				else if( factor =='alignment') //cas lien alignement 
				{
					var value = myResults.getValueByIndices(i-1,j-1);
					if (value != "none")
					{
						td.innerHTML = "<a target='_blank' href='./resultsCell.pl?typePage=alignement&amp;url="+value+ "'>View alignment</a>  "; // ajouter la valeur à la cellule
					}
					else
					{
						td.innerHTML = "<SMALL> No alignment</SMALL> ";
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
				var value = document.createTextNode('name');
				td.appendChild(value); // ajouter la valeur à la cellule
			}
			
			tr.appendChild(td); // ajoute colonne à la ligne
		}
	}
	tar.appendChild(table);
}

/**
 * 
 */
function selectAll(rowsNumber)
{
	for (var i=1;i<rowsNumber;i++)
	{
		document.getElementById('checkbox'+i).checked = true;	
	}
}

/**
 * 
 */
function deSelectAll(rowsNumber)
{
	for (var i=1;i<rowsNumber;i++)
	{
		document.getElementById('checkbox'+i).checked = false;	
	}
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


/**
 * Export the selected sequences as CSV
 * @param {String} id Identifier of the job
 */
function exportCSV(id)
{	
	var tab=new Array();
	tab = getChecked();
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
