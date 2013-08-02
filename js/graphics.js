
function main() // instantiation  de classe results.js
{
	myResults = new results();
	rowsNumber = myResults.getSequencesNamesList().length; //récupération de la longueur du tableau à afficher(ligne)
	columnsNumber = myResults.getFactorsNamesList().length;//nom de collone
	createGrid('table',rowsNumber+1,columnsNumber+1);// fonction permettant de créer le tableau par rapport au résultats
}	
function showCellInfo(i,j)// focntion permettant d'afficher la valeur d'une cellule
{
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
		window.open("./resultsRow.pl?nameSeq="+nameTemp+"&mfei="+ factorsTemp.mfei+"&pvalue="+factorsTemp.p_value+"&position="+factorsTemp.position+"&Vienna="+factorsTemp.Vienna+"&DNASequence="+factorsTemp.DNASequence+"&image="+factorsTemp.image+"&mfe="+factorsTemp.mfe+"&amfe="+factorsTemp.amfe+"&self_contain="+factorsTemp.self_contain);
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
function colorOver(a,b) // gérer couleur 
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
function colorOut(a,b) // gérer couleur en dehors de la sélection
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
function createGrid(id,rowsNumber,columnsNumber) // création du tableau avec les résultats
{
	var tar=document.getElementById(id); // div "table"
	var table=document.createElement('table');

	
	table.border='1';
	var tbdy=document.createElement('tbody');
	table.appendChild(tbdy); //ajouter au tableau table
	for (var i=0;i<rowsNumber;i++) // parcours nombre de ligne
	{
		var tr=document.createElement('tr'); 
		tbdy.appendChild(tr); //ajouter ligne
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
