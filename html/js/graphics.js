/**
 * Instantiation  de classe results.js
 */
function main(id, displayOrphanHairpins,pipeline_type)
{
	myResults = new results(id);
	rowsNumber = myResults.getSequencesNamesList().length; //récupération de la longueur du tableau à afficher(ligne)
	columnsNumber = myResults.getFactorsNamesList().length;//nom de collone
    miRNAFeaturesNumber = myResults.getMiRNAFactorsNamesList().length;
	createGrid('table',rowsNumber+1,columnsNumber+1,displayOrphanHairpins,miRNAFeaturesNumber,pipeline_type);// fonction permettant de créer le tableau par rapport au résultats
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
	if ((i!=0))
	{var identifier = myResults.getIdentifierByIndex(i-1);
		window.open("./displayCandidate.pl?jobID="+id_job+"&id="+identifier);
		var factor = myResults.getFactorsNamesList()[j-1];
	}
}

/**
 * Gérer couleur
 */
function colorOver(a,b)
{
	document.getElementById('checkboxSelect'+a).setAttribute('bgcolor','#EDEDED');
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
	document.getElementById('checkboxSelect'+a).setAttribute('bgcolor','white');
	for (var j=0;j<columnsNumber+1;j++)
	{
		document.getElementById('cell-'+a+'-'+j).setAttribute('bgcolor','white');
	}
}

/**
 * création du tableau avec les résultats
 */
function createGrid(id,rowsNumber,columnsNumber,displayOrphanHairpins,miRNAFeaturesNumber,pipeline_type)
{
	var tar=document.getElementById(id); // div "table"
	//tar.appendChild(div);

    var rowspanValue = 2;
    if ( pipeline_type == 'abinitio')
    {
        rowspanValue = 1;
    }

	var table=document.createElement('table');
	var tbdy=document.createElement('tbody');
	tbdy.id = 'cases';
	table.appendChild(tbdy); //ajouter au tableau table
	for (var i=0;i<rowsNumber;i++) // parcours nombre de ligne
	{
        var tr=document.createElement('tr');
        var td=document.createElement('td');
        tbdy.appendChild(tr); //ajouter ligne
        if (i!=0)
        {
            var checkbox = document.createElement('input');
            checkbox.type = "checkbox";
            checkbox.name = "name";
            checkbox.value = "value";
            checkbox.id = "checkbox"+i;
            checkbox.setAttribute("class" , "allCheckbox");
            td.id = "checkboxSelect"+i;
            tr.appendChild(td);
            td.appendChild(checkbox);
        }else {

            var checkbox = document.createElement('input');
            checkbox.type = "checkbox";
            checkbox.name = "name";
            checkbox.value = "value";
            checkbox.id = "checkboxSelect";
            checkbox.setAttribute('onclick',"checkboxSelect()");
            td.appendChild(checkbox);
            td.setAttribute("rowspan", rowspanValue);
            tr.appendChild(td);
        }

        for (var j=0;j<columnsNumber;j++) // parcours nombre de colonnes
        {
            if (i==0)
            {
                var td=document.createElement('th');
                td.setAttribute("class" , "factors");
            }
            else
            {
                var td=document.createElement('td');
                if (j==0) {
                    td.setAttribute("class" , "names");
                }
                td.setAttribute('onclick',"showCellInfo("+i+","+j+")");
                td.setAttribute('onmouseover',"colorOver("+i+","+j+")");
                td.setAttribute('onmouseout',"colorOut("+i+","+j+")");
            }
            td.setAttribute('id',"cell-"+i+"-"+j);
            td.width='100';

            if ((i==0)&&(j!=0))
            {
                var value = myResults.getFactorsNamesList()[j-1]; // ajouter critère
                //if ((myResults.getFactorsNamesList()[j-1]) == 'quality' ) {value = value + '<h5>fff</h5>'} ; 
                if  (value.toString() == 'image')  {
                    td.innerHTML = '<h3>2D<br>structure</h3>' ;
                }else if (value.toString() == 'precursor_name') 
                {
                    td.innerHTML = '<h3>miRbase name</h3>' ;
                    td.setAttribute("rowspan", rowspanValue);
                }
                else if (value.toString() == 'mirna_sequence') 
                {
                    var td_mirna = document.createElement('th');
                    td_mirna.innerHTML = '<h3>miRNA </h3>' ;
                    td_mirna.setAttribute("colspan", miRNAFeaturesNumber);
                    tr.appendChild(td_mirna); // ajoute colonne à la ligne
                    tbdy.appendChild(tr); // add the first header line to the table

                    var tr=document.createElement('tr');
                    tbdy.appendChild(tr);

                    td.innerHTML = '<h3>sequence</h3>' ;
                }
                else if (value.toString() == 'mirna_length') 
                {
                    td.innerHTML = '<h3>LENGTH</h3>' ;
                }
                else if (value.toString() == 'weight') 
                {
                    td.innerHTML = '<h3>WEIGHT</h3>' ;
                }
                else if (value.toString() == 'strand')
                {
                    td.innerHTML = '<h3>+/-</h3>' ;
                    td.setAttribute("rowspan", rowspanValue);
                }
                else if (value.toString() == 'nb_reads')
                {
                    td.innerHTML = '<h3>READS</h3>' ;
                    td.setAttribute("rowspan", rowspanValue);
                }
                else if (value.toString() == 'reads_distribution')
                {
                    td.innerHTML = '<h3>reads<br>distribution</h3>' ;
                    td.setAttribute("rowspan", rowspanValue);
                }
                else if (value.toString() == 'quality') 
                {
                    td.innerHTML = '<h3 onclick ="sortingTable(\'all2\',false)"   style="text-transform:uppercase;">'+value.toString()+'</h3>';
                    td.setAttribute("rowspan", rowspanValue);
                }
                else if (value.toString() == 'alignment')
                {
                    td.innerHTML = '<h3 >conserved<br>miRNA</h3>';
                }
                else if ((value.toString() == 'mfe') || (value.toString() == 'mfei') || (value.toString() == 'amfe') )
                {
                    td.innerHTML = '<h3 style="text-transform:uppercase;" >'+value.toString()+'</h3>' ;
                    td.setAttribute("rowspan", rowspanValue);
                }
                else 
                {
                    td.innerHTML = '<h3 style="text-transform:uppercase;" >'+value.toString()+'</h3>' ; // ajouter la valeur à la cellule
                    td.setAttribute("rowspan", rowspanValue);
                }

            }
            if ((i!=0)&&(j==0)) // ajouter nom séquence 
            {
                var value = myResults.getSequencesNamesList()[i-1];
                td.innerHTML = "<a class='name'>"+value.substr(0,29)  + "</a>"; // ajouter la valeur à la cellule
            }
            if ((i!=0)&&(j!=0))  // cas cellule
            {
                var factor = myResults.getFactorsNamesList()[j-1];

                if ( factor =='quality') 
                {
                    var value = myResults.getValueByFactor(i-1,'quality'); // appel fonction qui définit la valeur à partir des indices
                    var string = repeat("<img src='/mirkwood/style/star.png' alt='star' style='width:15px; height:15px;' /> 	 ", parseInt(value))

                    td.innerHTML = string;
                }
                else if ( factor =='reads_distribution')
                {
                    var value = myResults.getValueByFactor(i-1,'reads_distribution'); // appel fonction qui définit la valeur à partir des indices
                    var string = repeat("<img src='/mirkwood/style/star.png' alt='star' style='width:15px; height:15px;' /> 	 ", parseInt(value))

                    td.innerHTML = string;
                }
                else if ( factor =='image') //cas lien image
                {
                    var value = myResults.getValueByFactor(i-1,'image');

                    td.innerHTML = "<a  href='" +value + "'><img src='/mirkwood/style/loupe.png' alt='preview' style='width:15px; height:15px;' /></a></a>  "; // ajouter la valeur à la cellule

                }
                else if( factor =='alignment') //cas lien alignement 
                {
                    var value = myResults.getValueByFactor(i-1,'alignment');
                    str = "<img src='/mirkwood/style/check.png' alt='yes' style='width:15px; height:15px;' /> 	 ";
                    if (value == 1)
                    {
                        var string = repeat(str, parseInt(value));
                        td.innerHTML =  string;
                    }
                    else if (value == 2) {
                        var string = repeat(str, parseInt(value));
                        td.innerHTML =  string;
                    }
                    else
                    {
                        td.innerHTML = "";
                    }
                }
                else if( factor =='mfei')
                {
                    var value = myResults.getValueByFactor(i-1,'mfei');
                    if( value < -0.8 )
                    {
                        td.innerHTML = "<font color='#FF00FF'>"+value+"</font>";
                    }
                    else
                    {
                        td.innerHTML = value;
                    }
                }
                else if( factor =='nb_reads')
                {
                    var color = myResults.getValueByFactor(i-1,'criteria_nb_reads');
                    var value = myResults.getValueByFactor(i-1,'nb_reads');
                    if ( color == 1)
                    {
                        td.innerHTML = "<font color='#FF00FF'>"+value+"</font>";
                    }
                    else
                    {
                        td.innerHTML = value;
                    }
                }
                else if( factor =='mirna_sequence')
                {
                    var value = myResults.getValueByFactor(i-1,'mirna_sequence');
                    td.innerHTML = "<span class='mirna_sequence'>"+value+"</span>";
                }
                else
                {
                    var value = document.createTextNode(myResults.getValueByFactor(i-1,factor)); // appel fonction qui définit la valeur à partir des indices 
                    td.appendChild(value); // ajouter la valeur à la cellule
                }
            }
            if ((i==0)&&(j==0))
            {
                td.innerHTML = '<h3 onclick ="sortingTable(\'all\',false)">Chr</h3>';
                td.setAttribute("rowspan", rowspanValue);
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
	var checkboxes = document.getElementsByClassName('allCheckbox');
	for (var i = 0, length = checkboxes.length; i < length; i++) {
		checkboxes[i].checked = true;
	}
}

/**
 * 
 */
function deSelectAll()
{	
	var checkboxes = document.getElementsByClassName('allCheckbox');
	for (var i = 0, length = checkboxes.length; i < length; i++) {
		checkboxes[i].checked = false;
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
	var checkboxes = document.getElementsByClassName('allCheckbox');
	for (var i = 0, length = checkboxes.length; i < length; i++) {
		if ( checkboxes[i].checked == true )
		{
			var re = /checkbox(\d+)/;
			var index = re.exec(checkboxes[i].id)[1];
			var identifier = myResults.getIdentifierByIndex(index-1);
			tab.push(identifier);
		}
	}

	return tab;
}

/**
 * Export the selected sequences
 * @param {String} id Identifier of the job
 */
function exportTo(id, webroot, pipeline, mirna_type)
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

	var tab=new Array();
	tab = getChecked();
	if(tab.length ==0)
	{
	   alert("You must select at least one element to export!! ");
	   return;
	}
	else if(tab.length == rowsNumber ){
		tab = [];
	}
	else if(tab.length > 200 ){
		alert("Please select either all candidates, or less than 200 candidates.");
		return;
	}
	window.location = "/"+webroot+"/exportResults.pl?" + "data=" + tab.join(',') + "&run_id=" + id + "&type=" + checked + "&pipeline=" + pipeline + "&mirna_type=" + mirna_type;
}


function sortingTable(id, displayOrphanHairpins)
{
	document.getElementById('hrefposition').style.color= 'blue';
	document.getElementById('hrefquality').style.color= 'black';
	if (id == 'all2')
	{
		document.getElementById('hrefposition').style.color= 'black';
		document.getElementById('hrefquality').style.color= 'blue';
	}
	if( document.getElementById('displayAll') != null & document.getElementById('dontDisplayOrphanHairpins') != null)
	{
		document.getElementById('displayAll').style.color= 'black';
		document.getElementById('dontDisplayOrphanHairpins').style.color= 'blue';
		if (displayOrphanHairpins)
		{
			document.getElementById('displayAll').style.color= 'blue';
			document.getElementById('dontDisplayOrphanHairpins').style.color= 'black';
		}
	}
	var table = document.getElementById('table');
	table.innerHTML = "";
	main(id, displayOrphanHairpins);
}


function sortBy(sortValue)
{
	var displayOrphanHairpins = true;   // true by default because of abinitio pipeline
	var element = document.getElementById('dontDisplayOrphanHairpins');
	if (element != null && element.className == 'on')
	{
		displayOrphanHairpins = false;
	}
	if (sortValue == 'position')
	{	
		document.getElementById('hrefquality').style.color= 'black';
		document.getElementById('hrefquality').className='off';
		document.getElementById('hrefposition').style.color= 'blue';
		document.getElementById('hrefposition').className='on';
		sortingTable('all',displayOrphanHairpins);
	} else 
	{	
		document.getElementById('hrefposition').style.color= 'black';
		document.getElementById('hrefposition').className='off';
		document.getElementById('hrefquality').style.color= 'blue';
		document.getElementById('hrefquality').className='on';
		sortingTable('all2',displayOrphanHairpins);
	}
	 	
}

