function results() // classe js results
{	
	//attributs
	this.dataXml=document.getElementById("all"); // chargement de tout l'XML
	this.SequencesXML = this.dataXml.querySelectorAll("Sequence"); // selection de tous les tags sequences
}
  
results.prototype.getSequencesNamesList = function() // méthode : retourner la liste de tous les nom
{
	var names = new Array();
	var seqs = this.SequencesXML;
	for (var i=0;i<seqs.length;i++)
	{
		names.push(seqs[i].getAttribute('name')) // récupération de l'attribut name de chaque sequence
	}
	return names; // retour tableau name
}

results.prototype.getSequenceXMLByNameAndPosition = function(name,pos) // récupération de tous les attributs d'une séquence avec son nom
{
	
	var seq = this.dataXml.querySelector("Sequence[name='"+name+"'][position='"+pos+"']");
	return seq;
}

results.prototype.getSequenceByNameFactors = function(name,pos) // retourne un objet contenant les critères de chaque séquence avec son nom
{
	var factors = new Object();
	var seq = this.getSequenceXMLByNameAndPosition(name,pos);
	factors.position = seq.getAttribute("position");
	factors.mfei = seq.getAttribute("mfei");
	factors.amfe = seq.getAttribute("amfe");
	factors.mfe = seq.getAttribute("mfe");
	factors.self_contain = seq.getAttribute("self_contain");
	factors.p_value = seq.getAttribute("p_value");
	factors.Vienna = seq.getAttribute("Vienna");
	factors.DNASequence = seq.getAttribute("DNASequence");	
	factors.image = seq.getAttribute("image");	
	return factors;
}

results.prototype.getSequenceByNameByFactor = function(name,factor) //////// retourne la valeur d'une célulle (nom séquence, critère)
{
	var seq = this.getSequenceXMLByNameAndPosition(name);
	var factor = seq.getAttribute(factor);
	return factor;
}

results.prototype.getValueByIndices = function(i,j)// retourne valeur de la cellule avec les indices des tableaux
{
	var value = this.SequencesXML[i].attributes[j+1].value;
	return value;
}

results.prototype.getValuesByFactorName = function(factor) //avoir toutes les valeurs d'une colonne
{
	var values = new Array();
	var seqs = this.SequencesXML;
	for (var i=0;i<seqs.length;i++)
	{
		values.push(seqs[i].getAttribute(factor))
	}
	return values;
}

results.prototype.getFactorsNamesList = function() //liste des critères
{
	var allNames = new Array("position","mfe","mfei","amfe","p_value","self_contain", "alignment","image" );
	var names = new Array();
	for(var i = 0;i < allNames.length;i++)
	{
		if (document.getElementsByTagName('Sequence')[0].getAttribute(allNames[i])!= null)
		names.push(allNames[i]);
	}
	return names;
}

results.prototype.getSequenceNameByIndex = function(ind) // retourner nom séquence à partir d'un indice
{
	return this.SequencesXML[ind].getAttribute("name");
}
results.prototype.getFactorNameByIndex = function(ind)	// retourner critère à partir d'un indice
{
	return this.getFactorsNamesList()[ind];
}


