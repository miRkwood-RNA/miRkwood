/**
 * Classe js results
 */
function results(id)
{	
	//attributs
	this.dataXml=document.getElementById(id); // chargement de tout l'XML
	this.SequencesXML = this.dataXml.querySelectorAll("Sequence"); // selection de tous les tags sequences
}

/**
 * Retourner la liste de tous les nom
 */
results.prototype.getSequencesNamesList = function()
{
	var names = new Array();
	var seqs = this.SequencesXML;
	for (var i=0;i<seqs.length;i++)
	{
		names.push(seqs[i].getAttribute('name')) // récupération de l'attribut name de chaque sequence
	}
	return names; // retour tableau name
}

/**
 * Récupération de tous les attributs d'une séquence avec son nom
 */
results.prototype.getSequenceXMLByNameAndPosition = function(name,pos)
{
	
	var seq = this.dataXml.querySelector("Sequence[name='"+name+"'][position='"+pos+"']");
	return seq;
}

/**
 * Retourne un objet contenant les critères de chaque séquence avec son nom
 */
results.prototype.getSequenceByNameFactors = function(name,pos)
{
	var factors = new Object();
	var seq = this.getSequenceXMLByNameAndPosition(name,pos);
	factors.position = seq.getAttribute("position");
	factors.quality = seq.getAttribute("quality");	
	factors.mfei = seq.getAttribute("mfei");
	factors.amfe = seq.getAttribute("amfe");
	factors.mfe = seq.getAttribute("mfe");
	factors.p_value = seq.getAttribute("p_value");
	factors.Vienna = seq.getAttribute("Vienna");
	factors.DNASequence = seq.getAttribute("DNASequence");	
	factors.image = seq.getAttribute("image");	
	return factors;
}

/**
 * Retourne la valeur d'une cellule (nom séquence, critère)
 */
results.prototype.getSequenceByNameByFactor = function(name,factor)
{
	var seq = this.getSequenceXMLByNameAndPosition(name);
	var factor = seq.getAttribute(factor);
	return factor;
}

/**
 * Retourne valeur de la cellule avec les indices des tableaux
 */
results.prototype.getValueByIndices = function(i,j)
{
if (navigator.appName === "Microsoft Internet Explorer") { alert('youpî') ;  };
	var value = this.SequencesXML[i].attributes[j+1].value;
	return value;
}

/**
 * Retourne valeur de la cellule avec un factor
 */
results.prototype.getValueByFactor= function(i,factor)
{

	var value = this.SequencesXML[i].getAttribute(factor);
	return value;
}

/**
 * Avoir toutes les valeurs d'une colonne
 */
results.prototype.getValuesByFactorName = function(factor)
{
	var values = new Array(); 
	var seqs = this.SequencesXML;
	for (var i=0;i<seqs.length;i++)
	{
		values.push(seqs[i].getAttribute(factor))
	}
	return values;
}

/**
 * liste des critères
 */
results.prototype.getFactorsNamesList = function()
{
	var allNames = new Array(  "position","quality" , "mfe","mfei","amfe","p_value", "alignment", "image"  );
	var names = new Array();
	for(var i = 0;i < allNames.length;i++)
	{
		if (document.getElementsByTagName('Sequence')[0].getAttribute(allNames[i])!= null)
		names.push(allNames[i]);
	}
	return names;
	
}

/**
 * Retourner nom séquence à partir d'un indice
 */
results.prototype.getSequenceNameByIndex = function(ind)
{
	return this.SequencesXML[ind].getAttribute("name");
}

/**
 * Retourner critère à partir d'un indice
 */
results.prototype.getFactorNameByIndex = function(ind)
{
	return this.getFactorsNamesList()[ind];
}

/**
 * Retourner indice colonne à partir d'un critere
 */
results.prototype.getIndexByFactorName = function(factor)
{
	var allNames = new Array(  "name", "position","quality" , "mfe","mfei","amfe","p_value", "alignment", "image"  );
	var names = new Array();
	var index = 0;
	for(var i = 0;i < allNames.length;i++)
	{
		if ((document.getElementsByTagName('Sequence')[0].getAttribute(allNames[i])!= null) )
		{
			index++;
			if (allNames[i]== factor) {
				index;break;
			}
		}	
	}
	return index-1;
}