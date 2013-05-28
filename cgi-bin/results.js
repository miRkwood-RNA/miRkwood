function results()
{
	this.dataXml=document.getElementById("toto");
	this.SequencesXML = this.dataXml.querySelectorAll("Sequence");
}

results.prototype.getSequecesNamesList = function()
{
	var names = new Array();
	var seqs = this.SequencesXML;
	for (var i=0;i<seqs.length;i++)
	{
		names.push(seqs[i].getAttribute('name'))
	}
	return names;
}

results.prototype.getSequenceXMLByName = function(name)
{
	
	var seq = this.dataXml.querySelector("Sequence[name='"+name+"']");
	return seq;
}

results.prototype.getSequenceByNameFactors = function(name)
{
	var factors = new Object();
	var seq = this.getSequenceXMLByName(name);
	factors.position = seq.getAttribute("position");
	factors.mfei = seq.getAttribute("mfei");
	factors.p_value = seq.getAttribute("p_value");	
	return factors;
}

results.prototype.getSequenceByNameByFactor = function(name,factor)
{
	var seq = this.getSequenceXMLByName(name);
	var factor = seq.getAttribute(factor);
	return factor;
}

results.prototype.getValueByIndices = function(i,j)
{
	var factor = this.SequencesXML[i].attributes[j+1].value;
	return factor;
}

results.prototype.getFactorsByFactorName = function(factor)
{
	var factors = new Array();
	var seqs = this.SequencesXML;
	for (var i=0;i<seqs.length;i++)
	{
		factors.push(seqs[i].getAttribute(factor))
	}
	return factors;
}

results.prototype.getFactorsNamesList = function()
{
	var names = new Array("position","mfei","p_value");
	return names;
}
