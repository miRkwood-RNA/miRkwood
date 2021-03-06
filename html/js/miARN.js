/**
 * 
 */
function showHideBlock() 
{
    if(document.getElementById('CDS').checked)
    {
        document.getElementById('menuDb').style.visibility='visible';
        document.getElementById('menuDb').style.display='block';
    }
    else
    {
        document.getElementById('menuDb').style.visibility='hidden';
        document.getElementById('menuDb').style.display='none';
    }
}

/**
 * 
 */
function ResetForm()
{ 
	document.forms.form.reset();
	document.getElementById('menuDb').style.visibility='hidden';
	document.getElementById('menuDb').style.display='none';
}

/**
 * 
 */
function wainting()
{	
    document.getElementById('upload').disabled="true";
}

/**
 * 
 */
function verifySequence()
{	
    var a = document.getElementById('seqArea').value && document.getElementById('seqFile').value;
    var s = !(document.getElementById('seqArea').value || document.getElementById('seqFile').value);
    if (a) alert ("Choose between sequence data and sequence file");
    if (s) alert ("You must provide sequences");
    return !(a||s);
}

/**
 * 
 */
function generateExample() {
    var exampleSeqs =
    ['>sample',
     'GGTATGTATTACTCTTATCATTCATTACTTTGTGCCACGTGCTATATATATTACTCCCTCTGTCCCATTACAGTTGGCCACAATTTTTTCGGCACGGAGATTAAGAAAATGCTTTTGTAAGGTATAAAATGTTAAGGGCCCACCTACTTTTTGAGTTTTATGTGTTAAATTTTTGCTTTTTTTATAAAGGGGCCAACTGTAATGGGACATTCCAAAATGGAAAAATGGCCAACTGTAATGGGACGGAGGGAGTAGTTAGACTCTTGGT',
    ].join('\n');
    document.getElementById('seqArea').value = exampleSeqs;
}

/**
 * Check if the user provides a BED file and either choose a model organism or provide a reference sequence.
 */
function verifyBEDForm()
{	
    var b = !document.getElementById('bedFile').value;
    if (b) {
        alert ("You must provide a BED file.");
        return false;
    }

    var ref = 0;
    if ( document.getElementById('seqArea').value ) ref+=1;
    if ( document.getElementById('species').value ) ref+=1;
    if ( document.getElementById('seqFile').value ) ref+=1;    

    if (ref == 0){
        alert ("You must either provide a reference sequence or choose an organism in the list");
        return false;
    }
    if (ref > 1){
        alert ("Choose between a model organism in the list and a sequence");
        return false;
    }

    return !b;
}

/**
 * 
 */
function resetElementByID(id)
{	
    document.getElementById(id).value = "";
}

