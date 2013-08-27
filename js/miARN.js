
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

function wainting()
{	
    document.getElementById('upload').disabled="true";
}

function verifySequence()
{	
    var a = document.form.seqArea.value && document.form.seqFile.value;
    var s = !(document.form.seqArea.value || document.form.seqFile.value);
    if (a) alert ("Choose between sequence data and sequence file");
    if (s) alert ("You must provide sequences");
    return !(a||s);
}

function generateExample() {
    var exampleSeqs =
    ['>contig15750',
     'aatgagtaagataaattgctaattaaatgcgacgagaggttcatacatgaagagaagagtgctcttattatgtagccaaggatgaattgcctaatgacagctcaagtcgtttaaaaaacgactctttgttggtttattaggcgttcatttcttgactgacttaatcggctttttttcatcatgttagatcttctcaacttgttacgagcatatcgttcaatattttcatagtcttcttgtaatatgactttgtcaagtcatttcatatagctacttatgtgtagctattattgtcataattattatatagattatatacttaaagagagacttgtaagggatttaagatgtttagataatcatgtaacattcttgtcaagttatgatcaagcattat',
     '>contig15916',
     'aaaaaacctcacatacagcccccgtatctctctctctctataattgataggctattttcttctctctctagaaatgagcttacatggcatgcagatccattgcttatttataggtatagatacagcagatatatattatttattcatatatgtgtatcgaggtatcggaagaagaaattttcattgttacggcggttttctgattcgcttggtgcaggtcgggaacggcttggccgacggtttcatatttgtctccactgtgtgaaacctcgtagcttgagtactgtcctgccttgcatcaactgaatctgaaccgatgtaaatgatctgtgaccggtgtaggagaattggatgaatattgttggagat',
    ].join('\n');
    document.getElementById('seqArea').value = exampleSeqs;
}
