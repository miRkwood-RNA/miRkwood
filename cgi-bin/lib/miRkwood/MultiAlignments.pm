package miRkwood::MultiAlignments;

# ABSTRACT: Module to create multiple alignments from 2 by 2 alignments

use strict;
use warnings;

# Feature set that will help create a multiple alignment from alignment 2 by 2. 
# The depart is a hash (key : string with a basic position of candidat in the alignment 2 by 2, 
# ex:"5-24"; value : table of hash with the data for each alignment 2 by 2). 
# A 2 dimensional array will be created, it will contain the nucleotide positions 
# for miRBase sequence depending on the candidate sequence alignment 2 by 2. 
# Finally, a 2 dimensional array is created, it will represent the multiple alignment.


#Modify the sequence of alignment 2*2. Because, the position of candidates in the alignment can be different of basic position. Gaps or nucl can be add. To add a nucl in the candidate, the X (it's just a letter random) is choose, it's not important for the suite.
#Parameters : $cdtSeq (sequence of candidate), $miRSeq (sequence of mirBase), $posBaseCdtBegin (it's the begin basic position for the aligment), $posBaseCdtEnd (it's the end basic position for the aligment), $tabDataAligtRef (table with the data for a alignment 2*2), $line (position of candidate in $tab2DRef), $tab2DRef (table 2D with the position for the sequence  mirBase of alignment 2*2 compared sequence candidate)
#Return : $cdtSeq, $miRSeq, @tab2D
sub modifyAligt{
    my ($cdtSeq, $miRSeq, $posBaseCdtBegin,$posBaseCdtEnd, $tabDataAligtRef, $line, $tab2DRef) = @_;
    my @tabDataAligt = @{$tabDataAligtRef};
    my @tab2D = @{$tab2DRef};
    my $posCdtBegin = $tabDataAligt[0][$line]{'begin_target'};
    my $posCdtEnd = $tabDataAligt[0][$line]{'end_target'};
    if ($posCdtBegin>$posBaseCdtBegin){
        for(my $i=0;$i<($posCdtBegin-$posBaseCdtBegin);$i++){
            $cdtSeq='X'.$cdtSeq;
            $miRSeq='-'.$miRSeq;
        }
    }
    if ($posCdtEnd<$posBaseCdtEnd){
        for(my $i=0;$i<($posBaseCdtEnd-$posCdtEnd);$i++){
            $cdtSeq=$cdtSeq.'X';
            $miRSeq=$miRSeq.'-';
        }
    }
    if ($posBaseCdtEnd<$posCdtEnd){
        my $nbSeq = @{$tabDataAligt[0]};
        for(my $j=0; $j<($posCdtEnd-$posBaseCdtEnd);$j++){
            for( my $l=0; $l<$nbSeq; $l++){
            $tab2D[$l][($posCdtEnd-$posBaseCdtBegin)+$j]=0;
            }
        }
    }
    return ($cdtSeq, $miRSeq, \@tab2D);
}

#Indicate by a star if the column has the same nucleotide for each line.
#Parameter : @tabAlgtMult (multiple alignment)
#return : @tabAlgtMult
sub starInAlgt{
    my $tabAlgtRef = $_[0];
    my @tabAlgtMult = @{$tabAlgtRef};
    my $nbLine = scalar(@tabAlgtMult);
    for (my $col=0;$col<@{$tabAlgtMult[0]};$col++){
        my $cpt = 1;
        while($cpt<$nbLine && $tabAlgtMult[$cpt][$col] eq $tabAlgtMult[0][$col]){
            $cpt++;
        }
        if ($cpt == $nbLine){
            $tabAlgtMult[$nbLine][$col]='*';
        }
        else{
            $tabAlgtMult[$nbLine][$col]=' ';
        }
    }
    return @tabAlgtMult;
}

#Write results in file in the directory alignments. The name file is "id_aln.txt"
#Parameter :  $tabNameRef (table with the name of mirBase sequence), $tabAlgtMultRef (table two dimensions containing the aligt multiple), $posBeginCdt (position of begin for the candidate), $posEndCdt (position of end for the candidate)
sub writeInFile{
    my ($tabNameRef, $tabAlgtMultRef, $posBeginCdt, $posEndCdt, $id_candidate) = @_;
    my @tabName = @{$tabNameRef};
    my @tabAlgtMult = @{$tabAlgtMultRef};
    my $cfg = miRkwood->CONFIG();
    my $aln_dir = miRkwood::Paths::get_dir_alignments_path_from_job_dir( $cfg->param('job.directory') );
    my $output = '';
    
    for(my $i=0; $i<@tabAlgtMult; $i++){
		$output .= "\t";
		if ($i==0){
			$output .= 'query       '.$posBeginCdt.'  ';
		}
		elsif ($i==@tabAlgtMult-1){
				if($posBeginCdt>9 && $posBeginCdt < 100){
					$output .= '                ';
				}
				elsif($posBeginCdt>99){
					$output .= '                 ';
				}
				else{
					$output .= '               ';
				}
		}
		else{
			if($posBeginCdt>9 && $posBeginCdt < 100){
				if($i>9){
					$output .= 'miRBase '.$i.'	     ';
				}
				else{
					$output .= 'miRBase '.$i.'       ';
				}
			}
			elsif($posBeginCdt>99){
				if($i>9){
					$output .= 'miRBase '.$i.'	    ';
				}
				else{
					$output .= 'miRBase '.$i.'	     ';
				}
			}
			else{
				if($i>9){
					$output .= 'miRBase '.$i.'	     ';
				}
				else{
					$output .= 'miRBase '.$i.'	      ';
				}
			}
		}
		for(my $j=0; $j<@{$tabAlgtMult[0]}; $j++){
			$output .= $tabAlgtMult[$i][$j];
		}
		if ($i==0){
			$output .= '  '.$posEndCdt."\n";
		}
		else{
			$output .= "\n";
		}
    }
    $output .= "\n";
    for(my $cpt=1; $cpt<@tabName; $cpt++){
        $output .= 'miRBase '.$cpt.': '.$tabName[$cpt]."\n";
    }
    $output .= "\n\n\n";
    open(my $FILE,'>>', $aln_dir."/${id_candidate}_aln.txt") or die"open: $!";
    print $FILE $output;
    close($FILE);
}


#Compares the end positions between the basic candidate and current, if the end position of current candidate is bigger, the basic sequence is update by adding of the last nucleotide of current candidate.
#Parameter : $posBase (end position of basic candidate), $posAlgt (end position of current candidate), $cdtBase (sequence of basic candidate), $cdtCurrent (sequence of current candidate)
#Return : $cdtBase
sub updateCdtBase{
    my ($posBase, $posAlgt, $cdtBase, $cdtCurrent) =@_;
    if ($posAlgt>$posBase){
	for (my $i=0; $i<($posAlgt-$posBase); $i++){
	    $cdtBase.=substr($cdtCurrent, length($cdtCurrent)-(($posAlgt-$posBase)+$i),1);
	}
	$posBase=$posAlgt;
    }
    return ($cdtBase, $posBase);
}

#Split a string to recover the sequence of candidate and mirBase in a alignment 2 to 2
#Parameter : $aligt (alignement 2 to 2)
#Return : $cdtSeq (sequence of candidate), $miRSeq (sequence of mirBase)
sub splitSeqCdtAligt{
    my (@args) = @_;
    my $aligt = shift @args;
    my @tabAligt = split(/\n/, $aligt);
    my $cdtSeq = $tabAligt[0];
    my $miRSeq = $tabAligt[2];
    return ($cdtSeq, $miRSeq);
}

#Split a string to recover the positions of candidate the alignments 2 to 2
#Parameter : $pos (string with position begin and end)
#Return : $cdtSeq (sequence of candidate), $miRSeq (sequence of mirBase)
sub splitPosCdtAligt{
    my (@args) = @_;
    my $pos = shift @args;
    my @tabPos = split(/-/, $pos);
    my $posBegin = $tabPos[0];
    my $posEnd = $tabPos[1];
    return ($posBegin, $posEnd);
}



#Allows to fill the table 2D with the position of the sequences from the alignment 2*2
#Parameter : $hashDataRef (hash table containing the data for the alignment 2*2, key : position of candidate in alignment 2*2, value : table of hash with data)
#return : @hashPosAlgtMult (hash table with like key : position of candidate in alignment 2*2, value : table 2D of alignment multiple)
sub fillTabTemp2D{
    my (@args) = @_;
    my $hashDataRef = shift @args;
    my $id_candidate = shift @args;
    my %hashData = %$hashDataRef;
    while (my ($key, $value) = each(%hashData)){
        my @tab = $value;
        my ($posBeginBase, $posEndBaseCurrent) = splitPosCdtAligt($key);
        my @tabTemp2D;
        my %hashMirSeq;
        my @tabNameMir;
        my $posBaseEndFinal;
        my ($cdtBase, $mirTp) = splitSeqCdtAligt($tab[0][0]{'alignment'});
        $tabNameMir[0] = 'Query';
        for(my $i=0; $i<@{$tab[0]}; $i++){
            my ($cdtSeq, $mirSeq) = splitSeqCdtAligt($tab[0][$i]{'alignment'});
            ($cdtBase, $posBaseEndFinal) = updateCdtBase($posEndBaseCurrent, $tab[0][$i]{'end_target'}, $cdtBase, $cdtSeq);
            $tabNameMir[$i+1] = $tab[0][$i]{'name'};
            ($cdtSeq, $mirSeq, my $tabTemp2DRef) = modifyAligt($cdtSeq, $mirSeq, $posBeginBase, $posEndBaseCurrent, \@tab, $i, \@tabTemp2D);
            @tabTemp2D = @$tabTemp2DRef;
            $hashMirSeq{$tab[0][$i]{'name'}}=$mirSeq;
            @tabTemp2D = setTabTemp($cdtSeq, $mirSeq, $i+1, \@tabTemp2D);
        }
        $hashMirSeq{'Query'}=$cdtBase;
        @tabTemp2D = setTabTemp($cdtBase, $cdtBase, 0, \@tabTemp2D);
            my @tabAlgtMult = setAligtMultiple(\%hashMirSeq, \@tabTemp2D, \@tabNameMir);
        @tabAlgtMult = starInAlgt(\@tabAlgtMult);
        writeInFile(\@tabNameMir, \@tabAlgtMult, $posBeginBase, $posBaseEndFinal, $id_candidate);
    }
    return;
}


#Filled the table 2D for each alignment 2*2. This table contains the nucleotide position of the sequence MirBase compared the candidate sequence in the alignement 2*2
#Parameters : $cdtSeq (sequence of candidate), $miRSeq (sequence of mir), $numMir (number of mir in the tab2D), $tab2DRef (table 2D with the position of sequence)
#Return : @tab2D
sub setTabTemp{
    my ($cdtSeq, $miRSeq, $numMir, $tab2DRef) = @_;
    my $posMiR=0;
    my $posCandidate=0;
    my @tab2D = @$tab2DRef;
    my $i=0;
    
    while($posCandidate<length($cdtSeq)){
        if(substr($miRSeq,$posMiR,1) eq "-"){
            $tab2D[$numMir][$i]=0;
            $i++;
        }
        elsif(substr($cdtSeq, $posCandidate, 1) eq "-"){
        }
        else{
            $tab2D[$numMir][$i]=$posCandidate+1;
            $i++;
        }
        $posCandidate++;
        $posMiR++;
    }

    return @tab2D;
}

#Add a gap at all sequences to the position $posInMir expect the sequence of the miR with a insertion. For this case, the function addNuclAlgt is called.
#If in the position to fill there is a nucleotide, it is shift of 1 column to the right.  
#Parameters : $col (position of alignment), $line (number of the sequence), $posInMir (position of nucleotide to add in the alignment), $hashMirSeqRef (hashtable with key = name of MirBase and value = seq), $tabAlgtMutl (table containing the alignment multiple), $tabNameMirRef (table with the name of mirBase sequences in the alignment 2*2), $nbLine (number of sequence in the alignment multiple), $tabTemp2DRef (table with the position for the alignment 2*2), $colCurrent (column where the insert is create), $gapColumn (number of column of gap already add
# Return @tabTemp2D and @tabAlgtMult
sub addGapColumn{
    my ($col, $line, $posInMir, $hashMirSeqRef, $tabAlgtMultRef, $tabNameMirRef, $nbLine, $tabTemp2DRef, $colCurrent, $gapColumn) = @_;
    my %hashMirSeq = %$hashMirSeqRef;
    my @tabAlgtMult = @$tabAlgtMultRef;
    my @tabNameMir = @$tabNameMirRef;
    my @tabTemp2D = @$tabTemp2DRef;
    
    for (my $i=0; $i<$nbLine; $i++){
        if($i==$line){
            @tabAlgtMult = addNuclAlgt($col, $line, $posInMir-1, $hashMirSeqRef, $tabAlgtMultRef, $tabNameMirRef);
            @tabAlgtMult = addNuclAlgt($col+1, $line, $posInMir, $hashMirSeqRef, \@tabAlgtMult, $tabNameMirRef);
        }
        elsif($i>$line && ($tabTemp2D[$i][$colCurrent]-($tabTemp2D[$i][$colCurrent-1]+$gapColumn)>1)){
            my $diff = $tabTemp2D[$i][$colCurrent]-($tabTemp2D[$i][$colCurrent-1]+$gapColumn);
            my $posMir = $tabTemp2D[$i][$colCurrent]-$diff+1;
            @tabAlgtMult = addNuclAlgt($col, $i, $posMir, $hashMirSeqRef, $tabAlgtMultRef, $tabNameMirRef);
        }
        else{
            if( defined( $tabAlgtMult[$i][$col] ) ){
                $tabAlgtMult[$i][$col+1]=$tabAlgtMult[$i][$col];
            }
            $tabAlgtMult[$i][$col]="-";
        }
    }
    return (\@tabTemp2D, \@tabAlgtMult);
}


#Add the nucleotide of mirBase sequence in the alignment for the position given ($posInMir-1)
#Parameters : $col (position of alignment), $line (number of the sequence), $posInMir (position of nucleotide to add in the alignment), $hashMirSeqRef (hashtable with key = name of MirBase and value = seq), $tabAlgtMutl (table containing the alignment multiple), $tabNameMirRef (table with the name of mirBase sequences in the alignment 2*2) 
#Return @tabAlgtMult : tab of the alignment multiple
sub addNuclAlgt{
    my ($col, $line, $posInMiR, $hashMirSeqRef, $tabAlgtMultRef, $tabNameMirRef) = @_;
    my %hashMirSeq = %$hashMirSeqRef;
    my @tabAlgtMult = @$tabAlgtMultRef;
    my @tabNameMir = @$tabNameMirRef;
    my $nameMiR = $tabNameMir[$line];
    my $seqMiR = $hashMirSeq{$nameMiR};
    $tabAlgtMult[$line][$col] = substr($seqMiR,$posInMiR-1,1);
    return @tabAlgtMult;
}

#Create the alignment Multiple from tabTemp2D containing this position of the sequence in the alignment 2*2. The alignment multiple is a matrix.
#Parameters : %hashMirSeq (key : name mirBase, value : sequence), @tabTemp2D (tab 2 dimensions), @tabNameMirSeq (tab with the name of seq of mirBase)
#Return matrix with Alignment Multiple : @tabAlgtMult
sub setAligtMultiple{
    my ($hashMirSeqRef, $tabTemp2DRef, $tabNameMirRef) = @_;
    my %hashMirSeq = %$hashMirSeqRef;
    my @tabTemp2D = @$tabTemp2DRef;
    my @tabNameMir = @$tabNameMirRef;
    my $nbLine = @tabTemp2D;
    my $nbCol = @{$tabTemp2D[0]};
    my @tabAlgtMult;
    my $gapColumn = 0;
    my $diffCase = 0;
    my $maxDiff = 0;
    my $cptColInsert = 0;
    for(my $i=0; $i<$nbCol; $i++){
        $maxDiff=0;
        $cptColInsert=0;
        for(my $j=0; $j<$nbLine; $j++){
            if ($i==0){
                $diffCase = $tabTemp2D[$j][$i]-0;
            }
            else{
                $diffCase = $tabTemp2D[$j][$i]-$tabTemp2D[$j][$i-1];
            }
            if ($diffCase<=0){
                $tabAlgtMult[$j][$i+$gapColumn]="-";
            }
            elsif ($diffCase>1 && $tabTemp2D[$j][$i-1]!=0 && $diffCase > $maxDiff){
                my $gapAdd = $cptColInsert;
                for(my $nbInsert=1; $nbInsert<$diffCase-$cptColInsert; $nbInsert++){
                    my ($tab2DRef,$tabAlgtRef) = addGapColumn($i+$gapColumn, $j, $tabTemp2D[$j][$i], \%hashMirSeq, \@tabAlgtMult, \@tabNameMir, $nbLine, \@tabTemp2D, $i, $gapAdd);
                    @tabTemp2D=@$tab2DRef;
                    @tabAlgtMult=@$tabAlgtRef;
                    $gapColumn++;
                    $gapAdd++;
                }
                if ($diffCase > $maxDiff){
                    $maxDiff = $diffCase;
                }
                $cptColInsert++;
            }
            else{
                @tabAlgtMult = addNuclAlgt($i+$gapColumn, $j, $tabTemp2D[$j][$i], \%hashMirSeq, \@tabAlgtMult, \@tabNameMir);
            }
        }
    }
    return @tabAlgtMult;
}

1;
