#!/bin/bash

# commande awk qui elimine les retours a la ligne dans une sequence fasta 
awk '{if (NR==1 && $0 ~/>/){print$0;next}if($0~/^>/){print"\n"$0;next}else{printf("%s",$0)}}' $1 > $2
rm $1
