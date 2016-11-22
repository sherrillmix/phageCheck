#!/bin/bash
set -e
#/media/THING1/local/genomeIndexes/blast/02252014_Virgin_VLP_DB/
#http://www.cell.com/cell/pdfExtended/S0092-8674(15)00003-3

cd virgin
makeblastdb -in vlpDB_nr.c95.fasta -dbtype prot
cd ..

for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/*.fasta;do
  outfile=work/$(basename $ii|sed 's/\.fasta$//').virgin.blast.gz
  if [ -e "$outfile" ];then
    echo $outfile already exists
    continue
  fi
  echo $ii -- $outfile
  blastx -query $ii -db virgin/vlpDB_nr.c95.fasta -num_threads 20 -outfmt 6|gzip > $outfile
done


