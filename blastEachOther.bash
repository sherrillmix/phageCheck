#!/bin/bash
set -e

echo Making individual db
for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/[^g]*.fasta;do
  cat $ii|makeblastdb -dbtype nucl -out dbs/`basename ${ii%.fasta}` -title `basename ${ii%.fasta}` 
done

#exclude large gg_ sample
for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/[^g]*.fasta;do
  for jj in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/[^g]*.fasta;do
    outfile=work/$(basename $ii|sed 's/\.fasta$//').$(basename $jj|sed 's/\.fasta$//').pair.blast.gz
    if [ -e "$outfile" ];then
      echo $outfile already exists
      continue
    fi
    echo $ii vs. $jj -- $outfile
    blastn -query $ii -db dbs/`basename ${jj%.fasta}` -num_threads 20 -outfmt 6|gzip > $outfile
  done
done
