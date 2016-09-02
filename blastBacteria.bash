#!/bin/bash
set -e
#download bacteria sequences here:
#http://bacteria.ensembl.org/info/website/ftp/index.html
#wget -r --no-parent --no-directories ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria//fasta/ -A dna.toplevel.fa.gz

echo Making bacteria blast db
if [ -e bacteria/bacteria.00.nin ];then
  cd bacteria
  find -name '*.gz' -exec zcat {} \; |makeblastdb -dbtype nucl -out bacteria -title bacteria
  cd ..
fi

for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/*.fasta;do
  outfile=work/$(basename $ii|sed 's/\.fasta$//').bacteria.blast.gz
  if [ -e "$outfile" ];then
    echo $outfile already exists
    continue
  fi
  echo $ii -- $outfile
  blastn -query $ii -db bacteria/bacteria -num_threads 20 -outfmt 6|gzip > $outfile
done
