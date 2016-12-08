#!/bin/bash
set -e
#download bacteria sequences here:
#http://bacteria.ensembl.org/info/website/ftp/index.html
#wget -r --no-parent --no-directories ftp://ftp.ensemblgenomes.org/pub/release-32/bacteria//fasta/ -A dna.toplevel.fa.gz

echo Making bacteria blast db
if [ ! -e green/green.00.nin ];then
  cd green
  #a bit of a monstrosity to add the file name to all the > in the various bacteria files
  zcat gg_13_5.fasta.gz|makeblastdb -dbtype nucl -out green -title green
  cd ..
fi

#exclude large gg_ sample
for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/[^g]*.fasta;do
  outfile=work/$(basename $ii|sed 's/\.fasta$//').green.blast.gz
  if [ -e "$outfile" ];then
    echo $outfile already exists
    continue
  fi
  echo $ii -- $outfile
  blastn -query $ii -db green/green -num_threads 20 -outfmt 6|gzip > $outfile
done
