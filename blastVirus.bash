#!/bin/bash
set -e
#download phage sequences here:
# http://www.ebi.ac.uk/genomes/phage.html
# http://www.ebi.ac.uk/genomes/phage.txt
# http://www.ebi.ac.uk/cgi-bin/sva/sva.pl?&do_batch=1

cd phage
makeblastdb -in phage.fasta -dbtype nucl
cd ..

for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/*.fasta;do
  outfile=work/$(basename $ii|sed 's/\.fasta$//').blast.gz
  if [ -e "$outfile" ];then
    echo $outfile already exists
    continue
  fi
  echo $ii -- $outfile
  blastn -query $ii -db phage/phage.fasta -num_threads 20 -outfmt 6|gzip > $outfile
done



#http://www.ncbi.nlm.nih.gov/nuccore/?term=Viruses[Organism]%20NOT%20cellular%20organisms[ORGN]%20NOT%20wgs[PROP]%20NOT%20gbdiv%20syn[prop]%20AND%20(srcdb_refseq[PROP]%20OR%20nuccore%20genome%20samespecies[Filter])
cd allVirus
makeblastdb -in allVirus.fa -dbtype nucl
cd ..

for ii in /media/THING1/alexandra/4Gut_Phage_and_Bacteria_DNA_Modifications/014Blast16s/*.fasta;do
  outfile=work/$(basename $ii|sed 's/\.fasta$//').virus.blast.gz
  if [ -e "$outfile" ];then
    echo $outfile already exists
    continue
  fi
  echo $ii -- $outfile
  blastn -query $ii -db allVirus/allVirus.fa -num_threads 20 -outfmt 6|gzip > $outfile
done
