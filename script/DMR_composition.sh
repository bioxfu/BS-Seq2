BED=$1
GENE=$2
TE=$3

bedtools intersect -a $BED -b $GENE -u |awk '{print "gene"}' > $BED.tmp
bedtools intersect -a $BED -b $GENE -v |intersectBed -a - -b $TE -u |awk '{print "TE"}'>> $BED.tmp
bedtools intersect -a $BED -b $GENE -v |intersectBed -a - -b $TE -v |awk '{print "intergenic"}'>> $BED.tmp

sort $BED.tmp|uniq -c|awk '{print $2"\t"$1}' > $BED.comp
Rscript script/DMR_composition.R $BED.comp $BED.comp.pdf

mv $BED.comp.pdf figure
rm $BED.tmp $BED.comp
