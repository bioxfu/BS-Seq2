FORW=$1
REV=$2
BED=$3

cat $FORW|awk '{print $3"\t"$4"\t"$4+length($10)}' > ${FORW}.bed
cat $REV |awk '{print $3"\t"$4"\t"$4+length($10)}' > ${REV}.bed
cat ${FORW}.bed ${REV}.bed | sort -k1,1 -k2,2n -k3,3n > $BED
rm ${FORW}.bed ${REV}.bed
