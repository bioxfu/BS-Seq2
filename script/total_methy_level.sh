SAM=$1
CUTOFF=$2

cat dep4/${SAM}_dep4_Meth.txt|awk -v name="$SAM" -v cutoff="$CUTOFF" '{if($7 > cutoff) x+=1} END {print name"\t"x-1"\t"NR-1"\t"(x-1)/(NR-1)}'
