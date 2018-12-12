SAM=$1

cat dep4/${SAM}_dep4_Meth.txt|awk -v name="$SAM" '{if($5!=0) x+=1} END {print name"\t"x-1"\t"NR-1"\t"(x-1)/(NR-1)}'
