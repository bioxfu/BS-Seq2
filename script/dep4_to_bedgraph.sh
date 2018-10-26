cat $1|grep -v 'meth_level'|awk '{if($3=="-"){$7=$7*-1} if($7!=0)print $1"\t"$2-1"\t"$2"\t"$7}'
