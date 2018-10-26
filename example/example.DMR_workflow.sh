WT=Col_0_1
MUT1_1=suvh13_1
MUT1_2=suvh13_2
MUT2_1=suvh1378_1
MUT2_2=suvh1378_2
MUT3_1=ros1_4_1
MUT3_2=ros1_4_2
MUT4_1=duf6_1
MUT4_2=duf6_2
NUM_LIB=9
MUT=( $MUT1_1 $MUT1_2 $MUT2_1 $MUT2_2 $MUT3_1 $MUT3_2 $MUT4_1 $MUT4_2 )

TXDB=/home/xfu/Gmatic7/gene/tair10/txdb/tair10_txdb.sqlite
GENE_ANNO=/home/xfu/Gmatic7/gene/tair10/tair10_gene_anno.tsv


## Retaining the cytosines that have depth >=4 in all libraries
mkdir dep4
./script/step1_retain.pl dep4 4 $NUM_LIB \
count/${WT}_methylome_all.txt ${WT} \
count/${MUT1_1}_methylome_all.txt ${MUT1_1} count/${MUT1_2}_methylome_all.txt ${MUT1_2} \
count/${MUT2_1}_methylome_all.txt ${MUT2_1} count/${MUT2_2}_methylome_all.txt ${MUT2_2} \
count/${MUT3_1}_methylome_all.txt ${MUT3_1} count/${MUT3_2}_methylome_all.txt ${MUT3_2} \
count/${MUT4_1}_methylome_all.txt ${MUT4_1} count/${MUT4_2}_methylome_all.txt ${MUT4_2} 

## Idenfiying DMRs
mkdir DMR
## For Arabidopsis (for other species, please modifiy the chrom size in the script)
printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step2_DMR_tair10.pl dep4/${WT}_dep4_Meth.txt dep4/{}_dep4_Meth.txt DMR {}_vs_${WT}"

## Calculating the methylation level
mkdir meth
printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step3_cal_meth.pl dep4 DMR/{}_vs_${WT}_hypo_list.txt meth {}_vs_${WT}_hypo_meth_level.txt"
printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step3_cal_meth.pl dep4 DMR/{}_vs_${WT}_hyper_list.txt meth {}_vs_${WT}_hyper_meth_level.txt"
find meth/*_meth_level.txt|xargs -I {} wc -l {}|sed 's/_meth_level.txt//'|sed 's/meth\///'|sed -r 's/_hyper/\thyper/'|sed -r 's/_hypo/\thypo/'|awk '{print $2"\t"$3"\t"$1-1}' > meth/DMR_number.txt

## Annotation of DMR
mkdir anno
printf '%s\n' "${MUT[@]}"|parallel --gnu "Rscript script/DMR_anno.R DMR/{}_vs_${WT}_hyper_list.txt $TXDB $GENE_ANNO anno/{}_vs_${WT}_hyper_list_anno.xls anno/{}_vs_${WT}_hyper_list_annoPie.pdf"
printf '%s\n' "${MUT[@]}"|parallel --gnu "Rscript script/DMR_anno.R DMR/{}_vs_${WT}_hypo_list.txt $TXDB $GENE_ANNO anno/{}_vs_${WT}_hypo_list_anno.xls anno/{}_vs_${WT}_hypo_list_annoPie.pdf"

## Hyper/Hypo table
cut -f1 anno/*hyper*.xls|grep -v 'geneId'|sort|uniq > anno/hyper_geneId
cut -f1 anno/*hypo*.xls|grep -v 'geneId'|sort|uniq > anno/hypo_geneId
Rscript script/hyper_hypo_table.R $GENE_ANNO
rm anno/hyper_geneId anno/hypo_geneId

## pack all the results
DATE=`date +"%Y%m%d"`
mkdir -p results/${DATE}
cp -r anno meth results/${DATE}
tar czf results/BS-Seq_${DATE}.tar.gz results/${DATE}

