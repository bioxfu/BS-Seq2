WT=Col_0
MUT1=suvh13
MUT2=suvh1378
MUT3=ros1_4
MUT4=duf6

WT_1=${WT}_1
WT_2=${WT}_2
MUT1_1=${MUT1}_1
MUT1_2=${MUT1}_2
MUT2_1=${MUT2}_1
MUT2_2=${MUT2}_2
MUT3_1=${MUT3}_1
MUT3_2=${MUT3}_2
MUT4_1=${MUT4}_1
MUT4_2=${MUT4}_2
NUM_LIB=10

MUT=( $MUT1_1 $MUT1_2 $MUT2_1 $MUT2_2 $MUT3_1 $MUT3_2 $MUT4_1 $MUT4_2 )
ALL=( $MUT1_1 $MUT1_2 $MUT2_1 $MUT2_2 $MUT3_1 $MUT3_2 $MUT4_1 $MUT4_2 $WT_1 $WT_2 )
SAM=( $MUT1 $MUT2 $MUT3 $MUT4 )

TXDB=/cluster/home/xfu/Gmatic7/gene/tair10/txdb/tair10_txdb.sqlite
GENE_ANNO=/cluster/home/xfu/Gmatic7/gene/tair10/tair10_gene_anno.tsv
IGV=/cluster/home/xfu/igv/genomes/tair10.genome
GENE=/cluster/home/xfu/Gmatic7/gene/tair10/tair10_gene.bed
TE=/cluster/home/xfu/Gmatic7/gene/tair10/tair10_TE.bed


## Retaining the cytosines that have depth >=4 in all libraries
mkdir dep4
./script/step1_retain.pl dep4 4 $NUM_LIB \
count/${WT_1}_methylome_all.txt   ${WT_1}   count/${WT_2}_methylome_all.txt   ${WT_2} \
count/${MUT1_1}_methylome_all.txt ${MUT1_1} count/${MUT1_2}_methylome_all.txt ${MUT1_2} \
count/${MUT2_1}_methylome_all.txt ${MUT2_1} count/${MUT2_2}_methylome_all.txt ${MUT2_2} \
count/${MUT3_1}_methylome_all.txt ${MUT3_1} count/${MUT3_2}_methylome_all.txt ${MUT3_2} \
count/${MUT4_1}_methylome_all.txt ${MUT4_1} count/${MUT4_2}_methylome_all.txt ${MUT4_2} 

## total methylation level
echo -e "sample\tmC\tC\tmethy_level" > dep4/total_methy_level.tsv
printf '%s\n' "${ALL[@]}"|xargs -I {} ./script/total_methy_level.sh {} 0 >> dep4/total_methy_level.tsv

## Idenfiying DMRs
mkdir DMR
## For Arabidopsis (for other species, please modifiy the chrom size in the script)
printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step2_DMR_tair10.pl dep4/${WT_1}_dep4_Meth.txt dep4/{}_dep4_Meth.txt DMR {}_vs_${WT_1}"
printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step2_DMR_tair10.pl dep4/${WT_2}_dep4_Meth.txt dep4/{}_dep4_Meth.txt DMR {}_vs_${WT_2}"
find DMR/*_list.txt|xargs -I {} wc -l {}|sed 's/_list.txt//'|sed 's/DMR\///'|sed -r 's/_hyper/\thyper/'|sed -r 's/_hypo/\thypo/'|awk '{print $2"\t"$3"\t"$1-1}' > DMR/DMR_number_replicates.tsv

## merge DMRs
printf '%s\n' "${SAM[@]}"|parallel --gnu "cat DMR/{}_*_vs_*_hypo_list.txt|cut -f1-3,21,22|grep -v 'avge_meth_level'|sortBed|mergeBed > DMR/{}_merged_DMR_hypo.bed"
printf '%s\n' "${SAM[@]}"|parallel --gnu "cat DMR/{}_*_vs_*_hyper_list.txt|cut -f1-3,21,22|grep -v 'avge_meth_level'|sortBed|mergeBed > DMR/{}_merged_DMR_hyper.bed"
find DMR/*.bed|xargs -I {} wc -l {}|sed 's/.bed//'|sed 's/DMR\///'|sed -r 's/_hyper/\thyper/'|sed -r 's/_hypo/\thypo/'|awk '{print $2"\t"$3"\t"$1}' > DMR/DMR_number_merged.tsv

## plot DMR distribution on chromosomes
printf '%s\n' "${SAM[@]}"|parallel --gnu "python script/plot_chrom.py {}"

## draw composition of DMR using pie graph
printf '%s\n' "${SAM[@]}"|parallel --gnu "./script/DMR_composition.sh DMR/{}_merged_DMR_hypo.bed $GENE $TE"
printf '%s\n' "${SAM[@]}"|parallel --gnu "./script/DMR_composition.sh DMR/{}_merged_DMR_hyper.bed $GENE $TE"

## draw Venn diagram of DMR
python script/venn_maker.py DMR/${MUT1}_merged_DMR_hyper.bed,DMR/${MUT2}_merged_DMR_hyper.bed,DMR/${MUT3}_merged_DMR_hyper.bed ${MUT1},${MUT2},${MUT3} figure/DMR_hyper_${MUT1}_${MUT2}_${MUT3}.png
python script/venn_maker.py DMR/${MUT1}_merged_DMR_hypo.bed,DMR/${MUT2}_merged_DMR_hypo.bed,DMR/${MUT3}_merged_DMR_hypo.bed ${MUT1},${MUT2},${MUT3} figure/DMR_hypo_${MUT1}_${MUT2}_${MUT3}.png

## build track for IGV
mkdir track
printf '%s\n' "${ALL[@]}"|parallel --gnu "script/dep4_to_bedgraph.sh dep4/{}_dep4_Meth.txt > track/{}_dep4_Meth.bedgraph"
printf '%s\n' "${ALL[@]}"|parallel --gnu "igvtools toTDF track/{}_dep4_Meth.bedgraph track/{}_dep4_Meth.tdf $IGV"

## Calculating the methylation level
#mkdir meth
#printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step3_cal_meth.pl dep4 DMR/{}_vs_${WT_1}_hypo_list.txt meth {}_vs_${WT_1}_hypo_meth_level.txt"
#printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step3_cal_meth.pl dep4 DMR/{}_vs_${WT_2}_hypo_list.txt meth {}_vs_${WT_2}_hypo_meth_level.txt"
#printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step3_cal_meth.pl dep4 DMR/{}_vs_${WT_1}_hyper_list.txt meth {}_vs_${WT_1}_hyper_meth_level.txt"
#printf '%s\n' "${MUT[@]}"|parallel --gnu "./script/step3_cal_meth.pl dep4 DMR/{}_vs_${WT_2}_hyper_list.txt meth {}_vs_${WT_2}_hyper_meth_level.txt"

## Annotation of DMR
#mkdir anno
#printf '%s\n' "${MUT[@]}"|parallel --gnu "Rscript script/DMR_anno.R DMR/{}_vs_${WT}_hyper_list.txt $TXDB $GENE_ANNO anno/{}_vs_${WT}_hyper_list_anno.xls anno/{}_vs_${WT}_hyper_list_annoPie.pdf"
#printf '%s\n' "${MUT[@]}"|parallel --gnu "Rscript script/DMR_anno.R DMR/{}_vs_${WT}_hypo_list.txt $TXDB $GENE_ANNO anno/{}_vs_${WT}_hypo_list_anno.xls anno/{}_vs_${WT}_hypo_list_annoPie.pdf"

## Hyper/Hypo table
#cut -f1 anno/*hyper*.xls|grep -v 'geneId'|sort|uniq > anno/hyper_geneId
#cut -f1 anno/*hypo*.xls|grep -v 'geneId'|sort|uniq > anno/hypo_geneId
#Rscript script/hyper_hypo_table.R $GENE_ANNO
#rm anno/hyper_geneId anno/hypo_geneId

## pack all the results
#DATE=`date +"%Y%m%d"`
#mkdir -p results/${DATE}
#cp -r anno meth track results/${DATE}
#tar czf results/BS-Seq_${DATE}.tar.gz results/${DATE}
