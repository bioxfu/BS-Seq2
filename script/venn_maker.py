import pybedtools
import sys, os

beds = sys.argv[1].split(',')
names = sys.argv[2].split(',')
output = sys.argv[3]
Rpath = sys.argv[4]
#beds = ['DMR/ros1_4_merged_DMR_hyper.bed', 'DMR/suvh13_merged_DMR_hyper.bed', 'DMR/suvh1378_merged_DMR_hyper.bed']
#names = ['ros1_4', 'suvh13', 'suvh1378']

pybedtools.contrib.venn_maker.venn_maker(beds, names=names, figure_filename=output, script_filename='venn.R', additional_args=None, run=False)
os.system("grep -v 'diagram' venn.R|grep -v 'name'|sed 's/library(VennDiagram)/library(gplots);library(eulerr);library(RColorBrewer)/'|sed -r 's/,$//'|sed 's/^)//' > venn_data.R")
os.system("cat venn_data.R script/template.R > venn2.R")
os.system(Rpath + " venn2.R " + output)
os.system('rm venn*')
