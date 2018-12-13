import pybedtools
import sys, os

beds = sys.argv[1].split(',')
names = sys.argv[2].split(',')
output = sys.argv[3]
#beds = ['DMR/ros1_4_merged_DMR_hyper.bed', 'DMR/suvh13_merged_DMR_hyper.bed', 'DMR/suvh1378_merged_DMR_hyper.bed']
#names = ['ros1_4', 'suvh13', 'suvh1378']

args = {}
args['2'] = ['imagetype="png"', 'col="transparent"', 'fill=c("cornflowerblue","green")', 'alpha=0.50', 'cex=1.5', 'cat.cex=1.5', 'euler.d=TRUE', 'scaled=TRUE', 'cat.col=c("darkblue","darkgreen")']
args['3'] = ['imagetype="png"', 'col="transparent"', 'fill=c("cornflowerblue","green","yellow")', 'alpha=0.50', 'cex=1.5', 'cat.cex=1.5', 'euler.d=TRUE', 'scaled=TRUE', 'cat.col=c("darkblue","darkgreen","orange")']
args['4'] = ['imagetype="png"', 'col="transparent"', 'fill=c("cornflowerblue","green","yellow","darkorchid1")', 'alpha=0.50', 'cex=1.5', 'cat.cex=1.5', 'euler.d=TRUE', 'scaled=TRUE', 'cat.col=c("darkblue","darkgreen","orange","darkorchid4")', 'label.col=c("orange","white","darkorchid4","white","white","white","white","white","darkblue","white","white","white","white","darkgreen","white")']

pybedtools.contrib.venn_maker.venn_maker(beds, names=names, figure_filename=output, script_filename='venn.R', additional_args=args[str(len(beds))], run=True)

os.system("rm venn.*")
os.system("rm "+output+"*.log")
