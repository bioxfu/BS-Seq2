source activate gmatic

if [ ! -d fastqc ]; then
	mkdir -p fastq fastqc/raw fastqc/clean clean tmp mapping stat count sam bed figure
fi
