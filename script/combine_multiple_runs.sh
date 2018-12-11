DRY=$1

if [[ "$DRY" == "--dry-run" ]] || [[ "$DRY" == "" ]]; then
	ls ./fastq/*/*/*.gz|sed -r 's/\/[^\/]+$//'|sort|uniq|./script/rush -k 'cat {}/*1.fq.gz > {}/all_runs_1.fastq.gz; rm {}/*1.fq.gz' $DRY
	ls ./fastq/*/*/*.gz|sed -r 's/\/[^\/]+$//'|sort|uniq|./script/rush -k 'cat {}/*2.fq.gz > {}/all_runs_2.fastq.gz; rm {}/*2.fq.gz' $DRY
else
	echo 'the parameter is --dry-run or empty'
fi

