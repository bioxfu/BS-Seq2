configfile: "config.yaml"

rule all:
	input:
		expand('fastqc/raw/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/raw/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('clean/{sample}_reads1.txt', sample=config['samples']),
		expand('clean/{sample}_reads2.txt', sample=config['samples']),
		expand('fastqc/clean/{sample}_pair1_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_pair2_fastqc.html', sample=config['samples']),
		expand('stat/fastqc_stat.tsv'),
		expand('mapping/{sample}.out.nodupl', sample=config['samples']),
		expand('mapping/{sample}.out.single_mates.nodupl', sample=config['samples']),
		expand('count/{sample}_methylome_all.txt', sample=config['samples']),

rule fastqc_raw_PE:
	input:
		config['path']+'/{sample}_R1.fastq.gz',
		config['path']+'/{sample}_R2.fastq.gz'
	output:
		'fastqc/raw/{sample}_R1_fastqc.html',
		'fastqc/raw/{sample}_R2_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -t 2 -o fastqc/raw {input}'

rule unzip:
	input:
		r1 = config['path']+'/{sample}_R1.fastq.gz',
		r2 = config['path']+'/{sample}_R2.fastq.gz'
	output:
		r1 = 'tmp/{sample}_R1.fastq',
		r2 = 'tmp/{sample}_R2.fastq'
	shell:
		'zcat {input.r1} > {output.r1}; zcat {input.r2} > {output.r2}'

rule trim_PE:
	input:
		r1 = 'tmp/{sample}_R1.fastq',
		r2 = 'tmp/{sample}_R2.fastq'
	output:
		prefix = 'clean/{sample}',
		r1 = 'clean/{sample}_reads1.txt',
		q1 = 'clean/{sample}_pair1.fastq',
		r2 = 'clean/{sample}_reads2.txt',
		q2 = 'clean/{sample}_pair2.fastq',
	shell:
		'touch {output.prefix}; brat_bw-2.0.1/trim -1 {input.r1} -2 {input.r2} -P {output.prefix} -q 30 -L 33 -m 0'

rule fastqc_clean_PE:
	input:
		'clean/{sample}_pair1.fastq',
		'clean/{sample}_pair2.fastq'
	output:
		'fastqc/clean/{sample}_pair1_fastqc.html',
		'fastqc/clean/{sample}_pair2_fastqc.html'
	params:
		conda = config['conda_path']
	shell:
		'{params.conda}/fastqc -t 2 -o fastqc/clean {input}'

rule fastqc_stat_PE:
	input:
		['fastqc/raw/{sample}_R1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/raw/{sample}_R2_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_pair1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_pair2_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule brat_bw:
	input:
		r1 = 'clean/{sample}_reads1.txt',
		r2 = 'clean/{sample}_reads2.txt'
	output:
		f1 = 'mapping/{sample}.out',
		f2 = 'mapping/{sample}.out.single_mates'
	params:
		index = config['index'],
	shell:
		"brat_bw-2.0.1/brat_bw -P {params.index} -1 {input.r1} -2 {input.r2} -pe -o {output.f1} -i 0 -a 1000 -m 2"

rule write_file_names:
	input:
		f1 = 'mapping/{sample}.out',
		f2 = 'mapping/{sample}.out.single_mates'
	output:
		f1 = 'mapping/{sample}_pair_results',
		f2 = 'mapping/{sample}_single_results'
	shell:
		'echo "{input.f1}" > {output.f1}; echo "{input.f2}" > {output.f2}'

rule brat_remove_dupl:
	input:
		f1 = 'mapping/{sample}_pair_results',
		f2 = 'mapping/{sample}_single_results'
	output:
		f1 = 'mapping/{sample}.out.nodupl',
		f2 = 'mapping/{sample}.out.single_mates.nodupl'
	params:
		ref_file = config['ref_file']
	shell:
		"if [ -f {output.f1} ]; then rm {output.f1} {output.f2}; fi; brat_bw-2.0.1/remove-dupl -r {params.ref_file} -p {input.f1} -1 {input.f2}"

rule write_file_names_2:
	input:
		f1 = 'mapping/{sample}.out.nodupl',
		f2 = 'mapping/{sample}.out.single_mates.nodupl'
	output:
		f1 = 'mapping/{sample}_pair_results.nodupl',
		f2 = 'mapping/{sample}_single_results.nodupl'
	shell:
		'echo "{input.f1}" > {output.f1}; echo "{input.f2}" > {output.f2}'

rule brat_acgt_count:
	input:
		f1 = 'mapping/{sample}_pair_results.nodupl',
		f2 = 'mapping/{sample}_single_results.nodupl'
	output:
		'count/{sample}_methylome'
	params:
		ref_file = config['ref_file']
	shell:
		"touch {output}; brat_bw-2.0.1/acgt-count -r {params.ref_file} -P {output} -p {input.f1} -s {input.f2} -B"

rule combine_forw_rev:
	input:
		forw = 'count/{sample}_methylome_forw.txt',
		rev = 'count/{sample}_methylome_rev.txt'
	output:
		'count/{sample}_methylome_all.txt'
	shell:
		'cat {input.forw} {input.rev}|sort -k1,1 -k2,2n|uniq > {output}'
