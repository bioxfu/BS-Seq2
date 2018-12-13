import sys
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome

dmr = {'Chr1':[], 'Chr2':[], 'Chr3':[], 'Chr4':[], 'Chr5':[], 'ChrM':[], 'ChrC':[]}

sample_name = sys.argv[1]

with open('DMR/'+sample_name+'_merged_DMR_hyper.bed') as f:
	for line in f:
		chrom, start, end = line.strip().split('\t')
		dmr[chrom].append((int(start), int(end), +1, "", "red"))

with open('DMR/'+sample_name+'_merged_DMR_hypo.bed') as f:
	for line in f:
		chrom, start, end = line.strip().split('\t')
		dmr[chrom].append((int(start), int(end), -1, "", "blue"))

# TAIR10
entries = [("Chr1", 30427671),
           ("Chr2", 19698289),
           ("Chr3", 23459830),
           ("Chr4", 18585056),
           ("Chr5", 26975502)]

max_len = 30427671 # Could compute this from the entries dict
telomere_length = 100000 # For illustration

chr_diagram = BasicChromosome.Organism()
#chr_diagram.page_size = (29.7*cm, 21*cm) # A4 landscape
chr_diagram.page_size = (12*cm, 24*cm) # A4 landscape

for name, length in entries:
	cur_chromosome = BasicChromosome.Chromosome(name)
	# Set the scale to the MAXIMUM length plus the two telomeres in bp,
	# want the same scale used on all five chromosomes so they can be
	# compared to each other
	cur_chromosome.scale_num = max_len + 2 * telomere_length

	# Add an opening telomere
	start = BasicChromosome.TelomereSegment()
	start.scale = telomere_length
	cur_chromosome.add(start)

	# location of the chromosome
	# features = [(5, 1000000, +1, "", "blue")]
	features = dmr[name]

	# Add a body - again using bp as the scale length here.
	body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
	body.scale = length
	cur_chromosome.add(body)

	# Add a closing telomere
	end = BasicChromosome.TelomereSegment(inverted=True)
	end.scale = telomere_length
	cur_chromosome.add(end)

	# This chromosome is done
	chr_diagram.add(cur_chromosome)

chr_diagram.draw("figure/"+sample_name+"_DMR_chrom.pdf", sample_name)
