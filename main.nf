#!/usr/bin/env nextflow

IONICE = 'ionice -c2 -n7'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
		'hg38': 'human',
		'rn5': 'rat',
		'rn6': 'rat',
		'mm9': 'mouse',
		'mm10': 'mouse']

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]


libraries = params.libraries.keySet()

make_excluded_regions_arg = {
	genome ->
	return params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}

has_blacklist = {
	genome ->
	params.blacklist.containsKey(genome)
}

get_blacklists = {
	genome ->
	params.blacklist[genome]
}


get_bwa_index = {
	genome ->
	params.bwa_index[genome]
}


get_genome = {
	library ->
	params.libraries[library].genome
}


get_tss = {
	genome ->
	params.tss[genome]
}


get_organism = {
	genome ->
	ORGANISMS[genome]
}

get_macs2_genome_size = {
	genome ->
	MACS2_GENOME_SIZE[genome]
}


get_chrom_sizes = {
	genome ->
	params.chrom_sizes[genome]
}


library_to_readgroups = {
	library ->
	params.libraries[library].readgroups.keySet()
}


library_and_readgroup_to_fastqs = {
	library, readgroup ->
	params.libraries[library].readgroups[readgroup]
}


trim_in = []
fastqc_in = []

for (library in libraries) {
	for (readgroup in library_to_readgroups(library)) {
		fastqs = library_and_readgroup_to_fastqs(library, readgroup)
		first_insert = fastqs['1']
		second_insert = fastqs['2']
		trim_in << [library, readgroup, file(first_insert), file(second_insert)]
		fastqc_in << file(first_insert)
		fastqc_in << file(second_insert)
	}
}

process trim {

	publishDir "${params.results}/trim", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	time '24h'

	input:
	set val(library), val(readgroup), file(fastq_1), file(fastq_2) from Channel.from(trim_in)

	output:
	set val(library), val(readgroup), file("${library}-${readgroup}.1.trimmed.fastq.gz"), file("${library}-${readgroup}.2.trimmed.fastq.gz") into map_in
	set file("${library}-${readgroup}.1.trimmed.fastq.gz"), file("${library}-${readgroup}.2.trimmed.fastq.gz") into fastqc_trim_in

	"""
	${IONICE} cta $fastq_1 $fastq_2 ${library}-${readgroup}.1.trimmed.fastq.gz ${library}-${readgroup}.2.trimmed.fastq.gz
	"""
}

fastqc_in = Channel.from(fastqc_in).mix(fastqc_trim_in).flatten()


process fastqc {

	publishDir "${params.results}/fastqc", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1

	input:
	file(fastq) from fastqc_in

	output:
	file(fastqc_out)

	script:
	fastqc_out = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')
	
	"""
	fastqc $fastq -o .
	"""

}


process map {

	memory '50 GB'
	cpus 12
	errorStrategy 'retry'
	maxRetries 1
	time '48h'

	publishDir "${params.results}/bwa", mode: 'rellink'

	input:
	set val(library), val(readgroup), file(fastq_1), file(fastq_2) from map_in

	output:
	set val(library), file("${library}-${readgroup}.bam") into merge_in

	"""
	bwa mem -I 200,200,5000 -M -t 12 ${get_bwa_index(get_genome(library))} ${fastq_1} ${fastq_2} | samtools sort -m 1g -@ 11 -O bam -T sort_tmp -o ${library}-${readgroup}.bam -
	"""

}


process merge {
	
	publishDir "${params.results}/merge", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	time '5h'

	input:
	set val(library), file(bams) from merge_in.groupTuple(sort: true)

	output:
	set val(library), file("${library}.bam") into markdup_in

	"""
	${IONICE} samtools merge ${library}.bam ${bams.join(' ')}
	"""
}


process mark_duplicates {

	publishDir "${params.results}/mark_duplicates", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	time '5h'

	input:
	set val(library), file(bam) from markdup_in

	output:
	set val(library), file("${library}.md.bam"), file("${library}.md.bam.bai"), file("${library}.metrics") into md_bams
	set val(library), file("${library}.md.bam"), file("${library}.md.bam.bai") into prune_in
	set val(library), file("${library}.md.bam"), file("${library}.md.bam.bai") into ataqv_md_in

	"""
	java -Xmx4g -Xms4g -jar \$PICARD_JAR MarkDuplicates I=$bam O=${library}.md.bam ASSUME_SORTED=true METRICS_FILE=${library}.metrics VALIDATION_STRINGENCY=LENIENT; samtools index ${library}.md.bam
	"""
}


process prune {

	memory '3 GB'
	time '5h'
	errorStrategy 'retry'
	maxRetries 3

	publishDir "${params.results}/prune", mode: 'rellink'

	input:
	set val(library), file(bam), file(bam_index) from prune_in

	output:
	set val(library), file("${library}.pruned.bam") into pruned_bams

	"""
	${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $bam ${AUTOSOMAL_REFERENCES[get_genome(library)].join(' ')} > ${library}.pruned.bam 
	"""

}


process bamtobed {

	time '4h'

	input:
	set val(library), file(bam) from pruned_bams
	
	output:
	set val(library), file("${library}.bed") into beds

	"""
	bedtools bamtobed -i $bam > ${library}.bed
	"""

}


process macs2 {

	publishDir "${params.results}/macs2", mode: 'rellink'
	time '5h'

	input:
	set val(library), file(bed) from beds

	output:
	set val(library), file("${library}_peaks.broadPeak") into blacklist_in
	set val(library), file("${library}_peaks.broadPeak") into ataqv_macs2_in
	set val(library), file("${library}_treat_pileup.bdg") into bigwig_in

	"""
	macs2 callpeak -t $bed --outdir . --SPMR -f BED -n $library -g ${get_macs2_genome_size(get_genome(library))} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
	"""

}


process blacklist_filter_peaks {

	publishDir "${params.results}/macs2", mode: 'rellink'
	time '1h'

	input:
	set val(library), file(peaks) from blacklist_in

	output:
	file("${library}_peaks.broadPeak.noblacklist")

	when:
	has_blacklist(get_genome(library))

	"""
	bedtools intersect -a $peaks -b ${get_blacklists(get_genome(library)).join(' ')} -v > ${library}_peaks.broadPeak.noblacklist
	"""

}


process bigwig {

	time '5h'
	publishDir "${params.results}/bigwig", mode: 'rellink'

	input:
	set val(library), file(bedgraph) from bigwig_in

	output:
	file("${library}.bw")

	"""
	LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph > sorted.bedgraph
	bedClip sorted.bedgraph ${get_chrom_sizes(get_genome(library))} clipped.bedgraph
	bedGraphToBigWig clipped.bedgraph ${get_chrom_sizes(get_genome(library))} ${library}.bw
	rm sorted.bedgraph clipped.bedgraph
	"""	
	
}


ataqv_in = ataqv_md_in.combine(ataqv_macs2_in, by: 0)

process ataqv {
	
	publishDir "${params.results}/ataqv", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	memory '5 GB'
	time '3h'
	
	input:
	set val(library), file(bam), file(bam_index), file(peaks) from ataqv_in
	
	output:
	set file("${library}.ataqv.json.gz"), file("${library}.ataqv.out") into ataqv_out

	"""
	${IONICE} ataqv --peak-file $peaks --name ${library} --metrics-file ${library}.ataqv.json.gz --tss-file ${get_tss(get_genome(library))} ${make_excluded_regions_arg(get_genome(library))} --ignore-read-groups ${get_organism(get_genome(library))} $bam > ${library}.ataqv.out
	"""	

}
