#!/usr/bin/env nextflow

nextflow.enable.dsl=2
IONICE = 'ionice -c2 -n7'

// Generic data
AUTOSOMAL_REFERENCES = ['hg19': (1..22).collect({it -> 'chr' + it}),
			'hg38': (1..22).collect({it -> 'chr' + it}),
			'rn5': (1..20).collect({it -> 'chr' + it}),
			'rn6': (1..20).collect({it -> 'chr' + it}),
			'rn7': (1..20).collect({it -> 'chr' + it}),
			'mm9': (1..19).collect({it -> 'chr' + it}),
			'mm10': (1..19).collect({it -> 'chr' + it})
]

ORGANISMS = ['hg19': 'human', 
		'hg38': 'human',
		'rn5': 'rat',
		'rn6': 'rat',
		'rn7': 'rat',
		'mm9': 'mouse',
		'mm10': 'mouse']

MACS2_GENOME_SIZE = [
    'rn4': 'mm',
    'rn5': 'mm',
    'rn6': 'mm',
    'rn7': 'mm',
    'mm9': 'mm',
    'mm10': 'mm',
    'hg19': 'hs',
    'hg38': 'hs'
]


def make_excluded_regions_arg (genome) {
	params.blacklist[genome].collect({'--excluded-region-file ' + it}).join(' ')
}

def has_blacklist (genome) {
	params.blacklist.containsKey(genome)
}

def get_blacklists (genome) {
	params.blacklist[genome]
}

def get_bwa_index (genome) {
	params.bwa_index[genome]
}

def get_genome (library) {
	params.libraries[library].genome
}

def get_tss (genome) {
	params.tss[genome]
}

def get_organism (genome) {
	ORGANISMS[genome]
}

def get_macs2_genome_size (genome) {
	MACS2_GENOME_SIZE[genome]
}

def get_chrom_sizes (genome) {
	params.chrom_sizes[genome]
}

def library_to_readgroups (library) {
	params.libraries[library].readgroups.keySet()
}

def library_and_readgroup_to_fastqs (library, readgroup) {
	params.libraries[library].readgroups[readgroup]
}


process trim {

	publishDir "${params.results}/trim", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	time '24h'

	input:
	tuple val(library), val(readgroup), path(fastq_1), path(fastq_2)

	output:
	tuple val(library), val(readgroup), path("${library}-${readgroup}.1.trimmed.fastq.gz"), path("${library}-${readgroup}.2.trimmed.fastq.gz")

	"""
	${IONICE} cta $fastq_1 $fastq_2 ${library}-${readgroup}.1.trimmed.fastq.gz ${library}-${readgroup}.2.trimmed.fastq.gz
	"""
}


process fastqc {

	publishDir "${params.results}/fastqc", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1

	input:
	path(fastq)

	output:
	path(fastqc_out)

	script:
	fastqc_out = fastq.getName().replaceAll('.fastq.gz', '_fastqc.zip')
	
	"""
	fastqc $fastq -o .
	"""

}


process bwa {

	memory '50 GB'
	cpus 12
	errorStrategy 'retry'
	maxRetries 1
	time '48h'

	publishDir "${params.results}/bwa", mode: 'rellink'

	input:
	tuple val(library), val(readgroup), path(fastq_1), path(fastq_2)

	output:
	tuple val(library), path("${library}-${readgroup}.bam")

	"""
	bwa mem -I 200,200,5000 -M -t 12 ${get_bwa_index(get_genome(library))} ${fastq_1} ${fastq_2} | samtools sort -m 1g -@ 11 -O bam -T sort_tmp -o ${library}-${readgroup}.bam -
	"""

}


process merge_mapped {
	
	publishDir "${params.results}/merge", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	time '5h'

	input:
	tuple val(library), path(bams)

	output:
	tuple val(library), path("${library}.bam")

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
	tuple val(library), path(bam)

	output:
	tuple val(library), path("${library}.md.bam"), path("${library}.md.bam.bai")

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
	tuple val(library), path(bam), path(bam_index)

	output:
	tuple val(library), path("${library}.pruned.bam")

	"""
	${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $bam ${AUTOSOMAL_REFERENCES[get_genome(library)].join(' ')} > ${library}.pruned.bam 
	"""

}


process bamtobed {

	time '4h'

	input:
	tuple val(library), path(bam)
	
	output:
	tuple val(library), path("${library}.bed")

	"""
	bedtools bamtobed -i $bam > ${library}.bed
	"""

}


process macs2 {

	publishDir "${params.results}/macs2", mode: 'rellink'
	time '5h'

	input:
	tuple val(library), path(bed)

	output:
	tuple val(library), path("${library}_peaks.broadPeak"), emit: peaks
	tuple val(library), path("${library}_treat_pileup.bdg"), emit: bedgraph

	"""
	macs2 callpeak -t $bed --outdir . --SPMR -f BED -n $library -g ${get_macs2_genome_size(get_genome(library))} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
	"""

}


process blacklist_filter_peaks {

	publishDir "${params.results}/macs2", mode: 'rellink'
	time '1h'

	input:
	tuple val(library), path(peaks)

	output:
	path("${library}_peaks.broadPeak.noblacklist")

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
	tuple val(library), path(bedgraph)

	output:
	path("${library}.bw")

	"""
	LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph > sorted.bedgraph
	bedClip sorted.bedgraph ${get_chrom_sizes(get_genome(library))} clipped.bedgraph
	bedGraphToBigWig clipped.bedgraph ${get_chrom_sizes(get_genome(library))} ${library}.bw
	rm sorted.bedgraph clipped.bedgraph
	"""	
	
}


process ataqv {
	
	publishDir "${params.results}/ataqv", mode: 'rellink'
	errorStrategy 'retry'
	maxRetries 1
	memory '5 GB'
	time '3h'
	
	input:
	tuple val(library), path(bam), path(bam_index), path(peaks)
	
	output:
	tuple path("${library}.ataqv.json.gz"), path("${library}.ataqv.out")

	"""
	${IONICE} ataqv --peak-file $peaks --name ${library} --metrics-file ${library}.ataqv.json.gz --tss-file ${get_tss(get_genome(library))} ${make_excluded_regions_arg(get_genome(library))} --ignore-read-groups ${get_organism(get_genome(library))} $bam > ${library}.ataqv.out
	"""	

}


workflow {
	libraries = params.libraries.keySet()
	library_readgroup_fq1_fq2 = []

	for (library in libraries) {
		for (readgroup in library_to_readgroups(library)) {
			fastqs = library_and_readgroup_to_fastqs(library, readgroup)
			first_insert = fastqs['1']
			second_insert = fastqs['2']
			library_readgroup_fq1_fq2 << [library, readgroup, file(first_insert), file(second_insert)]
		}
	}

	library_readgroup_fq1_fq2 = Channel.from(library_readgroup_fq1_fq2)

	trimmed = trim(library_readgroup_fq1_fq2)
	library_readgroup_fq1_fq2.mix(trimmed).map({it -> [it[2], it[3]]}).flatten() | fastqc
	md_bams = bwa(trimmed).groupTuple() | merge_mapped | mark_duplicates
	peak_calling = prune(md_bams) | bamtobed | macs2
	blacklist_filter_peaks(peak_calling.peaks)
	bigwig(peak_calling.bedgraph)
	ataqv(md_bams.combine(peak_calling.peaks, by: 0))



}