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
    container 'library://porchard/default/cta:20220113'
    errorStrategy 'retry'
    maxRetries 1
    time '24h'
    maxForks 10
    memory '8 GB'

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
    container 'library://porchard/default/general:20220107'
    errorStrategy 'retry'
    maxRetries 1
    maxForks 10
    memory '8 GB'
    time '6h'

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


process multiqc {

    publishDir "${params.results}/multiqc/${before_or_after_trim}", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'
    memory '8 GB'
    time '6h'

    input:
    tuple val(before_or_after_trim), path(x)

    output:
    path('multiqc_data')
    path('multiqc_report.html')

    """
    multiqc .
    """

}


process bwa {

    publishDir "${params.results}/bwa", mode: 'rellink'
    container 'library://porchard/default/bwa:0.7.15'
    memory '50 GB'
    cpus 12
    errorStrategy 'retry'
    maxRetries 1
    time '48h'

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
    container 'library://porchard/default/general:20220107'
    errorStrategy 'retry'
    maxRetries 1
    time '5h'
    maxForks 10

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
    container 'library://porchard/default/general:20220107'
    errorStrategy 'retry'
    maxRetries 1
    time '5h'
    maxForks 15
    memory '6 GB'

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path("${library}.md.bam"), path("${library}.md.bam.bai")

    """
    java -Xmx4g -Xms4g -jar \$PICARD_JAR MarkDuplicates I=$bam O=${library}.md.bam ASSUME_SORTED=true METRICS_FILE=${library}.metrics VALIDATION_STRINGENCY=LENIENT
    samtools index ${library}.md.bam
    """
}


process prune {

    publishDir "${params.results}/prune", mode: 'rellink'
    container 'library://porchard/default/general:20220107'
    memory '3 GB'
    time '5h'
    errorStrategy 'retry'
    maxRetries 3
    maxForks 10

    input:
    tuple val(library), path(bam), path(bam_index)

    output:
    tuple val(library), path("${library}.pruned.bam")

    """
    ${IONICE} samtools view -h -b -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 $bam ${AUTOSOMAL_REFERENCES[get_genome(library)].join(' ')} > ${library}.pruned.bam 
    """

}


process bamtobed {

    container 'library://porchard/default/general:20220107'
    time '4h'
    maxForks 10

    input:
    tuple val(library), path(bam)

    output:
    tuple val(library), path("${library}.bed")

    """
    bedtools bamtobed -i $bam > ${library}.bed
    """

}


process macs2_broad {

    publishDir "${params.results}/macs2/broad", mode: 'rellink'
    container 'library://porchard/default/general:20220107'
    time '5h'
    memory '8 GB'

    input:
    tuple val(library), path(bed)

    output:
    tuple val(library), path("${library}_peaks.broadPeak"), emit: peaks
    tuple val(library), path("${library}_treat_pileup.bdg"), emit: bedgraph

    """
    macs2 callpeak -t $bed --outdir . --SPMR -f BED -n $library -g ${get_macs2_genome_size(get_genome(library))} --nomodel --shift -100 --seed 762873 --extsize 200 -B --broad --keep-dup all
    """

}


process macs2_narrow {

    publishDir "${params.results}/macs2/narrow", mode: 'rellink'
    container 'library://porchard/default/general:20220107'
    time '5h'
    memory '8 GB'

    input:
    tuple val(library), path(bed)

    output:
    tuple val(library), path("${library}_peaks.narrowPeak"), emit: peaks
    tuple val(library), path("${library}_treat_pileup.bdg"), emit: bedgraph
    tuple val(library), path("${library}_summits.bed"), emit: summits

    """
    macs2 callpeak -t $bed --qvalue 0.05 --outdir . --SPMR -f BED -n $library -g ${get_macs2_genome_size(get_genome(library))} --nomodel --shift -37 --seed 762873 --extsize 73 -B --keep-dup all --call-summits
    """

}


process extend_summits {

    publishDir "${params.results}/macs2/narrow", mode: 'rellink'
    container 'library://porchard/default/general:20220107'
    time '1h'
    memory '8 GB'

    input:
    tuple val(library), path(bed)

    output:
    tuple val(library), path("${library}_summits_extended.bed")

    """
    extend-summits.py --extension 150 $bed ${get_chrom_sizes(get_genome(library))} > ${library}_summits_extended.bed
    """

}


process blacklist_filter {

    publishDir "${params.results}/macs2/${macs2_subdir}", mode: 'rellink'
    container 'library://porchard/default/general:20220107'
    time '1h'

    input:
    tuple val(library), path(bed)

    output:
    path("$outfile")

    when:
    has_blacklist(get_genome(library))

    script:
    outfile = bed.getName() + '.noblacklist'
    macs2_subdir = outfile.contains('broad') ? 'broad' : 'narrow'

    """
    bedtools intersect -a $bed -b ${get_blacklists(get_genome(library)).join(' ')} -v > $outfile
    """

}


process bigwig {

    publishDir "${params.results}/bigwig", mode: 'rellink'
    container 'library://porchard/default/general:20220107'
    time '5h'
    memory '8 GB'

    input:
    tuple val(library), path(bedgraph)

    output:
    tuple val(genome), path("${library}.bw")

    script:
    genome = get_genome(library)

    """
    LC_COLLATE=C sort -k1,1 -k2n,2 $bedgraph > sorted.bedgraph
    bedClip sorted.bedgraph ${get_chrom_sizes(genome)} clipped.bedgraph
    bedGraphToBigWig clipped.bedgraph ${get_chrom_sizes(genome)} ${library}.bw
    rm sorted.bedgraph clipped.bedgraph
    """

}


process plot_signal_at_tss {

    publishDir "${params.results}/bigwig/plot", mode: 'rellink', overwrite: true
    container 'library://porchard/default/general:20220107'
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    tag "${genome}"

    input:
    tuple val(genome), path(bw)

    output:
    path("*.png") optional true

    """
    plot-signal-at-tss.py --genes ${params.plot_signal_at_genes.join(' ')} --tss-file ${get_tss(genome)} --bigwigs ${bw.join(' ')} --prefix ${genome}.
    """

}


process ataqv {

    publishDir "${params.results}/ataqv", mode: 'rellink'
    container 'library://porchard/default/ataqv:1.3.0'
    errorStrategy 'retry'
    maxRetries 1
    memory '5 GB'
    time '3h'

    input:
    tuple val(library), path(bam), path(bam_index), path(peaks)

    output:
    tuple val(genome), path("${library}.ataqv.json.gz"), emit: for_viewer
    path("${library}.ataqv.out")

    script:
    genome = get_genome(library)

    """
    export TERM=xterm-256color && ataqv --peak-file $peaks --name ${library} --metrics-file ${library}.ataqv.json.gz --tss-file ${get_tss(genome)} ${make_excluded_regions_arg(genome)} --ignore-read-groups ${get_organism(genome)} $bam > ${library}.ataqv.out
    """

}


process ataqv_viewer {

    publishDir "${params.results}/ataqv-viewer", mode: 'rellink'
    container 'library://porchard/default/ataqv:1.3.0'
    errorStrategy 'retry'
    maxRetries 1
    memory { 10.GB * task.attempt }
    time '1h'
    tag "${genome}"

    input:
    tuple val(genome), path(json)

    output:
    path("ataqv-viewer-${genome}")

    """
    export TERM=xterm-256color && mkarv ataqv-viewer-${genome} ${json.join(' ')}
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
    (library_readgroup_fq1_fq2.mix(trimmed).map({it -> [it[2], it[3]]}).flatten() | fastqc).map({it -> [it.getName().contains('trimmed') ? 'after-trim' : 'before_trim', it]}).groupTuple() | multiqc
    md_bams = bwa(trimmed).groupTuple() | merge_mapped | mark_duplicates
    beds = prune(md_bams) | bamtobed

    narrow_peak_calling = macs2_narrow(beds)
    extended_summits = extend_summits(narrow_peak_calling.summits)
    broad_peak_calling = macs2_broad(beds)

    broad_peak_calling.peaks.mix(narrow_peak_calling.peaks).mix(extended_summits) | blacklist_filter

    bigwig(broad_peak_calling.bedgraph).groupTuple() | plot_signal_at_tss
    ataqv(md_bams.combine(broad_peak_calling.peaks, by: 0)).for_viewer.groupTuple() | ataqv_viewer

}