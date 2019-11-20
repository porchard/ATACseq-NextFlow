# NextFlow pipeline for ATAC-seq data

## Dependencies
If you have Singularity installed, you can use the config provided here ('Singularity') to build a container with all the dependencies.

Otherwise, you'll need to have the following installed:
1. cta
2. bedtools
3. bwa
4. picardtools
5. fastqc
6. samtools
7. ataqv

I've used this pipeline with NextFlow v. 19.04.1

## Configuration
Paths to various generic files (e.g., bwa indices) must be included in the nextflow.config file -- check that file and change paths accordingly. These include:

1. Blacklist bed files for each genome
2. Chrom size files for each genome
3. BWA indices
4. TSS files (BED6 files denoting TSS positions)

You'll also need to set the params.results variable -- either in the nextflow.config file itself, or on the command line when you run the pipeline ('--results /path/to/results').

Lastly, you'll need to include information about each ATAC-seq library, including the genome that each library should be mapped to and the paths to the fastq files for each readgroup. Organize this information in a JSON file, as in library-config.json.

## Running
Once you have all of the above information, you can run the pipeline as follows (in this case, indicating the path to the results on the command line):

```bash
nextflow run -with-singularity /path/to/Singularity.simg -with-trace -with-report -with-dag -with-timeline -params-file library-config.json --results /path/to/results /path/to/main.nf
```
