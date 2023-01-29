params.reference = "hg19"
params.fastq_dir = "data/fastq"
params.bam_dir = "data/bam"
params.peak_dir = "results/peaks"

process alignReads {
  publishDir params.bam_dir
  input:
    file("${params.fastq_dir}/*_R1.fastq.gz"), file("${params.fastq_dir}/*_R2.fastq.gz")
  output:
    bam = "${params.bam_dir}/${sampleID}.bam"
  script:
    """
    bwa mem -t 8 ${params.reference}.fa ${input1} ${input2} | \
    samtools sort -o ${bam}
    """
}

process callPeaks {
  publishDir params.peak_dir
  input:
    bam from alignReads
  output:
    bed = "${params.peak_dir}/${sampleID}.bed"
  script:
    """
    macs2 callpeak -t ${bam} -f BAM -g hs -n ${sampleID} --outdir ${params.peak_dir}
    """
}

process annotatePeaks {
  input:
    bed from callPeaks
  output:
    annotated = "${params.peak_dir}/${sampleID}_annotated.bed"
  script:
    """
    bedtools intersect -a ${input} -b hg19_refseq.bed -wa > ${annotated}
    """
}

process runChIPseq {
  input:
    set sampleID from params.fastq_dir
  output:
    file("${params.peak_dir}/${sampleID}_annotated.bed")
  into {
    sampleID = fileBaseName(input)
  }
  call annotatePeaks
}

runChIPseq.maxForks = 8

