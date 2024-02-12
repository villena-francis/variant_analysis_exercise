#Create genome index
bwa index genome/*.fa

mkdir -p alignment

for sample in normal tumor
do
    #Aligment of raw seq data from each sample to the reference genome
    bwa mem -R "@RG\tID:exercise\tSM:${sample}" genome/*.fa \
        raw_data/${sample}_1.fq.gz raw_data/${sample}_2.fq.gz \
        > alignment/${sample}.sam
    #Convert the SAM file to BAM and fill in mate information
    samtools fixmate -O bam alignment/${sample}.sam \
             alignment/${sample}_fixmate.bam
    #Sort the BAM file by coordinates
    samtools sort -O bam -o alignment/${sample}_sorted.bam \
             alignment/${sample}_fixmate.bam
    #Remove PCR duplicates and generate a refined BAM file
    samtools rmdup -S alignment/${sample}_sorted.bam \
             alignment/${sample}_refined.bam
    #Index the refined BAM file
    samtools index alignment/${sample}_refined.bam