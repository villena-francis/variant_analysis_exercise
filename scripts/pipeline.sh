#Create genome index
bwa index genome/*.fa

mkdir -p alignment
mkdir -p calling

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
    #Calculate the most likely genotype
    bcftools mpileup -Ou -f genome/*.fa \
             alignment/${sample}_refined.bam | bcftools call -vmO z -o \
             calling/${sample}_rawcalls.vcf.gz
    #Sort and remove duplicates
    bcftools sort -Ou calling/${sample}_rawcalls.vcf.gz | bcftools norm -d all -O z -o \
             calling/${sample}_sorted_dedup.vcf.gz
    #Create the index for the vcf file
    bcftools index calling/${sample}_sorted_dedup.vcf.gz
    #Retrieve only calls with high quality
    bcftools filter -i 'QUAL>=20&&DP>=20' calling/${sample}_sorted_dedup.vcf.gz
done

#Intersection of detected variants in the two samples from the patient
bcftools isec -i 'DP>=20' calling/normal_sorted_dedup.vcf.gz calling/tumor_sorted_dedup.vcf.gz \
         -p calling/intersection