echo -e "\e[34m##############################################\e[0m"
echo -e "\e[34m##### Starting variant analysis pipeline #####\n\e[0m"

echo -e "\e[34mCreating genome index...\e[0m"

bwa index genome/*.fa

echo -e "\e[34mDone\n\e[0m"

mkdir -p alignment
mkdir -p calling

for sample in normal tumor
do
    echo -e "\e[34mAligning raw seq data from each ${sample} sample to the reference genome...\e[0m"
    
    bwa mem -R "@RG\tID:exercise\tSM:${sample}" genome/*.fa \
        raw_data/${sample}_1.fq.gz raw_data/${sample}_2.fq.gz \
        > alignment/${sample}.sam
    
    echo -e "\e[34mDone\n\e[0m" 
    echo -e "\e[34mConverting the ${sample} SAM file to BAM and filling in mate information...\e[0m"
    
    samtools fixmate -O bam alignment/${sample}.sam \
             alignment/${sample}_fixmate.bam
    
    echo -e "\e[34mDone\n\e[0m"
    echo -e "\e[34mSorting the ${sample} BAM file by coordinates...\e[0m"
    
    samtools sort -O bam -o alignment/${sample}_sorted.bam \
             alignment/${sample}_fixmate.bam

    echo -e "\e[34mDone\n\e[0m"       
    echo -e "\e[34mRemoving PCR duplicates and generate a refined ${sample} BAM file...\e[0m"
   
    samtools rmdup -S alignment/${sample}_sorted.bam \
             alignment/${sample}_refined.bam
    
    echo -e "\e[34mDone\n\e[0m"
    echo -e "\e[34mIndexing the refined ${sample} BAM file...\e[0m"
    
    samtools index alignment/${sample}_refined.bam

    echo -e "\e[34mDone\n\e[0m"
    echo -e "\e[34mCalculating the most likely genotype for ${sample} sample...\e[0m"
    
    bcftools mpileup -Ou -f genome/*.fa -T intervals/*.bed\
             alignment/${sample}_refined.bam | bcftools call -vmO z -o \
             calling/${sample}_rawcalls.vcf.gz

    echo -e "\e[34mDone\n\e[0m" 
    echo -e "\e[34mCSorting and remove duplicates...\e[0m"
    
    bcftools sort -Ou calling/${sample}_rawcalls.vcf.gz | bcftools norm -d all -O z -o \
             calling/${sample}_sorted_dedup.vcf.gz

    echo -e "\e[34mDone\n\e[0m" 
    echo -e "\e[34mCreating the index for the vcf ${sample} file...\e[0m"
    
    bcftools index calling/${sample}_sorted_dedup.vcf.gz

    echo -e "\e[34mDone\n\e[0m"
    echo -e "\e[34mRetrieving only calls with high quality in ${sample} sample...\e[0m"
   
    bcftools filter -i 'QUAL>=20&&DP>=20' calling/${sample}_sorted_dedup.vcf.gz

    echo -e "\e[34mDone\n\e[0m"
done

echo -e "\e[34mMaking intersection of detected variants in the normal and tumor samples from the patient...\e[0m"

bcftools isec -i 'DP>=20' calling/normal_sorted_dedup.vcf.gz calling/tumor_sorted_dedup.vcf.gz \
         -p calling/intersection

echo -e "\e[34mDone\n\e[0m"
echo -e "\e[34m#### variant analysis pipeline completed #####\e[0m"
echo -e "\e[34m##############################################\e[0m"