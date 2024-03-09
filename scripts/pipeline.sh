blue="\e[34m"
endblue="\e[0m"

echo -e "${blue}##############################################${endblue}"
echo -e "${blue}##### Starting variant analysis pipeline #####\n${endblue}"

echo -e "${blue}Creating genome index...${endblue}"

bwa index genome/*.fa

echo -e "${blue}Done\n${endblue}"

mkdir -p alignment
mkdir -p calling

for sample in normal tumor
do
    echo -e "${blue}Aligning raw seq data from each ${sample} sample to the reference genome...${endblue}"
    
    bwa mem -R "@RG\tID:exercise\tSM:${sample}" genome/*.fa \
        raw_data/${sample}_1.fq.gz raw_data/${sample}_2.fq.gz \
        > alignment/${sample}.sam
    
    echo -e "${blue}Done\n${endblue}" 
    echo -e "${blue}Converting the ${sample} SAM file to BAM and filling in mate information...${endblue}"
    
    samtools fixmate -O bam alignment/${sample}.sam \
             alignment/${sample}_fixmate.bam
    
    echo -e "${blue}Done\n${endblue}"
    echo -e "${blue}Sorting the ${sample} BAM file by coordinates...${endblue}"
    
    samtools sort -O bam -o alignment/${sample}_sorted.bam \
             alignment/${sample}_fixmate.bam

    echo -e "${blue}Done\n${endblue}"       
    echo -e "${blue}Removing PCR duplicates and generate a refined ${sample} BAM file...${endblue}"
   
    samtools rmdup -S alignment/${sample}_sorted.bam \
             alignment/${sample}_refined.bam
    
    echo -e "${blue}Done\n${endblue}"
    echo -e "${blue}Indexing the refined ${sample} BAM file...${endblue}"
    
    samtools index alignment/${sample}_refined.bam

    echo -e "${blue}Done\n${endblue}"
    echo -e "${blue}Calculating the most likely genotype for ${sample} sample...${endblue}"
    
    bcftools mpileup -Ou -f genome/*.fa -T intervals/*.bed\
             alignment/${sample}_refined.bam | bcftools call -vmO z -o \
             calling/${sample}_rawcalls.vcf.gz

    echo -e "${blue}Done\n${endblue}" 
    echo -e "${blue}CSorting and remove duplicates...${endblue}"
    
    bcftools sort -Ou calling/${sample}_rawcalls.vcf.gz | bcftools norm -d all -O z -o \
             calling/${sample}_sorted_dedup.vcf.gz

    echo -e "${blue}Done\n${endblue}" 
    echo -e "${blue}Creating the index for the vcf ${sample} file...${endblue}"
    
    bcftools index calling/${sample}_sorted_dedup.vcf.gz

    echo -e "${blue}Done\n${endblue}"
    echo -e "${blue}Retrieving only calls with high quality in ${sample} sample...${endblue}"
   
    bcftools filter -i 'QUAL>=20&&DP>=20' calling/${sample}_sorted_dedup.vcf.gz

    echo -e "${blue}Done\n${endblue}"
done

echo -e "${blue}Making intersection of detected variants in the normal and tumor samples from the patient...${endblue}"

bcftools isec -i 'QUAL>=20&&DP>=20' calling/normal_sorted_dedup.vcf.gz calling/tumor_sorted_dedup.vcf.gz \
         -p calling/intersection

echo -e "${blue}Done\n${endblue}"
echo -e "${blue}#### variant analysis pipeline completed #####${endblue}"
echo -e "${blue}##############################################${endblue}"