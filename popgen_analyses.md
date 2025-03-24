## Scripts for the population genomic analyses of the European polecat 

# Raw read processing of samples ##
Scripts used to align reads and downsample bams for input to GATK
```

#!/bin/bash -e

source package fa33234e-dceb-4a58-9a78-7bcf9809edd7 #BWA aligner
source package 638df626-d658-40aa-80e5-14a275b7464b #samtools

## align reads to reference ##
SAMPLE=$1
REF=$2
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/data/raw/polecat_data
#bwa index -a bwtsw $1

while read -r i; do
srun bwa mem -t 16 -R "@RG\tID:${SAMPLE}\tPL:illumina\tSM:${SAMPLE}" ${REF} ${SAMPLE}_1P.fastq.gz ${SAMPLE}_2P.fastq.gz |
samtools sort - -o ${SAMPLE}_sortedmputorius.bam
#done < samples.txt

source package 638df626-d658-40aa-80e5-14a275b7464b #samtools

SAMPLE=$1

## convert sam file to bam file then sort ##
#srun samtools view -b ${SAMPLE}.sam -o ${SAMPLE}.bam  
samtools sort ${SAMPLE}.bam -o ${SAMPLE}_sorted.bam -@8

## markduplicate reads ##
source picardtools-2.1.1


wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/sortedbamfiles/
outd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123

while read -r i; do
java -jar /ei/software/testing/picardtools/2.1.1/x86_64/bin/picard.jar MarkDuplicates I=${wd}/${i}_sorted.bam O=${outd}/${i}_mkd.bam M=${outd}/$>
done <samples.txt



## calculate coverage of individuals ##
source package 638df626-d658-40aa-80e5-14a275b7464b #samtools

FILE=$1

#samtools coverage ${FILE }-o ${FILE}coverage  

## downsample higher coverage samples ##

#calculate proportion to downsample using script calc_cov.R

TARGET_DIR=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123
# List of files
bam_files=(${TARGET_DIR}/*_mkd.bam)
coverage=${TARGET_DIR}/proportion_to_downsample.txt

# Get file to be processed by this task
this_bam=${bam_files[$((SLURM_ARRAY_TASK_ID-1))]}
#this_bam=${this_bam_gz%.gz}

echo "Processing file: ${this_bam} on $HOSTNAME"

# Decompress the BAM file
#gunzip -c $this_bam_gz > $this_bam

# Proportion to downsample for array job
PROP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $coverage | awk '{print $2}')

echo "At this proportion: $PROP"

# Downsample the BAM file
samtools view -s $PROP -b $this_bam > ${this_bam%.bam}_downsample.bam

# Remove the decompressed BAM file to save space
gzip $this_bam

```
## GATK variant calling
```
## first step is to run gatk haplotype caller on samples in order to get an estimate which can be then used to recalibrate bam files and rerun gatk haplotype caller ##

#define variables
samplename=$1
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.vcf
REF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/bcftools_231123/MusPutFur1.0_bionano.fasta
RECBAM=${wd}/$samplename\_recal.bam
GVCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.g.vcf.gz

source package a2643e0d-b838-4ca0-a31a-6f4282b5ad56  #gatk
# first run of haplotype caller. each sample is run seperately and a seperate VCF will be generated #

srun gatk HaplotypeCaller -R ${REF} -I ${wd}/${samplename}_mkd.bam -O $VCF

# PCRfree model #

srun gatk HaplotypeCaller -R ${REF} --pcr-indel-model NONE -I ${wd}/${samplename}_mkd.bam -O $VCF

# Filter variant set per sample

samplename=$1
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123/
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.vcf
REF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/bcftools_231123/MusPutFur1.0_bionano.fasta
FILT=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\_filtered.vcf
source package a2643e0d-b838-4ca0-a31a-6f4282b5ad56  #gatk


srun gatk VariantFiltration  -R $REF -V $VCF -O $FILT --filter-name "LowMQ" -filter "MQ<40.00" --filter-name 'LowQual' -filter "QUAL<30" --filter-name 'LowCov' -filter "DP<7"  --filter-name 'HighCov'  -filter "DP>20" 


#Can then run the base recalibrator per sample to generate the recalibration table 

#define variables
samplename=$1
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/markdups_bam_91123
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.vcf
REF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/bcftools_231123/MusPutFur1.0_bionano.fasta
FILT=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\_filtered.vcf
RECTABLE=$samplename\recal_data.table

source package a2643e0d-b838-4ca0-a31a-6f4282b5ad56  #gatk
srun gatk  BaseRecalibrator -R ${REF} -I ${wd}/${samplename}_mkd.bam --known-sites ${FILT} -O ${RECTABLE} 

# Apply the recalibration table to the BAM files 

srun gatk  ApplyBQSR -R ${REF} -I ${wd}/${samplename}_mkd.bam  --bqsr-recal-file ${RECTABLE} -O ${RECBAM}



# second run of haplotype caller after the base recalibration has been applied to bam file on known sites #
srun gatk HaplotypeCaller  -R ${REF} -I $RECBAM -O $GVCF -ERC GVCF

###ERC is Estimate the confidence that no SNP exists at the site by contrasting all reads with the REF base vs. all reads with any non-reference base.


# second run of haplotype caller using pcr-free model after the base recalibration has been applied to bam file on known sites #
srun gatk HaplotypeCaller  -R ${REF} --pcr-indel-model NONE -I $RECBAM -O $GVCF -ERC GVCF



# Consolidate GVCFs  - This joint calls genotypes and first requires building a datastore to consolidate all GVCFS using GenomicsDBImport. #

samplename=$1
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/
workspace=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/all_samplesGenomicsDB #have to name the workspace for database to be made
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.vcf
REF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/bcftools_231123/MusPutFur1.0_bionano.fasta
RECBAM=${wd}/$samplename\_recal.bam
GVCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.g.vcf.gz
scaffolds=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/scaffolds.intervals #has to have the suffix .intervals or .list
map=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/gvcf_sample_list

##here the ${map} file is laid out as
sample1	sample1.gvcf
sample2	sample2.gvcf

##here the ${intervals} file is a list of all scaffolds/chromosomes in ref genome


source package a2643e0d-b838-4ca0-a31a-6f4282b5ad56  #gatk 4.1.4.0
srun gatk --java-options "-Xmx900g" GenomicsDBImport -R $REF  --intervals ${scaffolds} --genomicsdb-workspace-path ${workspace} --sample-name-map ${map}  --tmp-dir ${wd} -max-num-intervals-to-import-in

#here, Xmx tells the software how much space it has to to build the database - make sure to tell java that it has less space than it does due to how greedy it is.


# Then we can do a joint genotype call for snps and indels using all per sample gvcfs #

samplename=$1
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk
workspace=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/all_samples_GenomicsDB #have to name the workspace for database to be made
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.vcf
REF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/bcftools_231123/MusPutFur1.0_bionano.fasta
RECBAM=${wd}/$samplename\_recal.bam
GVCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/$samplename\.g.vcf.gz
scaffolds=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/scaffolds.intervals #has to have the suffix .intervals or .list
map=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/gvcf_sample_list

export SINGULARITYENV_TILEDB_DISABLE_FILE_LOCKING=1 #this is required to disable file locking in gatk v 4.1.4

source package a2643e0d-b838-4ca0-a31a-6f4282b5ad56  #gatk

gatk GenotypeGVCFs -R $REF  -V gendb://all_samples_GenomicsDB  -O all_samples_genotyped.vcf ##this took 7 days to run

#select just snps
gatk --java-options "-Xmx550g" SelectVariants  -R $REF  -V all_samples_genotyped.vcf -O all_samples_genotyped_downsample_snps.vcf -select-type SNP


```
## Quality filtering and selecting for biallelic SNPs only
```
source package 638df626-d658-40aa-80e5-14a275b7464b # Samtools - 1.15.1

#Filter SNPs with a genotype quality under 20

bcftools view -i 'QUAL>=20' ${VCF} -O z -o  all_samples_genotyped_downsample_snps_Q20.vcf

##filter for AC=>3
bcftools view --min-ac 3:minor ${VCF}  -O z -o  all_samples_genotyped_downsample_snps_Q20_AC3.vcf

##select only biallelic SNPs
#bcftools view -m 2 -M 2 ${VCF} -Oz -o all_samples_genotyped_downsample_snps_Q20_AC3_BI.vcf
```
## Pruning
```
source package 638df626-d658-40aa-80e5-14a275b7464b # Samtools - 1.15.1


VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/all_samples_genotyped_downsample_snps_Q20_AC3_BI.vcf
OUT=${VCF%.vcf}_pruned.vcf
bcftools +prune ${VCF} -m 0.6 -w 50kb -Oz -o  ${OUT}
```

# Relatedness
```
source vcftools-0.1.16
vcftools --vcf ${wd}/all_samples_genotyped_downsample_snps_Q20_AC3_BI_pruned.vcf --relatedness2
```

# PCA
```
source plink-1.90.beta6.6 
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_pruned.vcf 

plink --vcf $VCF  --allow-extra-chr --const-fid  --make-bed  --maf 0.1 --pca --out all_downsampled
```

# Isolation by distance 
see mantel_test.R

## Genomic diversity estimates

#ROH
```
source plink-1.90.beta6.6 
# calculate roh for each individual
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_pruned.vcf 

plink --vcf ${VCF} --het  --geno 0.5 --maf 0.1 --homozyg --homozyg-density 50 --homozyg-gap 1000 --homozyg-kb 300 --homozyg-snp 50 --homozyg-window-het 4 --homozyg-window-missing 5 --homozyg-window-snp 50 --homozyg-window-threshold 0.05  --aec --double-id

```
#Genomic diversity
```
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_pruned.vcf 

INDS=($(for i in *_pop
do
       echo $(basename ${i})
done))

# calculate pi within each population

for IND in ${INDS[@]}
do
    OUT=mput_${IND}
	vcftools --vcf $VCF --maf 0.1 --window-pi 20000000 --keep ${IND} --out ${OUT}_20mb
done
```
#Heterozygosity
```
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_pruned.vcf 
OUT_HET=${VCF%._het}

vcftools --vcf ${VCF}  --het --out $OUT_HET
```


#FST
```
source vcftools-0.1.16
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered
VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_pruned.vcf 
 
vcftools --vcf $VCF --weir-fst-pop welsh_pop  --weir-fst-pop english_pop  --out musput_pairwisefst_welsh_english
vcftools --vcf $VCF --weir-fst-pop welsh_pop  --weir-fst-pop euro_pop  --out musput_pairwisefst_welsh_euro
vcftools --vcf $VCF --weir-fst-pop english_pop  --weir-fst-pop euro_pop  --out musput_pairwisefst_euro_english
```

#Private SNPs
```
source package 638df626-d658-40aa-80e5-14a275b7464b # Samtools - 1.15.1
VCF1=all_samples_genotyped_downsample_snps_Q20_AC3_BI_unrelated_englishonly.vcf_mono.gz 
VCF2=all_samples_genotyped_downsample_snps_Q20_AC3_BI_unrelated_euroonly.vcf_mono.gz 
VCF3=all_samples_genotyped_downsample_snps_Q20_AC3_BI_unrelated_welshonly.vcf_mono.gz 

bcftools isec ${VCF1} ${VCF2} -p euro_english_dir

bcftools isec ${VCF1} ${VCF3} -p euro_welsh_dir

#bcftools isec ${VCF2}.gz ${VCF3}.gz -p euro_welsh_dir

# Define the data
private_to_welsh <- c(66198,512671,0)
private_to_english <- c(257909,0,580562)
private_to_euro <- c(0,498434,1791873)

# Shared values
shared_by_both_12 <- 2268113
shared_by_both_13 <- 1821640
shared_by_both_23 <- 1945460
shared_all <- 1782650

# Areas
area1 <- private_to_welsh[1] + private_to_welsh[2] + shared_by_both_12 + shared_by_both_13 + shared_all
area2 <- private_to_english[1] +  private_to_english[3] + shared_by_both_12 + shared_by_both_23 + shared_all
area3 <- private_to_euro[2] +  private_to_euro[3] + shared_by_both_13 + shared_by_both_23 + shared_all
print(area1)
# Overlaps
n12 <- shared_by_both_12 + shared_all
n13 <- shared_by_both_13 + shared_all
n23 <- shared_by_both_23 + shared_all
n123 <- shared_all
```
# SnpEff - annotation of called variants 
```
source package ae02bb7b-09d8-4325-a480-6d7c17b4e8af # snpeff 5.1

# build snpeff database for domestic ferret #
# need to have a directory (here named MusPutFur1.0_bionano) where there is the genic information in genes.gff file and gene sequences in a sequence.fa file (this needs to be indexed and a .fai file in the same directory)

snpEff build -c ./M_putorius.cfg -gff3 -v MusPutFur1.0_bionano

# run snpeff on population variants #
welsh_VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_welshonly.vcf_mono.gz
english_VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_englishonly.vcf_mono.gz
euro_VCF=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/gatk/bcftools_filtered/downsampled_data_280724/rerun_paper_jan25/all_samples_genotyped_downsample_snps_Q20_AC3_BI_euroonly.vcf_mono.gz

english_out=${english_VCF%.gz}_lofsnpeff.vcf
welsh_out=${welsh_VCF%.gz}_lofsnpeff.vcf
euro_out=${euro_VCF%.gz}_lofsnpeff.vcf


snpEff -c M_putorius.cfg  MusPutFur1.0_bionano -lof ${english_VCF} > ${english_out};
snpEff -c M_putorius.cfg  MusPutFur1.0_bionano -lof ${welsh_VCF} > ${welsh_out};
snpEff -c M_putorius.cfg  MusPutFur1.0_bionano -lof ${euro_VCF} > ${euro_out};

# run snpsift - example only for welsh population #
welsh_VCF=welsh_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf


#java -jar /tgac/software/testing/bin/core/../..//snpeff/4.3/x86_64/bin/SnpSift.jar  filter "ANN[0].EFFECT has 'missense_variant'"  ${welsh_VCF} > missense_${welsh_VCF}

#java -jar /tgac/software/testing/bin/core/../..//snpeff/4.3/x86_64/bin/SnpSift.jar  filter "ANN[0].EFFECT has 'nonsense_variant'"  ${welsh_VCF} > nonsense_${welsh_VCF}

#java -jar /tgac/software/testing/bin/core/../..//snpeff/4.3/x86_64/bin/SnpSift.jar  filter "( ANN[*].EFFECT = 'nonsynonymous_variant' )" ${welsh_VCF} > nonsyn_${welsh_VCF}

#java -jar /tgac/software/testing/bin/core/../..//snpeff/4.3/x86_64/bin/SnpSift.jar  filter "( ANN[*].EFFECT = 'synonymous_variant' )" ${welsh_VCF} > syn_${welsh_VCF}

java -jar /tgac/software/testing/bin/core/../..//snpeff/4.3/x86_64/bin/SnpSift.jar  filter "( exists LOF[*].GENE )" ${welsh_VCF} > lof_${welsh_VCF}

# Query the VCF file to get genotype information for each sample

source package 638df626-d658-40aa-80e5-14a275b7464b #samtools 1.15.1

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' lof_${euro_VCF} > euro_genotype_lof_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' lof_${welsh_VCF} > welsh_genotype_lof_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' lof_${english_VCF} > english_genotype_lof_info.txt

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' missense_${euro_VCF} > euro_missense_genotype_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' missense_${welsh_VCF} > welsh_missense_genotype_info.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[\t%SAMPLE=%GT]\n' missense_${english_VCF} > english_missense_genotype_info.txt

# extract whether individuals are HOM or HET for each mutation
#chrom, pos, ref, alt store the chromosome, position, reference allele, and alternative allele for each variant.
#for(i=5; i<=NF; i++) loops over each sample's genotype field (which starts from the 5th column).
#split($i, sample_info, "=") splits each sample's entry into sample and genotype parts.
#The if condition checks the genotype type (heterozygous 0/1 or homozygous derived 1/1) and categorizes accordingly.
#It prints sample name, SNP identifier (chrom:pos), and genotype status (HET or HOM).

awk 'BEGIN {FS="\t"; OFS="\t"}
{
    chrom=$1; pos=$2; ref=$3; alt=$4;
    for(i=5; i<=NF; i++) {
        split($i, sample_info, "=");
        sample = sample_info[1];
        genotype = sample_info[2];
        if (genotype == "0/1" || genotype == "1/0") {
            print sample, chrom ":" pos, "HET";
        } else if (genotype == "1/1") {
            print sample, chrom ":" pos, "HOM";
        }
    }
}' genotype_info.txt > categorized_genotypes.txt

#Sample and Genotype Parsing:
#sample = $1: Extracts the sample name from the first column.
#genotype = $3: Extracts the genotype status (HET or HOM) from the third column.
#Count Initialization:
#if (!(sample in het_count)) { ... }: Initializes the count for het and hom mutations if this is the first time the sample is encountered.
#Counting:
#het_count[sample]++: Increments the count for het mutations for the sample.
#hom_count[sample]++: Increments the count for hom mutations for the sample.

awk '
{
    sample = $1;
    genotype = $3;

    # Initialize counts if sample is encountered for the first time
    if (!(sample in het_count)) {
        het_count[sample] = 0;
        hom_count[sample] = 0;
    }

    # Increment the count based on genotype
    if (genotype == "HET") {
        het_count[sample]++;
    } else if (genotype == "HOM") {
        hom_count[sample]++;
    }
}
END {
    # Print the results
    print "Sample\tHET_Count\tHOM_Count";
    for (sample in het_count) {
        print sample "\t" het_count[sample] "\t" hom_count[sample];
    }
}
```
# HAL Snps - to determine the ancestral allele and the derived allele for variants
```
source package c8a0203d-707a-44e8-acdd-2c6ec17f8d7d #hal 2.2

##running halliftover

#mammalhal=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/masked_cactus_genomes/branchlength_cactus/branchlength_cactus_aln.hal
mammalhal=/ei/projects/3/3383c516-0548-45b8-89b4-74e075d066bd/scratch/masked_cactus_genomes/branchlength_cactus/branchlength_cactus_aln.hal

halSnps ${mammalhal} ferret putorius,eversmanii,weasel --tsv  ferret_polecat_eversmanni_nivalis_snps_160225.tsv

# intersect hal_VCF with annotated snpeff vcfs #
wd=/ei/projects/6/61ea8b33-6b11-4e22-8ca7-67abbd719d07/scratch/mustelids/polecat_data/SnpEff
welshlofVCF=${wd}/lof_welsh_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf
englishlofVCF=${wd}/lof_english_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf
eurolofVCF=${wd}/lof_euro_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf
welshmissenseVCF=${wd}/missense_welsh_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf
englishmissenseVCF=${wd}/missense_english_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf
euromissenseVCF=${wd}/missense_euro_snps_Q20_AC3_BI_unrelated_lofsnpeff.vcf
VCF_HAL=${wd}/ferret_polecat_eversmanni_nivalis_snps_160225.bed 

bedtools intersect -a ${VCF_HAL} -b ${welshlofVCF}  -wa -wb > lof_welsh_variants.bed
bedtools intersect -a ${VCF_HAL} -b ${englishlofVCF}  -wa -wb  > lof_english_variants.bed
bedtools intersect -a ${VCF_HAL} -b ${eurolofVCF}  -wa -wb > lof_euro_variants.bed

bedtools intersect -a ${VCF_HAL} -b ${welshmissenseVCF}  -wa -wb > missense_welsh_variants.bed
bedtools intersect -a ${VCF_HAL} -b ${englishmissenseVCF}  -wa -wb > missense_english_variants.bed
bedtools intersect -a ${VCF_HAL} -b ${euromissenseVCF}  -wa -wb > missense_euro_variants.bed

# determine ancestral or derived allele in polecat #
see ancestral_allele.py
```
# GONE - estimating recent historic effective population size 
```
# input files were prepped for input to exclude mono allelic sites and to only include the largest scaffolds, that don't over lap with sex chromosomes and that had over 10,000 SNPs #

# filter to contain only the required scaffolds 
source package 638df626-d658-40aa-80e5-14a275b7464b
VCF_english=all_samples_genotyped_downsample_snps_Q20_AC3_BI_englishhigh.vcf
out_english=${VCF_english%.gz}_100GONE
bcftools view ${VCF_english} --regions Super-Scaffold_1,Super-Scaffold_2,Super-Scaffold_3,Super-Scaffold_4,Super-Scaffold_5,Super-Scaffold_6,Super-Scaffold_7,Super-Scaffold_8,Super-Scaffold_9,Super-Scaffold_10,Super-Scaffold_11,Super-Scaffold_12,Super-Scaffold_13,Super-Scaffold_14,Super-Scaffold_15,Super-Scaffold_16,Super-Scaffold_17,Super-Scaffold_18,Super-Scaffold_19,Super-Scaffold_20,Super-Scaffold_21,Super-Scaffold_22,Super-Scaffold_23,Super-Scaffold_24,Super-Scaffold_25,Super-Scaffold_26,Super-Scaffold_27,Super-Scaffold_28,Super-Scaffold_29,Super-Scaffold_30,Super-Scaffold_31,Super-Scaffold_32,Super-Scaffold_33,Super-Scaffold_35,Super-Scaffold_36,Super-Scaffold_37,Super-Scaffold_38,Super-Scaffold_39,Super-Scaffold_40,Super-Scaffold_41,Super-Scaffold_42,Super-Scaffold_43,Super-Scaffold_44,Super-Scaffold_46,Super-Scaffold_52,Super-Scaffold_55,Super-Scaffold_56,Super-Scaffold_57,Super-Scaffold_58,Super-Scaffold_59,Super-Scaffold_60,Super-Scaffold_61,Super-Scaffold_62,Super-Scaffold_63,Super-Scaffold_65,Super-Scaffold_66,Super-Scaffold_67,Super-Scaffold_69,Super-Scaffold_70,Super-Scaffold_71,Super-Scaffold_72,Super-Scaffold_74,Super-Scaffold_76,Super-Scaffold_77,Super-Scaffold_78,Super-Scaffold_79,Super-Scaffold_80,Super-Scaffold_81,Super-Scaffold_82,Super-Scaffold_83,Super-Scaffold_85,Super-Scaffold_87,Super-Scaffold_91,Super-Scaffold_92,Super-Scaffold_95,Super-Scaffold_96,Super-Scaffold_98,Super-Scaffold_99,Super-Scaffold_102,Super-Scaffold_105,Super-Scaffold_111,Super-Scaffold_118,Super-Scaffold_121,Super-Scaffold_126,Super-Scaffold_128,Super-Scaffold_129,Super-Scaffold_132,Super-Scaffold_139,Super-Scaffold_141,Super-Scaffold_146,Super-Scaffold_148,Super-Scaffold_149,Super-Scaffold_154,Super-Scaffold_156,Super-Scaffold_163,Super-Scaffold_165,Super-Scaffold_167,Super-Scaffold_189,Super-Scaffold_190,Super-Scaffold_197 -o ${out_english}

##remove monomorphic snps
#VCF1=$1
#OUT1=${VCF1%.}_mono

bcftools filter -e 'AC=0'  ${VCF1} -o ${OUT1}

#prep files and directories for running
#!/bin/bash

IND=($(for i in {0..99}
do
	echo "$i"
done))

# loop to make 100 files for each population, copy all the required files (input and gone files) and run gone within each of the population replicates
for INDS in ${IND[@]}
do
	pop=english_high
	mkdir ${pop}_${INDS}
	mkdir ${pop}_${INDS}/PROGRAMMES
	cp -r PROGRAMMES/* ${pop}_${INDS}/PROGRAMMES
	cp ${pop}.map ${pop}_${INDS}/
	cp ${pop}.ped ${pop}_${INDS}/
	cp INPUT_PARAMETERS_FILE ${pop}_${INDS}/
	cp script_GONE.sh ${pop}_${INDS}/
done
sbatch GONE2.sh ${pop}

# run as an array #

# Define the directories
directories=()

#loop through from english_0 to english_99 and add each directory path

for ((i=0; i<99; i++)); do
        directories+=("$GONEwd/english_high_$i");
done

#calculate index of directory to process based on slurm id
dir_index=$((SLURM_ARRAY_TASK_ID - 1))
dir="${directories[$dir_index]}"

# Check if the directory exists
if [ -d "$dir" ]; then
    echo "Running script in directory: $dir"
    cd "$dir" || exit
    ${dir}/script_GONE.sh english_high # Replace "your_script.sh" with the name of your script
else
    echo "Directory not found: $dir"
fi

## script_GONE.sh ##
#!/bin/bash
#script_GONE.sh

########################################################

#Check number of arguments
if [ $# -ne 1 ]  
then
	echo "Usage: $0 <FILE>" 
	exit 1
fi

### Set arguments

FILE=$1  ### Data file name for files .ped and .map

### Take input parameters from file INPUT_PARAMETERS_FILE

source INPUT_PARAMETERS_FILE

###################### FILES NEEDED ########################

### data.ped
### data.map

### EXECUTABLES FILES NEEDED IN DIRECTORY PROGRAMMES:

### MANAGE_CHROMOSOMES2
### LD_SNP_REAL3
### SUMM_REP_CHROM3
### GONE (needs gcc/7.2.0)
### GONEaverages
### GONEparallel.sh

################### Remove previous output files ##################

if [ -f "OUTPUT_$FILE" ]
then
rm OUTPUT_$FILE
fi

if [ -f "Ne_$FILE" ]
then
rm Ne_$FILE
fi

################### Create temporary directory ##################

if [ -d "TEMPORARY_FILES" ]
then
rm -r TEMPORARY_FILES
fi

mkdir TEMPORARY_FILES

################### Obtain sample size, number of chromosomes, number of SNPs ##################

cp $FILE.map data.map
cp $FILE.ped data.ped

tr '\t' ' ' < data.map > KK1
cut -d ' ' -f1 < KK1 > KK2

grep -w "" -c data.ped > NIND

tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1 > NCHR

SAM=$(grep -w "" -c $FILE.ped)

NCHR=$(tail -n 1 data.map | tr '\t' ' ' | cut -d ' ' -f1)

for((i=1;i<=$NCHR;i++))
do
grep -wc "$i" < KK2 > NCHR$i
done

if [ -f "SNP_CHROM" ]
then
rm SNP_CHROM
fi

for((i=1;i<=$NCHR;i++))
do
cat NCHR$i >> SNP_CHROM
done

rm KK*

################### Divide ped and map files into chromosomes ##################

echo "DIVIDE .ped AND .map FILES IN CHROMOSOMES"
echo "DIVIDE .ped AND .map FILES IN CHROMOSOMES" > timefile

num=$RANDOM
echo "$num" > seedfile

./PROGRAMMES/MANAGE_CHROMOSOMES2>>out<<@
-99
$maxNSNP
@

rm NCHR*
rm NIND
rm SNP_CHROM
###mv checkfile TEMPORARY_FILES/

################### LOOP CHROMOSOMES ##################
### Analysis of linkage disequilibrium in windows of genetic
### distances between pairs of SNPs for each chromosome

if [ $maxNCHROM != -99 ]
then
NCHR=$maxNCHROM
fi

echo "RUNNING ANALYSIS OF CHROMOSOMES ..."
echo "RUNNING ANALYSIS OF CHROMOSOMES" >> timefile

options_for_LD="$SAM $MAF $PHASE $NGEN $NBIN $ZERO $DIST $cMMb"

if [ $threads -eq -99 ]
then
threads=$(getconf _NPROCESSORS_ONLN)
fi

START=$(date +%s)

cp chromosome* TEMPORARY_FILES/

###### LD_SNP_REAL3 #######

### Obtains values of c, d2, etc. for pairs of SNPs in bins for each chromosome
for ((n=1; n<=$NCHR; n++)); do echo $n; done | xargs -I % -P $threads bash -c "./PROGRAMMES/LD_SNP_REAL3 % $options_for_LD"

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "CHROMOSOME ANALYSES took $DIFF seconds"
echo "CHROMOSOME ANALYSES took $DIFF seconds" >> timefile

######################## SUMM_REP_CHROM3 #########################
### Combination of all data gathered from chromosomes into a single output file

### Adds results from all chromosomes

for ((n=1; n<=$NCHR; n++))
do
cat outfileLD$n >> CHROM
echo "CHROMOSOME $n" >> OUTPUT
sed '2,3d' outfileLD$n > temp
mv temp outfileLD$n
cat parameters$n >> OUTPUT
done

mv outfileLD* TEMPORARY_FILES/
rm parameters*

./PROGRAMMES/SUMM_REP_CHROM3>>out<<@
$NGEN	NGEN
$NBIN	NBIN
$NCHR	NCHR
@

mv chrom* TEMPORARY_FILES/

echo "TOTAL NUMBER OF SNPs" >> OUTPUT_$FILE
cat nsnp >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
echo "HARDY-WEINBERG DEVIATION" >> OUTPUT_$FILE
cat outfileHWD >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
cat OUTPUT >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
echo "INPUT FOR GONE" >> OUTPUT_$FILE
echo -e "\n" >> OUTPUT_$FILE
cat outfileLD >> OUTPUT_$FILE

rm nsnp
rm OUTPUT
rm CHROM

############################# GONE.cpp ##########################
### Obtain estimates of temporal Ne from GONE

echo "Running GONE"
echo "Running GONE" >> timefile
START=$(date +%s)

./PROGRAMMES/GONEparallel.sh -hc $hc outfileLD $REPS

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "GONE run took $DIFF seconds"
echo "GONE run took $DIFF seconds" >> timefile

echo "END OF ANALYSES"
echo "END OF ANALYSES" >> timefile

mv outfileLD_Ne_estimates Output_Ne_$FILE
mv outfileLD_d2_sample Output_d2_$FILE
rm outfileLD
rm data.ped
rm data.map
rm out
mv outfileLD_TEMP TEMPORARY_FILES/

###################################################################
```
