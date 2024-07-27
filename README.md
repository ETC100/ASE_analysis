# ASE_analysis
The allelic specific expression for Muscovy ducks  

## What is ASE (allele specific expression) ?
Allelic specific expression (ASE) refers to the phenomenon where alleles of a gene are expressed at different levels in a diploid organism. In a diploid organism, each individual inherits two alleles (one from each parent) for each gene. ASE occurs when one allele is expressed at a higher level than the other allele, leading to allele-specific differences in gene expression.  

ASE can be influenced by various factors, including genetic variations such as single nucleotide polymorphisms (SNPs), epigenetic modifications, transcription factor binding, and chromatin structure. It is an important aspect of gene regulation and can contribute to phenotypic diversity and disease susceptibility.  

ASE can be studied using techniques such as RNA sequencing (RNA-seq) that can differentiate between alleles based on their sequence variations. By analyzing the expression levels of alleles within an individual, researchers can uncover allele-specific patterns of gene expression and gain insights into the regulatory mechanisms underlying gene expression variability.  

## How to identify ASE?
I recommend this reference: Perspectives on Allele-Specific Expression, it provides detailed information for causes and identification of ASE.
The first step is to find suitable method to obtain paternal-maternal origin read counts table. Here is my sharing of experiences.
1. Use hisat2/star directly
When species with distant genetic distances are used for hybridization, this method can be utilized. For example, when studying mules (hybrids of horses and donkeys), one can directly merge the genomes of both species and align them using an RNA-seq aligner. Subsequently, one can remove suboptimal alignments and proceed with counting.
2. Split bam according to SNP
This method is currently the most widely used method in situations involving the same species but different strains. For instance, when studying mice, one can first align the RNA-seq data using an RNA-seq aligner to obtain a BAM file, and then use strain-specific SNP information to phase the BAM file.


The second step is to find a viable statistical method for ASE identification.
The main methods includes:
1. Binomial distribution test
Here, we have totally 6 biological replicates. Each sample have maternal lineage for A group and paternal lineage for B group.   
```python
from statsmodels.stats.proportion import binom_test
from statsmodels.stats.proportion import binom_test_reject_interval
gene_id = line_list[0]
    gene_name = line_list[1]
    chr_num = line_list[2]
    MAT1, PAT1, MAT2, PAT2, MAT3, PAT3, MAT4, PAT4, MAT5, PAT5, MAT6, PAT6 = line_list[3], line_list[4], line_list[5], line_list[6], \
                                                                              line_list[7], line_list[8], line_list[9], line_list[10], \
                                                                              line_list[11], line_list[12], line_list[13], line_list[14]
    data = [MAT1, PAT1, MAT2, PAT2, MAT3, PAT3, MAT4, PAT4, MAT5, PAT5, MAT6, PAT6] ## data[0:5] is A group, and data[6:14] is B group.

    p_value_list = []
    for i in range(6):
        p_value = binom_test(data[2*i], data[2*i] + data[2*i + 1], prop=0.5, alternative='two-sided') ## calculate the p-value one by one
        p_value_list.append(p_value)
    print(gene_id, gene_name, chr_num, sep='\t', end='\t')
    for i in range(6):
        print(data[2*i] / (data[2*i] + data[2*i + 1]), p_value_list[i], sep='\t', end='\t') ## output the maternal/paternal reads ratio and p-value.
    print()
```
Then, you can consider if each gene is ASE gene according to the p-value and ratio.

2. Bayes test
  
3. Proportional distribution test
This method has been used in Qllelic for ASE identification based on techinical replicates. You can refer to the following link. Noted that, Qllelic can also be used for biological replicates, however, the author dont recommend to do so.
https://github.com/gimelbrantlab/Qllelic
```python
from statsmodels.stats.proportion import proportions_ztest
from statsmodels.stats.proportion import proportion_confint
def DeMA_judge2(line_list):
    gene_id = line_list[0]
    gene_name = line_list[1]
    chr_num = line_list[2]
    PAT1, MAT1, PAT2, MAT2, PAT3, MAT3, PAT4, MAT4, PAT5, MAT5, PAT6, MAT6 = line_list[3], line_list[4], line_list[5], line_list[6], \
                                                                              line_list[7], line_list[8], line_list[9], line_list[10], \
                                                                              line_list[11], line_list[12], line_list[13], line_list[14]
    data = [PAT1, PAT2, PAT3, PAT4, PAT5, PAT6, MAT1, MAT2, MAT3, MAT4, MAT5, MAT6]
    try:
        int_list = [int(num) for num in data]
        pats = sum(int_list[0: 5])
        mats = sum(int_list[6: 11])
        zstat1, p_value1 = proportions_ztest(pats, pats + mats, value=0.5, alternative='two-sided', prop_var=False)
        p = pats / (pats + mats)
        se = np.sqrt(p * (1 - p) / (pats + mats))
        # zstat1, p_value1 = proportion_confint(pats, pats + mats, alpha=0.05, method='normal')
        if (p - zstat1 * se - 0.5) * (p + zstat1 * se - 0.5) > 0:
            print(gene_id, gene_name, chr_num, p, zstat1, se, p - zstat1 * se, p + zstat1 * se, p_value1, sep='\t') ## output the p-value and confidence interval based on Z test.
    except ZeroDivisionError:
        pass
```
In fact, I found the p-values calculated by propotional test were always small, while (upper_interval - 0.5)(lower_interval - 0.5) < 0 established for most genes. Their confidence interval crossed 0.5 and means that we cannot determine the propotion of one lineage is more or less than 0.5.
However. when I used binom confidence interval instead, most genes can pass the confidence interval check. Undoubtly, This is unbelievable.

# What is Allelic Methylation Regions (AMR)
Please let me to copy the answer for chatGPT, in fact, I known it too much to define it again.

Allelic methylation refers to the methylation of one allele (copy) of a gene or genomic region while the other allele remains unmethylated. This phenomenon can result in allele-specific DNA methylation patterns, where one allele is methylated and the other is not.
Allelic methylation regions are specific genomic regions where this allele-specific methylation occurs. These regions can play a role in gene regulation and expression, as methylation patterns can influence gene activity.
Studying allelic methylation regions can provide insights into epigenetic regulation and gene expression variability between alleles. Researchers use various techniques such as bisulfite sequencing and methylation-specific PCR to identify and analyze allelic methylation regions in the genome.

# Identification of AMR
Unfortunately, the identification of AMR is not as easy as ASE, because there are base change during library construction. Base C become U,  and amplified and identified as base T at last.  
So, aligner or SNP caller based on WGBS data (we assume you don't use RRBS) probably have higher demand on error checking.
We recommend Bismark for alignling, completing the big jump from fastq to bam.
https://github.com/FelixKrueger/Bismark
```
conda install bismark ## install bismark
bismark_genome_preparation --path_to_aligner /usr/bin/bowtie2/ --verbose /data/genomes/homo_sapiens/GRCh37/ # Transfer the base C to T in your reference genome
bismark [options] --genome <genome_folder> {-1 <mates1> -2 <mates2> | <singles>} # align the data to transferred genome and output bam, paired mode recommended.
```
As a result, you will get mapped and unmapped fastq as output. You'd better map the unmapped fastq to the reference again using single mode, and use samtools merge bam from paired and single modes.  
Okay, let's split the bam into maternal and paternal bam. I recommand SNPsplit software.  
https://github.com/FelixKrueger/SNPsplit?tab=readme-ov-file  
We need a bam and a vcf. The question is: How to get and rearrange the VCF?  
1.For model animals such as mouse: https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz  
2.For others: You can sequence as many as samples, and construct a breed specific SNP database yourself.  
3.There is also someone get the haplotype genomes for each sample and mapped reads on them to identify the parental SNP. That's okay, but too costly.
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.13+htslib-1.13
##bcftoolsCommand=mpileup -f Mus_musculus.GRCm39.dna.toplevel.fa.gz -b samples -g 10 -a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD -E -Q 0 -pm 3 -F 0.25 -d 500
##reference=Mus_musculus.GRCm39.dna.toplevel.fa.gz
##contig=<ID=1,length=195154279>
##contig=<ID=2,length=181755017>
##contig=<ID=3,length=159745316>
##contig=<ID=4,length=156860686>
##contig=<ID=5,length=151758149>
##contig=<ID=6,length=149588044>
##contig=<ID=7,length=144995196>
##contig=<ID=8,length=130127694>
##contig=<ID=9,length=124359700>
##contig=<ID=10,length=130530862>
##contig=<ID=11,length=121973369>
##contig=<ID=12,length=120092757>
##contig=<ID=13,length=120883175>
##contig=<ID=14,length=125139656>
##contig=<ID=15,length=104073951>
##contig=<ID=16,length=98008968>
##contig=<ID=17,length=95294699>
##contig=<ID=18,length=90720763>
##contig=<ID=19,length=61420004>
##contig=<ID=X,length=169476592>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">
##INFO=<ID=AD,Number=R,Type=Integer,Description="Total allelic depths (high-quality bases)">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=MinDP,Number=1,Type=Integer,Description="Minimum per-sample depth in this gVCF block">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled Genotype Quality">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callCommand=call -mAv -f GQ,GP -p 0.99; Date=Wed Aug 11 21:20:03 2021
##bcftools_normCommand=norm --fasta-ref Mus_musculus.GRCm39.dna.toplevel.fa.gz -m +indels; Date=Fri Aug 13 11:11:49 2021
##FORMAT=<ID=FI,Number=1,Type=Integer,Description="High confidence (1) or low confidence (0) based on soft filtering values">
##FILTER=<ID=LowQual,Description="Low quality variants">
##VEP="v104" time="2021-08-30 23:27:00" cache="mus_musculus/104_GRCm39" ensembl-funcgen=104.59ae779 ensembl-variation=104.6154f8b ensembl=104.1af1dce ensembl-io=104.1d3bb6e assembly="GRCm39" dbSNP="150" gencode="GENCODE M27" regbuild="1
.0" sift="sift"
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Am
ino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|SIFT|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">
##bcftools_viewVersion=1.13+htslib-1.13
##bcftools_viewCommand=view -i 'FORMAT/FI[*] = 1' mgp_REL2021_snps.vcf.gz; Date=Sat Dec 18 19:08:09 2021
##bcftools_annotateVersion=1.13+htslib-1.13
##bcftools_annotateCommand=annotate -x INFO/VDB,INFO/SGB,INFO/RPBZ,INFO/MQBZ,INFO/MQBZ,INFO/MQSBZ,INFO/BQBZ,INFO/SCBZ,INFO/FS,INFO/MQOF,INFO/AC,INFO/AN,FORMAT/SP,FORMAT/ADF,FORMAT/ADR,FORMAT/GP; Date=Sat Dec 18 19:08:09 2021
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##bcftools_viewCommand=view -a -Oz -o final_mgp_REL2021_snps.vcf.gz; Date=Sat Dec 18 19:08:09 2021
##bcftools_annotateCommand=annotate -x INFO/MQ0F -Oz -o final_mgp_REL2021_snps.vcf.gz mgp_REL2021_snps.vcf.gz; Date=Mon Dec 20 07:12:23 2021
#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  129P2_OlaHsd    129S1_SvImJ 129S5SvEvBrd    A_J AKR_J   B10.RIII    BALB_cByJ   BALB_cJ BTBR_T+_Itpr3tf_J   BUB_BnJ C3H_HeH C3H_HeJ C57BL_10J   C57BL_10SnJ C57BL_6NJ   C57BR_cdJ   C57L_J  C58_J   CAST_EiJ    CBA_J   CE_J    CZECHII_EiJ DBA_1J  DBA_2J  FVB_NJ  I_LnJ   JF1_MsJ KK_HiJ  LEWES_EiJ   LG_J    LP_J    MAMy_J  MOLF_EiJ    NOD_ShiLtJ  NON_LtJ NZB_B1NJ    NZO_HlLtJ   NZW_LacJ    PL_J    PWK_PhJ QSi3    QSi5    RF_J    RIIIS_J SEA_GnJ SJL_J   SM_J    SPRET_EiJ   ST_bJ   SWR_J   WSB_EiJ ZALENDE_EiJ
1   3050050 .   C   G   186.728 PASS    DP=3533;AD=3488,42;DP4=3041,447,44,1;MQ=60;CSQ=G|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||,A|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||;AC=2;AN=104  GT:PL:DP:AD:GQ:FI   0/0:0,57,168:19:19,0:71:1   0/0:0,108,255:36:36,0:122:1 0/0:0,30,154:10:10,0:44:1   0/0:0,199,255:66:66,0:127:1 0/0:0,84,222:28:28,0:98:1   0/0:0,208,255:69:69,0:127:1 0/0:0,163,255:54:54,
0:127:1 0/0:0,223,255:74:74,0:127:1 0/0:0,138,255:46:46,0:127:1 0/0:0,93,254:31:31,0:107:1  0/0:0,9,85:3:3,0:23:0   0/0:0,111,242:37:37,0:125:1 0/0:0,120,255:40:40,0:127:1 0/0:0,154,233:51:51,0:127:1 0/0:0,223,25
5:74:74,0:127:1 0/0:0,255,255:171:171,0:127:0   0/0:0,255,255:213:212,0:127:0   0/0:0,181,255:60:60,0:127:1 0/0:0,163,255:54:54,0:127:1 0/0:0,96,255:32:32,0:110:1  0/0:0,190,255:63:63,0:127:1 0/0:0,255,255:500:499,0:127:
0   0/0:0,93,255:31:31,0:107:1  0/0:0,96,255:32:32,0:110:1  0/0:0,51,233:17:17,0:65:1   0/0:0,90,255:30:30,0:104:1  0/0:0,255,255:220:220,0:127:0   0/0:0,255,255:106:106,0:127:1   0/0:0,196,255:65:65,0:127:0 0/0:
0,18,178:6:6,0:32:1 0/0:0,93,231:31:31,0:107:1  0/0:0,69,255:23:23,0:83:1   0/0:0,255,255:208:208,0:127:0   0/0:0,96,255:32:32,0:110:1  0/0:0,81,255:27:27,0:95:1   0/0:0,111,255:37:37,0:125:1 0/0:0,199,255:66:66,
0:127:1 0/0:0,144,255:49:48,0:127:1 0/0:0,57,255:19:19,0:71:1   0/0:0,255,255:366:366,0:127:0   0/0:0,93,255:31:31,0:107:1  0/0:0,78,255:26:26,0:92:1   0/0:0,69,255:23:23,0:83:1   0/0:0,187,255:62:62,0:127:1 0/0:
0,166,255:55:55,0:127:1 0/0:0,54,224:18:18,0:68:1   0/0:0,12,122:4:4,0:26:0 1/1:252,126,0:42:0,42:105:1 0/0:0,151,255:50:50,0:127:1 0/0:0,114,255:38:38,0:127:1 0/0:0,226,255:75:75,0:127:1 0/0:0,39,201:13:13,0:53:1
1   3050069 .   C   T   5422.31 PASS    DP=3674;AD=2868,796;DP4=2447,421,659,147;MQ=60;CSQ=T|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||,A|intergenic_variant|MODIFIER|||||||||||||||||||SNV||||||||,G|interge
nic_variant|MODIFIER|||||||||||||||||||SNV||||||||;AC=45;AN=104 GT:PL:DP:AD:GQ:FI   1/1:198,72,0:24:0,24:67:1   1/1:255,120,0:40:0,40:115:1 1/1:149,27,0:9:0,9:22:1 0/0:0,208,255:69:69,0:127:1 1/1:253,90,0:35:1,34:85:1   0/0:0,226,255:75:75,0:127:1 0/0:0,166,255:55:55,0:127:1 0/0:0,199,255:66:66,0:127:1 1/1:255,160,0:53:0,53:127:1 1/1:255,96,0:32:0,32:91:1   0/0:0,12,119:4:4,0:10:0 0/0:0,129,255:43:43,0:127:1 0/0:0,123,255:41:41,
0:121:1 0/0:0,148,255:50:49,0:127:1 0/0:0,223,255:74:74,0:127:1 0/1:200,0,255:175:138,37:127:0  0/1:104,0,255:235:190,44:105:0  0/0:0,163,255:54:54,0:127:1 0/1:179,0,226:55:28,27:127:0    0/0:0,93,255:31:31,0:91:1   0/1:
255,0,255:65:35,30:127:0    0/0:0,255,255:493:467,26:127:0  0/0:0,96,255:32:32,0:94:1   0/0:0,105,255:35:35,0:103:1 1/1:234,72,0:24:0,24:67:1   1/1:255,90,0:30:0,30:85:1   0/0:0,255,255:228:219,8:127:0   0/0:0,255,25
5:109:109,0:127:1   0/0:0,181,255:61:60,0:127:0 0/0:0,24,170:9:8,0:22:1 1/1:239,96,0:32:0,32:91:1   1/1:252,63,0:22:0,21:58:1   0/0:0,255,255:228:226,1:127:0   1/1:255,105,0:35:0,35:100:1 1/1:255,84,0:28:0,28:79:1   0/1:255,0,178:36:17,19:127:0    0/0:0,208,255:70:69,0:127:1 0/0:0,138,255:46:46,0:127:1 1/1:255,75,0:25:0,25:70:1   0/0:0,255,255:387:368,18:127:0  1/1:255,123,0:41:0,41:118:1 1/1:255,87,0:29:0,29:82:1   0/0:0,81,255
:27:27,0:79:1   0/0:0,205,255:68:68,0:127:1 0/0:0,172,255:57:57,0:127:1 1/1:222,57,0:19:0,19:52:1   1/1:142,18,0:6:0,6:13:0 0/0:0,114,255:38:38,0:112:1 1/1:255,148,0:49:0,49:127:1 1/1:255,117,0:39:0,39:112:1 0/0:
0,208,255:70:69,0:127:1 1/1:251,48,0:16:0,16:43:1
```
VCF is always like this, Let's rearrange it as SNPsplit input.
```python
## use SNPsplit_genome_preparation to masked the reference genome and replace the SNP on it for situation 1 and 2. AS for situation 3, no need.
## filter out 1/1 and 0/0 firstly, they have no use for phasing.
## rearrange.py, let's output the file as stdout directly.
import sys

data_vcf = "path/to/data.vcf"
print('ID', 'Chr', 'Position', 'Strand', 'Ref/SNP', sep='\t')
file = open(data_vcf)
line = file.readable()
count = 1
while line:
    line = file.readline()
    if line.startswith('#'): continue
    if len(line) == 0: break
    temp = line.split('\t')
    chrom, pos, ref, alt = temp[0], temp[1], temp[3], temp[4]
    print('snp' + str(count), chrom, pos, 1, ref+'/'+alt, sep='\t')
    count += 1
## python rearrange.py > data_rearrange.txt
```
I don't really recommond you to output the result to stdout. This is just an example.
```
## run SNPsplit
SNPsplit --bisulfite --snp_file data_rearrange.txt your_bam_file.bam -o path/to/output
## noticed that, don't add --paired or --single_end args, because your bam have both single and paired mapping results.
```
What's more, we can call SNP using WGBS data.
Bis-SNP: https://people.csail.mit.edu/dnaase/bissnp2011/  
BS-SNPer: https://github.com/hellbelly/BS-Snper  
When you only have breed specific SNP data with no NGS resequencing, and you want to get more SNP for phasing, then, WGBS-based SNP calling is used.
Bis-SNP is quite hard to install, because the authors were no longer to maintain it. Use conda to install Bis-SNP firstly, and degrade the version of gatk to Bis-SNP required, nor Bis-SNP will report an error.
There is also a pipeline for WGBS reads phasing and counting: https://github.com/BRL-BCM/allelic_epigenome.  
This process is not strict with the reads phasing. As long as the corresponding position corresponds to the read information, it asign the read to maternal/paternal lineages, allowing reads with errors in other sites. This is different from SNPsplit. When you are in above situation 3, if the quality of your haploid genome is too poor, SNPsplit may not provide any useful output.
