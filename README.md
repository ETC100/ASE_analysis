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
